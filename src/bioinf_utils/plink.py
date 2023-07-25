import subprocess
from pathlib import Path
import os
from numpy import var

import pandas

PLINK_BIN = "plink"
PLINK2_BIN = "plink2"


class VariantFilters:
    def __init__(
        self,
        max_major_freq=0.01,
        max_missing_rate=0.1,
        lists_of_vars_to_keep_paths: list[Path] = None,
    ):
        self.max_major_freq = max_major_freq
        self.max_missing_rate = max_missing_rate

        if lists_of_vars_to_keep_paths is None:
            lists_of_vars_to_keep_paths = []
        self.lists_of_vars_to_keep_paths = lists_of_vars_to_keep_paths

    def create_cmd_arg_list(self):
        args = []
        if self.max_major_freq is not None:
            args.extend(["--maf", str(self.max_major_freq)])
        if self.max_missing_rate is not None:
            args.extend(["--geno", "0.1"])
        if self.lists_of_vars_to_keep_paths:
            args.append("--extract")
            args.extend(map(str, self.lists_of_vars_to_keep_paths))
        return args


def run_cmd(cmd, stdout_path, stderr_path):
    stdout_fhand = stdout_path.open("wt")
    stderr_fhand = stderr_path.open("wt")
    try:
        subprocess.run(
            cmd,
            stdout=stdout_fhand,
            stderr=stderr_fhand,
            check=True,
        )
    except subprocess.CalledProcessError:
        print("stdout")
        stdout_fhand = stdout_path.open("rt")
        print(stdout_fhand.read())
        stdout_fhand.close()
        print("stderr")
        stderr_fhand = stderr_path.open("rt")
        print(stderr_fhand.read())
        stderr_fhand.close()
        raise
    stderr_fhand.close()
    stdout_fhand.close()


def create_plink_bfiles(vcf_path, out_base_path):
    if not vcf_path.exists():
        raise ValueError(f"VCF file not found: {vcf_path}")

    cmd = [PLINK_BIN]
    cmd.extend(["--vcf", str(vcf_path)])
    cmd.extend(["--out", str(out_base_path)])
    cmd.extend(
        [
            "--allow-no-sex",
            "--make-bed",
            "--allow-extra-chr",
            "--double-id",
            "--vcf-half-call",
            "missing",
            "--set-missing-var-ids",
            "@:# ",
        ]
    )

    stderr_path = Path(str(out_base_path) + ".bfiles_creation.stderr")
    stdout_path = Path(str(out_base_path) + ".bfiles_creation.stdout")
    run_cmd(cmd, stdout_path, stderr_path)

    out_paths = {}
    dir_ = out_base_path.parent
    base_name = out_base_path.name
    for path in dir_.iterdir():
        if base_name in path.name:
            if path.suffix == ".fam":
                out_paths["fam"] = path
            if path.suffix == ".bed":
                out_paths["bed"] = path
            if path.suffix == ".bim":
                out_paths["bim"] = path
    return {
        "out_paths": out_paths,
        "bfiles_base_path": out_base_path,
        "stdout_path": stdout_path,
        "stderr_path": stderr_path,
    }


def write_covariate_file(covariates: pandas.DataFrame, covars_fhand):
    sep = " "
    items = ["FID", "IID"] + list(covariates.columns)
    covars_fhand.write(sep.join(items))
    covars_fhand.write("\n")

    for sample_id, row in covariates.iterrows():
        to_write = sep.join([sample_id, sample_id] + list(map(str, row)))
        covars_fhand.write(f"{to_write}\n")
    covars_fhand.flush()


def _read_mds_projections_file(mds_path):
    fhand = mds_path.open("rt")
    projections = []
    ids = []
    columns = None
    for line in fhand:
        if not columns:
            columns = [item.strip() for item in line.strip().split()][3:]
            continue
        items = [item.strip() for item in line.strip().split()]
        ids.append(items[0])
        projections.append([float(number) for number in items[3:]])
    projections = pandas.DataFrame(projections, index=ids, columns=columns)
    fhand.close()
    return projections


def do_mds(
    bfiles_base_path,
    out_base_path,
    n_dims=10,
    variant_filters: VariantFilters | None = None,
):
    # Identity-by-descent with --genome
    cmd = [PLINK_BIN]
    cmd.extend(["--bfile", str(bfiles_base_path)])
    cmd.extend(["--genome"])
    cmd.append("--allow-no-sex")
    cmd.append("--allow-extra-chr")
    cmd.extend(["-out", str(out_base_path)])

    stderr_path = Path(str(out_base_path) + ".ibd.stderr")
    stdout_path = Path(str(out_base_path) + ".ibd.stdout")
    run_cmd(cmd, stdout_path, stderr_path)

    # MDS
    cmd = [PLINK_BIN]
    cmd.extend(["--bfile", str(bfiles_base_path)])
    cmd.extend(["--read-genome", str(out_base_path) + ".genome"])
    cmd.append("--allow-extra-chr")
    cmd.append("--allow-no-sex")
    cmd.extend(["--cluster", "--mds-plot", str(n_dims)])
    cmd.extend(["-out", str(out_base_path)])

    if variant_filters is not None:
        cmd.extend(variant_filters.create_cmd_arg_list())

    stderr_path = Path(str(out_base_path) + ".mds.stderr")
    stdout_path = Path(str(out_base_path) + ".mds.stdout")
    run_cmd(cmd, stdout_path, stderr_path)
    mds_path = Path(str(out_base_path) + ".mds")

    projections = _read_mds_projections_file(mds_path)

    return {"mds_path": mds_path, "projections": projections}


def do_pca(
    bfiles_base_path,
    out_base_path,
    variant_filters: VariantFilters,
    n_dims=10,
    approx=False,
):
    if not variant_filters.max_major_freq:
        raise ValueError("It is really important to filter out the low freq variants")
    if variant_filters.max_missing_rate < 0.9:
        raise ValueError(
            "Consider that PCA will fill missing values with mean, so use variants with few missing data"
        )

    out_dir = out_base_path if out_base_path.is_dir() else out_base_path.parent

    cmd = [PLINK2_BIN]
    cmd.extend(["--bfile", str(bfiles_base_path)])
    cmd.append("--allow-extra-chr")
    cmd.extend(["-out", str(out_base_path)])

    cmd.append("--pca")
    cmd.append(str(n_dims))
    if approx:
        cmd.append("approx")

    if variant_filters is not None:
        cmd.extend(variant_filters.create_cmd_arg_list())

    stderr_path = Path(str(out_base_path) + ".pca.stderr")
    stdout_path = Path(str(out_base_path) + ".pca.stdout")
    run_cmd(cmd, stdout_path, stderr_path)
    eigenvec_path = Path(str(out_base_path) + ".eigenvec")
    projections = _read_mds_projections_file(eigenvec_path)

    return {"eigenvec_path": eigenvec_path, "projections": projections}


def do_gwas(
    bfiles_base_path,
    phenotypes_path,
    test_type,
    out_base_path,
    covars_path=None,
    allow_no_covars=False,
    variant_filters: VariantFilters | None = None,
):
    out_base_path = Path(str(out_base_path) + ".gwas")

    cmd = [PLINK2_BIN]
    cmd.extend(["--bfile", str(bfiles_base_path)])

    cmd.append("--allow-extra-chr")

    # Note: --allow-no-sex no longer has any effect.  (Missing-sex samples are
    # automatically excluded from association analysis when sex is a covariate, and
    # treated normally otherwise.)
    # cmd.append("--allow-no-sex")

    cmd.extend(["--adjust", "cols=+qq"])

    if test_type == "linear":
        cmd.append("--linear")
    elif test_type == "logistic":
        cmd.append("--logistic")
    else:
        raise ValueError(
            f"Unknown test type ({test_type}, it should be linear (quantitative) or logistic (qualitative))"
        )

    if covars_path:
        cmd.extend(["--covar", str(covars_path)])
    else:
        if allow_no_covars:
            cmd.append("allow-no-covars")
        else:
            raise ValueError("If allow_no_covars is False should should provide covars")

    cmd.extend(["--pheno", str(phenotypes_path)])

    if variant_filters is not None:
        cmd.extend(variant_filters.create_cmd_arg_list())

    cmd.extend(["-out", str(out_base_path)])

    stderr_path = Path(str(out_base_path) + ".gwas.stderr")
    stdout_path = Path(str(out_base_path) + ".gwas.stdout")

    run_cmd(cmd, stdout_path, stderr_path)

    stdout_fhand = stdout_path.open("rt")
    stdout = stdout_fhand.read()
    stdout_fhand.close()

    if "Zero valid tests; --adjust skipped." in stdout:
        pvalues = None
    else:
        test_str = "logistic.hybrid" if test_type == "logistic" else test_type
        adjusted_pvalues_path = Path(
            str(out_base_path) + f".PHENO1.glm.{test_str}.adjusted"
        )
        pvalues = pandas.read_csv(
            str(adjusted_pvalues_path), delim_whitespace=True, index_col="ID"
        )

    col_mapping = {
        "#CHROM": "chrom",
        "ID": "variant_id",
        "UNADJ": "pval",
        "GC": "genomic_control_corrected_pval",
        "QQ": "pval_quantile",
        "BONF": "bonferroni_pval",
        "HOLM": "holm_bonferroni_pval",
        "SIDAK_SS": "sidak_single_step_pval",
        "SIDAK_SD": "sidak_step_down_pval",
        "FDR_BH": "benjamini_hochberg_pval",
        "FDR_BY": "benjamini_yekutieli_pval",
    }
    if pvalues is not None:
        pvalues.columns = [col_mapping.get(col, col) for col in pvalues.columns]
        result = {"adjusted_pvalues": pvalues}
    else:
        result = {}
    return result


def write_phenotype_file(phenotypes: dict, fhand, quantitative=False):
    for acc, phenotype in phenotypes.items():
        if quantitative:
            phenotype = int(phenotype)
            if phenotype not in (0, 1, 2, -9):
                raise ValueError(
                    f"Phenotypes should be 0, 1, 2, -9, but there is: {phenotype} for acc {acc}"
                )
        fhand.write(f"{acc}\t{acc}\t{phenotype}\n")
    fhand.flush()


def create_ld_indep_variants_file(
    bfiles_base_path,
    out_base_path,
    pruned_vars_list_path,
    variant_filters: VariantFilters,
    window_size=50,
    step_size=5,
    sizes_are_in_number_of_vars=True,
    r2_threshold=0.5,
):
    pruned_vars_list_path = Path(pruned_vars_list_path)

    if not variant_filters.max_major_freq:
        raise ValueError("It is really important to filter out the low freq variants")

    current_working_dir = os.getcwd()
    working_dir = out_base_path if out_base_path.is_dir() else out_base_path.parent
    os.chdir(working_dir)

    cmd = [PLINK2_BIN]
    cmd.extend(["--bfile", str(bfiles_base_path)])

    cmd.append("--allow-extra-chr")
    # Note: --allow-no-sex no longer has any effect.  (Missing-sex samples are
    # automatically excluded from association analysis when sex is a covariate, and
    # treated normally otherwise.)
    # cmd.append("--allow-no-sex")

    if variant_filters is not None:
        cmd.extend(variant_filters.create_cmd_arg_list())

    cmd.append("--indep-pairwise")
    if sizes_are_in_number_of_vars:
        cmd.extend([str(window_size), str(step_size), str(r2_threshold)])
    else:
        cmd.extend([f"{window_size}kb", str(step_size), str(r2_threshold)])

    out_base_path = Path(str(out_base_path) + ".ld")

    stderr_path = Path(str(out_base_path) + ".stderr")
    stdout_path = Path(str(out_base_path) + ".stdout")

    run_cmd(cmd, stdout_path, stderr_path)

    plink_out_path = working_dir / "plink2.prune.in"
    if not pruned_vars_list_path.exists() or not os.path.samefile(
        plink_out_path, pruned_vars_list_path
    ):
        os.rename(plink_out_path, pruned_vars_list_path)

    os.chdir(current_working_dir)
