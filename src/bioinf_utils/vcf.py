import tempfile
import os
import gzip

import allel

from bioinf_utils import subprocess, sk_allel
from bioinf_utils.sk_allel import vcf_is_empty


def filter_vcf_region(vcf_path, chrom, start, end, out_path):
    cmd = [
        "tabix",
        "-h",
        str(vcf_path),
        f"{chrom}:{start}-{end}",
        "|",
        "bgzip",
        ">",
        str(out_path),
    ]
    subprocess.run_in_sh_script(cmd)


def filter_var_by_obs_het(vcf_path, max_allowed_var_obs_het, out_path, tmp_dir=None):
    h5 = sk_allel.read_vcf(str(vcf_path))
    genotypes = allel.GenotypeArray(h5[sk_allel.GT_FIELD])
    obs_het_per_var = sk_allel.calc_obs_het_per_var(genotypes)
    mask = obs_het_per_var > max_allowed_var_obs_het

    chroms = h5[sk_allel.CHROM_FIELD][mask]

    poss = h5[sk_allel.POS_FIELD][mask]

    if chroms.size:
        with tempfile.NamedTemporaryFile(
            suffix=".positions.txt", mode="wt", dir=tmp_dir
        ) as positions_fhand:
            for chrom, pos in zip(chroms, poss):
                positions_fhand.write(f"{chrom}\t{pos}\n")

            positions_fhand.flush()

            cmd = [
                "vcftools",
                "--recode",
                "--gzvcf",
                str(vcf_path),
                "--exclude-positions",
                positions_fhand.name,
                "--stdout",
                "|",
                "bgzip",
                "--stdout",
                ">",
                str(out_path),
            ]
            subprocess.run_in_sh_script(cmd)
    else:
        if os.path.exists(out_path):
            os.remove(out_path)
        os.symlink(vcf_path, out_path)


def filter_with_vcftools(
    input_vcf,
    out_path,
    samples_to_keep=None,
    chroms_and_poss_to_exclude: tuple[list[str], list[int]] = None,
    regions_to_include: list[tuple[str, int, int]] = None,
    min_mac=0,
    tmp_dir=None,
):
    cmd = ["vcftools", "--recode", "--gzvcf", str(input_vcf)]

    if chroms_and_poss_to_exclude:
        positions_fhand = tempfile.NamedTemporaryFile(
            suffix=".positions.txt", mode="wt", dir=tmp_dir
        )
        chroms, poss = chroms_and_poss_to_exclude
        for chrom, pos in zip(chroms, poss):
            positions_fhand.write(f"{chrom}\t{pos}\n")
        positions_fhand.flush()

        cmd.extend(
            [
                "--exclude-positions",
                positions_fhand.name,
            ]
        )
    else:
        positions_fhand = None

    if regions_to_include:
        regions_fhand = tempfile.NamedTemporaryFile(
            suffix=".regions.bed", mode="wt", dir=tmp_dir
        )
        for chrom, start, end in regions_to_include:
            regions_fhand.write(f"{chrom}\t{start}\t{end}\n")
        regions_fhand.flush()
        cmd.extend(
            [
                "--bed",
                regions_fhand.name,
            ]
        )

    if samples_to_keep:
        samples_fhand = tempfile.NamedTemporaryFile(
            suffix=".samples.txt", mode="wt", dir=tmp_dir
        )
        for sample in samples_to_keep:
            samples_fhand.write(f"{sample}\n")
        samples_fhand.flush()
        cmd.extend(["--keep", samples_fhand.name])
    else:
        samples_fhand = None

    if min_mac:
        cmd.extend(["--mac", str(min_mac)])

    cmd.extend(
        [
            "--stdout",
            "|",
            "bgzip",
            "--stdout",
            ">",
            str(out_path),
        ]
    )
    subprocess.run_in_sh_script(cmd)

    if positions_fhand:
        positions_fhand.close()
    if samples_fhand:
        samples_fhand.close()


def count_vars(vcf_path):
    return sum(1 for line in gzip.open(vcf_path, "rt") if not line.startswith("#"))
