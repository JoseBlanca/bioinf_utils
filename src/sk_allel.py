import config

from pathlib import Path

import h5py
import allel
import numpy
import pandas

import utils_vcf

GT_FIELD = "calldata/GT"
CHROM_FIELD = "variants/CHROM"
POS_FIELD = "variants/POS"
SAMPLE_FIELD = "samples"
MISSING_GT = -1


class EmptyVCFError(ValueError):
    pass


def read_vcf(vcf_path: Path, fields=None):
    try:
        data = allel.read_vcf(str(vcf_path), fields=fields)
        if data is None:
            if utils_vcf.vcf_is_empty(vcf_path):
                raise EmptyVCFError(f"Empty VCF: {vcf_path}")
    except RuntimeError:
        if utils_vcf.vcf_is_empty(vcf_path):
            raise EmptyVCFError(f"Empty VCF: {vcf_path}")
        raise RuntimeError(f"Error reading VCF: {vcf_path}")
    return data


def get_h5(h5_path: Path | h5py.File | dict):
    if isinstance(h5_path, dict):
        return h5_path
    if isinstance(h5_path, h5py.File):
        return h5_path
    h5 = h5py.File(h5_path, "r")
    return h5


def load_genotypes(h5: Path | h5py.File):
    h5 = get_h5(h5)
    genotypes = allel.GenotypeArray(h5[GT_FIELD])
    return genotypes


def calc_allelic_freqs(genotypes, subpop_idxs: dict[list]):
    counts_for_pops = genotypes.count_alleles_subpops(subpop_idxs)
    freqs = {}
    for pop, counts in counts_for_pops.items():
        freqs[pop] = counts.to_frequencies()
    return freqs


def get_samples(h5):
    return [sample.decode() for sample in get_h5(h5)[SAMPLE_FIELD][:]]


def get_chrom_and_pos(h5):
    h5 = get_h5(h5)
    chroms = h5[CHROM_FIELD][:].astype("U")
    poss = h5[POS_FIELD][:]
    return chroms, poss


ANNOTATION_COLUMNS = (
    "allele",
    "predicted_effect",
    "impact",
    "gene_name",
    "gene_id",
    "feature_type",
    "feature_id",
    "transcript_type",
    "rank",
    "dna_mutation(HGVS.c)",
    "prot_mutation(HGVS.p)",
    "cDNA.pos/cDNA.length",
    "CDS.pos/CDS.length",
    "AA.pos/AA.length",
    "distance_to_feature",
    "errors",
)


def get_snpeff_info(genotyping_data):
    h5 = get_h5(genotyping_data["h5"])
    # Functional annotations: 'Allele | Annotation | Annotation_Impact | Gene_Name | Gene_ID | Feature_Type
    # | Feature_ID | Transcript_BioType | Rank | HGVS.c | HGVS.p | cDNA.pos / cDNA.length | CDS.pos / CDS.len
    result = {}
    annotation = h5["variants/ANN"][:].astype("U")
    if not annotation[0]:
        return {}

    chroms_and_poss = tuple(zip(genotyping_data["chroms"], genotyping_data["poss"]))
    annotation = pandas.DataFrame(
        list(numpy.char.split(annotation, "|")),
        columns=ANNOTATION_COLUMNS,
        index=chroms_and_poss,
    )
    result["annotation"] = annotation

    # Predicted loss of function effects for this variant. Format: 'Gene_Name | Gene_ID | Number_of_transcripts_in_gene | Percent_of_transcripts_affected'
    result["loss_of_function"] = h5["variants/LOF"][:]
    # Predicted nonsense mediated decay effects for this variant. Format: 'Gene_Name | Gene_ID | Number_of_transcripts_in_gene | Percent_of_transcripts_affected
    result["nonsense"] = h5["variants/NMD"][:]
    return result


def calc_obs_het_per_var(genotypes: allel.GenotypeArray):
    fill = numpy.nan
    n_het = numpy.asarray(genotypes.count_het(axis=1))
    n_called = numpy.asarray(genotypes.count_called(axis=1))

    with allel.util.ignore_invalid():
        obs_het = numpy.where(n_called > 0, n_het / n_called, fill)

    return obs_het


if __name__ == "__main__":
    print(get_samples(config.H5_ALLEL_CORE_SELECTION))
    res = get_snpeff_info(config.H5_ALLEL_CORE_SELECTION)
    lof = res["loss_of_function"]
    non = res["nonsense"]
    print(lof.shape)
    mask = numpy.logical_or(lof != b"", non != b"")
    lof = lof[mask]
    non = non[mask]
    print(lof)
    print(non)
    print(lof.shape)
