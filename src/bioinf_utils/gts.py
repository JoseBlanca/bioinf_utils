import numpy
import pandas

import utils_sk_allel


def create_haplotypic_gts_from_gts(gts):
    different_gts = numpy.unique(gts)
    if utils_sk_allel.MISSING_GT in different_gts:
        raise ValueError("There are missing GTs in the genotype matrix")

    n_snps, n_samples, ploidy = gts.shape
    gts2d = gts.reshape(n_snps, n_samples * ploidy)
    unique_haplotypes, haplotypic_gts_1d = numpy.unique(
        gts2d, axis=1, return_inverse=True
    )
    haplotypic_gts = haplotypic_gts_1d.reshape(1, n_samples, ploidy)
    return {"haplotypes": unique_haplotypes, "gts": haplotypic_gts}


def generate_tuple_gt_dframe(
    gts, samples=None, chroms=None, poss=None
) -> pandas.DataFrame:
    gts = numpy.rec.fromarrays(
        [gts[:, :, ploidy_idx] for ploidy_idx in range(gts.shape[2])]
    )
    kwargs = {}
    if samples is not None:
        kwargs["columns"] = samples
    if chroms is not None and poss is not None:
        chroms_and_poss = list(zip(chroms, poss))
        kwargs["index"] = chroms_and_poss
    gts = pandas.DataFrame(gts.tolist(), **kwargs)

    # ensure that hets are correctly sorted
    different_gts = numpy.unique(gts)
    het_gts = [gt for gt in different_gts if len(set(gt)) > 1]
    for het_gt in het_gts:
        sorted_het_gt = tuple(sorted(het_gt))
        if sorted_het_gt != het_gt:
            gts = gts.applymap(lambda cell: sorted_het_gt if cell == het_gt else cell)
    return gts
