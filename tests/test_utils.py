import itertools
import random
from functools import reduce
import operator

import numpy
import pandas


def generate_genotypes(num_vars, num_samples, ploidy, genotype_freqs):
    pass


def generate_gts_from_allelic_freqs_assuming_hw(allelic_freqs, num_samples, ploidy=2):
    num_vars, num_alleles = allelic_freqs.shape
    alleles = list(range(num_alleles))

    possible_gts = list(itertools.combinations_with_replacement(alleles, ploidy))
    gts = numpy.empty((num_vars, num_samples, ploidy), dtype=int)
    for var_idx in range(num_vars):
        allelic_freqs_this_var = allelic_freqs[var_idx, :]

        gt_freqs = []
        for gt in possible_gts:
            gt_freqs.append(
                reduce(operator.mul, [allelic_freqs_this_var[allele] for allele in gt])
            )

        for sample_idx in range(num_samples):
            this_var_this_sample_gt = random.choices(possible_gts, gt_freqs)[0]
            gts[var_idx, sample_idx, :] = this_var_this_sample_gt
    return gts


def generate_normal_phenotypes_assuming_dominance(gts, phenotype_weights_per_var):
    num_vars, num_samples, ploidy = gts.shape
    is_dominant_gt = numpy.any(gts == 0, axis=2).astype(int)
    print(is_dominant_gt)
    phenotype_weights = (
        numpy.array(phenotype_weights_per_var)
        .repeat(num_samples)
        .reshape((num_vars, num_samples))
    )
    phenotypes = is_dominant_gt * phenotype_weights
    print(phenotypes)
    phenotypes = phenotypes.sum(axis=0)
    print(phenotypes)
    # now we normalize
    normalized_phenotypes = (phenotypes - numpy.mean(phenotypes)) / numpy.std(
        phenotypes
    )
    return normalized_phenotypes


def test_gts_from_allelic_freqs():
    freqs_snp1 = [0.1, 0.9]
    freqs_snp2 = [0.5, 0.5]
    freqs_snp3 = [0.9, 0.1]
    freqs_snp4 = [0.5, 0.5]
    freqs_snp5 = [0.5, 0.5]
    allelic_freqs = numpy.array(
        [freqs_snp1, freqs_snp2, freqs_snp3, freqs_snp4, freqs_snp5]
    )
    gts = generate_gts_from_allelic_freqs_assuming_hw(allelic_freqs, num_samples=10)
    phenotype_weights_per_snp = [0.5, 0.2, 0.1, 0.05, 0]
    phenos = generate_normal_phenotypes_assuming_dominance(
        gts, phenotype_weights_per_snp
    )
