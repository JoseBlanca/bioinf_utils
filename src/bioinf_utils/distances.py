import itertools
from io import StringIO

import numpy
import pandas
import scipy


class Distances:
    def __init__(self, dist_vector, names):
        if not isinstance(dist_vector, numpy.ndarray):
            dist_vector = numpy.array(dist_vector)
        self.names = numpy.array(names)
        self._dist_vector = dist_vector

    @property
    def dist_vector(self):
        return self._dist_vector

    @property
    def square_dists(self):
        dists = scipy.spatial.distance.squareform(self.dist_vector)
        dists = pandas.DataFrame(dists, index=self.names, columns=self.names)
        return dists

    @property
    def triang_list_of_lists(self):
        dist_vector = iter(self.dist_vector)
        length = 0
        dists = []
        while True:
            dist_row = list(itertools.islice(dist_vector, length))
            if length and not dist_row:
                break
            dist_row.append(0)
            dists.append(dist_row)
            length += 1
        return dists

    def write_mega(self, fhand):
        fhand.write("\n#MEGA\n!Format DataType=distance;\n\n")
        for name in self.names:
            fhand.write(f"#{name}\n")
        fhand.write("\n")

        sep = "   "
        for idx, row in enumerate(self.square_dists.values):
            row = row[:idx]
            if not row.size:
                continue
            fhand.write(sep)
            fhand.write(sep.join(map(lambda x: f"{x:.5f}", row)))
            fhand.write("\n")


def calc_dists_from_pca(pca_result, metric="euclidean"):
    projections = pca_result["projections"]
    dists = scipy.spatial.distance.pdist(projections, metric=metric)
    return Distances(names=projections.index, dist_vector=dists)


def calc_kosman_dist(indi1_gt, indi2_gt):
    """It calculates the distance between two individuals using the Kosman dist

    The Kosman distance is explain in DOI: 10.1111/j.1365-294X.2005.02416.x
    """

    indi1 = indi1_gt
    indi2 = indi2_gt

    if indi1.shape[1] != 2:
        raise ValueError("Only diploid are allowed")

    alleles_comparison1 = indi1 == indi2.transpose()[:, :, None]
    alleles_comparison2 = indi2 == indi1.transpose()[:, :, None]

    result = numpy.add(
        numpy.any(alleles_comparison2, axis=2).sum(axis=0),
        numpy.any(alleles_comparison1, axis=2).sum(axis=0),
    )

    result2 = numpy.full(result.shape, fill_value=0.5)
    result2[result == 0] = 1
    result2[result == 4] = 0
    return result2.sum(), result2.shape[0]


def _get_sample_gts(gts, sample_i, sample_j, is_missing_gt):
    indi1 = gts[:, sample_i]
    is_missing_1 = is_missing_gt[:, sample_i]

    indi2 = gts[:, sample_j]
    is_missing_2 = is_missing_gt[:, sample_j]

    is_called = numpy.logical_not(numpy.logical_or(is_missing_1, is_missing_2))

    indi1 = indi1[is_called]
    indi2 = indi2[is_called]

    assert issubclass(indi1.dtype.type, numpy.integer)
    assert issubclass(indi2.dtype.type, numpy.integer)

    return indi1, indi2


def calc_indi_pairwise_dists(
    variations,
    pairwise_dist_funct,
    pop1_samples=None,
    pop2_samples=None,
    min_num_snps=None,
):
    gts = variations[variation.GT_FIELD][:]
    gts.flags.writeable = False
    is_missing_gt = variation.matrix.methods.is_missing(gts, axis=2)

    if pop1_samples is None:
        n_samples = gts.shape[1]
        num_dists_to_calculate = int((n_samples**2 - n_samples) / 2)
        dists = numpy.zeros(num_dists_to_calculate)
        n_snps_matrix = numpy.zeros(num_dists_to_calculate)
    else:
        shape = (len(pop1_samples), len(pop2_samples))
        dists = numpy.zeros(shape)
        n_snps_matrix = numpy.zeros(shape)

    if pop1_samples is None:
        sample_combinations = itertools.combinations(range(n_samples), 2)
    else:
        pop1_sample_idxs = [
            idx
            for idx, sample in enumerate(variations.samples)
            if sample in pop1_samples
        ]
        pop2_sample_idxs = [
            idx
            for idx, sample in enumerate(variations.samples)
            if sample in pop2_samples
        ]
        sample_combinations = itertools.product(pop1_sample_idxs, pop2_sample_idxs)

    index = 0
    for sample_i, sample_j in sample_combinations:
        indi1_gt, indi2_gt = _get_sample_gts(gts, sample_i, sample_j, is_missing_gt)
        abs_dist, n_snps = pairwise_dist_funct(indi1_gt, indi2_gt)

        if min_num_snps is not None:
            n_snps[n_snps < min_num_snps] = numpy.nan
        with numpy.errstate(invalid="ignore"):
            dist = abs_dist / n_snps

        if pop1_samples is None:
            dists[index] = dist
            n_snps_matrix[index] = n_snps
            index += 1
        else:
            dists_samplei_idx = pop1_sample_idxs.index(sample_i)
            dists_samplej_idx = pop2_sample_idxs.index(sample_j)
            dists[dists_samplei_idx, dists_samplej_idx] = dist
            n_snps_matrix[dists_samplei_idx, dists_samplej_idx] = n_snps

    if dists.ndim == 2:
        dists = pandas.DataFrame(dists, index=pop1_samples, columns=pop2_samples)
    else:
        dists = pandas.Series(dists)
    return dists


def test_to_mega():
    dists = Distances(["a", "b", "c"], [0.5, 0.3, 0.2])
    fhand = StringIO()
    dists.write_mega(fhand)
    text = fhand.getvalue()
    assert "#a\n#b\n#c" in text
    assert "   0.50000\n   0.30000   0.20000\n"


if __name__ == "__main__":
    test_to_mega()
    dframe = pandas.read_csv(
        open(
            "candidate_root_genes.extended_2500pb.tier1_low_qual_085_vars.dists", "rt"
        ),
        sep="\t",
        index_col=0,
        header=0,
    )
    dists = Distances(dframe.columns, scipy.spatial.distance.squareform(dframe.values))
    fhand = open(
        "candidate_root_genes.extended_2500pb.tier1_low_qual_085_vars.meg", "wt"
    )
    print(dists.write_mega(fhand))
