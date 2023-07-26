import itertools

import numpy
import pandas

from bioinf_utils.sk_allel import GT_FIELD, get_samples


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


MISSING_INT = -1
MISSING_FLOAT = float("nan")
MISSING_STR = ""
MISSING_BYTE = b""
MISSING_BOOL = False


class _MissingValues:
    def __init__(self):
        self._missing_values = {
            int: MISSING_INT,
            "Integer": MISSING_INT,
            float: MISSING_FLOAT,
            "Float": MISSING_FLOAT,
            str: MISSING_STR,
            "String": MISSING_STR,
            numpy.int8: MISSING_INT,
            numpy.int16: MISSING_INT,
            numpy.int32: MISSING_INT,
            numpy.float16: MISSING_FLOAT,
            numpy.float32: MISSING_FLOAT,
            numpy.bool_: MISSING_BOOL,
            numpy.bytes_: MISSING_BYTE,
            bool: MISSING_BOOL,
        }

    def __getitem__(self, dtype):
        str_dtype = str(dtype)
        if dtype in self._missing_values:
            return self._missing_values[dtype]
        elif isinstance(dtype, str):
            if "str" in dtype:
                return MISSING_STR
            elif "int" in dtype:
                return MISSING_INT
            elif "float" in dtype:
                return MISSING_FLOAT
            elif dtype[0] == "S":
                return MISSING_BYTE
            elif dtype[:2] == "|S":
                return MISSING_BYTE
        elif "int" in str_dtype:
            return MISSING_INT
        elif "float" in str_dtype:
            return MISSING_FLOAT
        elif "bool" in str_dtype:
            return MISSING_BOOL
        elif str_dtype[:2] == "|S":
            return MISSING_BYTE
        elif str_dtype[:2] == "<U":
            return MISSING_STR
        else:
            raise ValueError("No missing type defined for type: " + str(dtype))


MISSING_VALUES = _MissingValues()


def _is_missing(matrix, axis=1):
    if axis is None:
        return matrix == MISSING_VALUES[matrix.dtype]
    else:
        return numpy.any(matrix == MISSING_VALUES[matrix.dtype], axis=axis)


def calc_indi_pairwise_dists(
    vars_h5,
    pairwise_dist_funct,
    pop1_samples=None,
    pop2_samples=None,
    min_num_snps=None,
):
    gts = vars_h5[GT_FIELD][:]
    gts.flags.writeable = False
    is_missing_gt = _is_missing(gts, axis=2)

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
        samples = set(get_samples(vars_h5))
        pop1_sample_idxs = [
            idx for idx, sample in enumerate(samples) if sample in pop1_samples
        ]
        pop2_sample_idxs = [
            idx for idx, sample in enumerate(samples) if sample in pop2_samples
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
