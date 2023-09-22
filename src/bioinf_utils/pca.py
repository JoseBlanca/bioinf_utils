import random

import numpy
import pandas
from scipy.spatial.distance import squareform

# from sklearn.preprocessing import StandardScaler
# from sklearn.decomposition import PCA

from .distances import calc_kosman_dist, calc_indi_pairwise_dists, Distances


def _center_matrix(matrix):
    means = matrix.mean(axis=0)
    return matrix - means


def _standarize_matrix(matrix):
    means = matrix.mean(axis=0)
    std_devs = matrix.std(axis=0)
    # center the matrix
    matrix = (matrix - means) / std_devs
    return matrix


def some_nan_in_numpy(array):
    return numpy.isnan(numpy.sum(array))


def do_pca(variations):
    "It does a Principal Component Analysis"

    # transform the genotype data into a 2-dimensional matrix where each cell
    # has the number of non-reference alleles per call

    matrix = variations.gts_as_mat012.T

    n_rows, n_cols = matrix.shape
    if n_cols < n_rows:
        # This restriction is in the matplotlib implementation, but I don't
        # know the reason
        msg = "The implementation requires more SNPs than samples"
        raise RuntimeError(msg)

    # Implementation based on the matplotlib PCA class
    cen_matrix = _center_matrix(matrix)
    # The following line should be added from a example to get the correct
    # variances
    # cen_scaled_matrix = cen_matrix / math.sqrt(n_rows - 1)
    cen_scaled_matrix = cen_matrix

    singular_vals, princomps = numpy.linalg.svd(cen_scaled_matrix, full_matrices=False)[
        1:
    ]
    eig_vals = singular_vals**2
    pcnts = eig_vals / eig_vals.sum() * 100.0
    projections = numpy.dot(princomps, cen_matrix.T).T

    return {
        "projections": projections,
        "var_percentages": pcnts,
        "princomps": princomps,
    }


def _make_f_matrix(matrix):
    """It takes an E matrix and returns an F matrix

    The input is the output of make_E_matrix

    For each element in matrix subtract mean of corresponding row and
    column and add the mean of all elements in the matrix
    """
    num_rows, num_cols = matrix.shape
    # make a vector of the means for each row and column
    # column_means = (numpy.add.reduce(E_matrix) / num_rows)
    column_means = (numpy.add.reduce(matrix) / num_rows)[:, numpy.newaxis]
    trans_matrix = numpy.transpose(matrix)
    row_sums = numpy.add.reduce(trans_matrix)
    row_means = row_sums / num_cols
    # calculate the mean of the whole matrix
    matrix_mean = numpy.sum(row_sums) / (num_rows * num_cols)
    # adjust each element in the E matrix to make the F matrix

    matrix -= row_means
    matrix -= column_means
    matrix += matrix_mean

    return matrix


def _get_square_dist(dists, squareform_checks=True):
    if len(dists.shape) == 1:
        return squareform(dists, checks=squareform_checks)
    else:
        return dists


def do_pcoa(dists: Distances):
    samples = dists.names
    dists = dists.dist_vector

    if numpy.any(numpy.isnan(dists)):
        raise ValueError("dists array has nan values")

    "It does a Principal Coordinate Analysis on a distance matrix"
    # the code for this function is taken from pycogent metric_scaling.py
    # Principles of Multivariate analysis: A User's Perspective.
    # W.J. Krzanowski Oxford University Press, 2000. p106.

    dists = _get_square_dist(dists)

    e_matrix = (dists * dists) / -2.0
    f_matrix = _make_f_matrix(e_matrix)

    eigvals, eigvecs = numpy.linalg.eigh(f_matrix)
    eigvecs = eigvecs.transpose()
    # drop imaginary component, if we got one
    eigvals, eigvecs = eigvals.real, eigvecs.real

    # convert eigvals and eigvecs to point matrix
    # normalized eigenvectors with eigenvalues

    # get the coordinates of the n points on the jth axis of the Euclidean
    # representation as the elements of (sqrt(eigvalj))eigvecj
    # must take the absolute value of the eigvals since they can be negative
    pca_matrix = eigvecs * numpy.sqrt(abs(eigvals))[:, numpy.newaxis]

    # output
    # get order to output eigenvectors values. reports the eigvecs according
    # to their cooresponding eigvals from greatest to least
    vector_order = list(numpy.argsort(eigvals))
    vector_order.reverse()

    eigvals = eigvals[vector_order]

    # eigenvalues
    pcnts = (eigvals / numpy.sum(eigvals)) * 100.0

    # the outputs
    # eigenvectors in the original pycogent implementation, here we name them
    # princoords
    # I think that we're doing: if the eigenvectors are written as columns,
    # the rows of the resulting table are the coordinates of the objects in
    # PCO space
    projections = []
    for name_i in range(dists.shape[0]):
        eigvect = [pca_matrix[vec_i, name_i] for vec_i in vector_order]
        projections.append(eigvect)
    projections = numpy.array(projections)
    dimensions = [f"dim-{idx}" for idx in range(1, projections.shape[1] + 1)]
    projections = pandas.DataFrame(projections, index=samples, columns=dimensions)

    pcnts = pandas.Series(pcnts, index=dimensions)

    return {"projections": projections, "var_percentages": pcnts}


def do_pcoa_from_dists(dists, samples, max_dims=5):
    multivar_result = do_pcoa(dists)

    multivar_result["projections"] = pandas.DataFrame(
        multivar_result["projections"][:, :max_dims], index=samples
    )

    return multivar_result


def do_pca(dframe, num_dims=None, standarize=True):
    if standarize:
        standarized_data = StandardScaler().fit_transform(dframe)

    pca = PCA(n_components=num_dims)
    principal_components = pca.fit_transform(standarized_data)
    result = {}
    old_dimensions = list(dframe.columns)
    new_dimensions = [f"dim_{idx + 1}" for idx in range(principal_components.shape[1])]
    result["projections"] = pandas.DataFrame(
        principal_components, index=list(dframe.index), columns=new_dimensions
    )
    result["var_percentages"] = pandas.Series(
        pca.explained_variance_ratio_ * 100, index=new_dimensions
    )
    result["princomps"] = pandas.DataFrame(
        pca.components_, index=new_dimensions, columns=old_dimensions
    )
    return result


def _select_seed_samples_for_embedding(
    vars, num_initial_samples, max_num_seed_expansions
):
    samples = numpy.array(vars.samples)
    num_samples = samples.size
    if not num_initial_samples:
        num_initial_samples = int(round(math.log2(num_samples) ** 2))
    seed_samples = numpy.array(random.sample(list(samples), k=num_initial_samples))

    cached_dists = None
    for _ in range(max_num_seed_expansions):
        seed_dists, cached_dists = _get_dists(
            vars,
            pairwise_dist_funct=calc_kosman_dist,
            pop1_samples=seed_samples,
            pop2_samples=vars.samples,
            cached_dists=None,
        )

        sample_idxs_with_max_dists_to_seeds = numpy.argmax(seed_dists.values, axis=1)
        most_distant_samples = numpy.unique(
            samples[sample_idxs_with_max_dists_to_seeds]
        )
        # print(most_distant_samples)
        dists_to_most_distant_samples, cached_dists = _get_dists(
            vars,
            pairwise_dist_funct=calc_kosman_dist,
            pop1_samples=most_distant_samples,
            pop2_samples=vars.samples,
            cached_dists=cached_dists,
        )
        samples_idxs_most_distant_to_most_distant_samples = numpy.argmax(
            dists_to_most_distant_samples.values, axis=1
        )
        samples_most_distant_to_most_distant_samples = numpy.unique(
            samples[samples_idxs_most_distant_to_most_distant_samples]
        )
        # print(samples_most_distant_to_most_distant_samples)
        old_num_seeds = seed_samples.size
        seed_samples = numpy.union1d(
            seed_samples, samples_most_distant_to_most_distant_samples
        )
        new_num_seeds = seed_samples.size
        if old_num_seeds == new_num_seeds:
            break
    return seed_samples, cached_dists


def _get_dists(
    vars, pairwise_dist_funct, pop1_samples, pop2_samples, cached_dists=None
):
    if cached_dists is None:
        samples_to_calc_dists_from = pop1_samples
    else:
        assert all(numpy.equal(pop2_samples, cached_dists.columns))
        samples_to_calc_dists_from = pop1_samples[
            numpy.logical_not(numpy.in1d(pop1_samples, cached_dists.index))
        ]

    if samples_to_calc_dists_from.size:
        new_dists = calc_indi_pairwise_dists(
            vars,
            pairwise_dist_funct=pairwise_dist_funct,
            pop1_samples=samples_to_calc_dists_from,
            pop2_samples=pop2_samples,
        )
    else:
        new_dists = None

    if cached_dists is None:
        cached_dists = new_dists
        dists = new_dists
    else:
        if new_dists is not None:
            cached_dists = pandas.concat([new_dists, cached_dists], axis="index")
        dists = cached_dists.reindex(
            index=pandas.Index(pop1_samples), columns=cached_dists.columns
        )
    return dists, cached_dists


def do_dists_and_pca_using_embeddding(
    vars, num_initial_samples=None, max_num_seed_expansions=5
):
    # following https://almob.biomedcentral.com/articles/10.1186/1748-7188-5-21
    # Sequence embedding for fast construction of guide trees for multiple sequence alignment

    seed_samples, cached_dists = _select_seed_samples_for_embedding(
        vars, num_initial_samples, max_num_seed_expansions
    )

    dists_for_embedding, _ = _get_dists(
        vars,
        pairwise_dist_funct=calc_kosman_dist,
        pop1_samples=seed_samples,
        pop2_samples=vars.samples,
        cached_dists=cached_dists,
    )

    multivar_result = do_pca(dists_for_embedding.T)

    return {"multivar_result": multivar_result, "dists": dists_for_embedding}


def do_pca_using_embeddding(vars, num_initial_samples=None, max_num_seed_expansions=5):
    # following https://almob.biomedcentral.com/articles/10.1186/1748-7188-5-21
    # Sequence embedding for fast construction of guide trees for multiple sequence alignment

    return do_dists_and_pca_using_embeddding(
        vars,
        num_initial_samples=num_initial_samples,
        max_num_seed_expansions=max_num_seed_expansions,
    )["multivar_result"]


def test_pcoa():
    dists = Distances([0.1, 0.9, 0.8, 0.7, 0.85, 0.05], names=["a", "b", "c", "d"])
    do_pcoa(dists)


if __name__ == "__main__":
    test_pcoa()
