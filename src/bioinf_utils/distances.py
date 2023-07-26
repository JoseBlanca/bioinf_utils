import itertools

import numpy
import pandas
import scipy


class Distances:
    def __init__(self, names, dist_vector):
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


def calc_dists_from_pca(pca_result, metric="euclidean"):
    projections = pca_result["projections"]
    dists = scipy.spatial.distance.pdist(projections, metric=metric)
    return Distances(names=projections.index, dist_vector=dists)
