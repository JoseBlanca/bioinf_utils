import itertools
from io import StringIO

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
