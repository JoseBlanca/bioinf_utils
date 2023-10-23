from itertools import islice
import random

import pytest

from bioinf_utils.seqs import generate_random_seqs, Seq, reverse_complement


def test_generate_seqs():
    random.seed(42)
    seqs = list(islice(generate_random_seqs(seq_length_distrib=10), 2))
    assert seqs[0].seq == "CATACCGATA"
    assert seqs[1].seq == "ACAACCACGA"

    with pytest.raises(ValueError):
        list(generate_random_seqs(seq_length_distrib=1, gc_content=1.5))

    seqs = list(islice(generate_random_seqs(seq_length_distrib=10, gc_content=0.1), 2))
    assert seqs[0].seq == "TTAAGAAATT"
    assert seqs[1].seq == "TTTGATTTTT"


def test_reverse_complement():
    seq = Seq(id="seq1", seq="ACT-GT")
    seq = reverse_complement(seq)
    assert seq.seq == "AC-AGT"
