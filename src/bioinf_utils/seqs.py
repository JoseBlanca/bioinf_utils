from collections.abc import Callable
from typing import Generator
from collections import namedtuple
import random
import math

Seq = namedtuple("Seq", ["id", "seq"])

NUCLEOTIDES = "ATCG"
REVERSE_COMPLEMENT = {
    "A": "T",
    "T": "A",
    "G": "C",
    "C": "G",
    "-": "-",
    "R": "Y",
    "Y": "R",
    "S": "S",
    "W": "W",
    "K": "M",
    "M": "K",
    "B": "V",
    "V": "B",
    "D": "H",
    "H": "D",
    "N": "N",
}


def uniform_seq_distrib(seq_lengh: int) -> Generator[int, None, None]:
    while True:
        yield seq_lengh


def generate_random_seqs(
    seq_length_distrib: Callable | int, gc_content: float = 0.5
) -> Generator[Seq, None, None]:
    if isinstance(seq_length_distrib, int):
        seq_length_distrib = uniform_seq_distrib(seq_length_distrib)

    if gc_content < 0.0 or gc_content > 1.0:
        raise ValueError(
            f"gc_content should be between 0 and 1, but it is: {gc_content}"
        )

    g_content = gc_content / 2.0
    a_content = (1 - gc_content) / 2.0
    weights = (a_content, a_content, g_content, g_content)
    assert math.isclose(sum(weights), 1.0)

    for idx, seq_len in enumerate(seq_length_distrib):
        seq = "".join((random.choices(NUCLEOTIDES, weights)[0] for _ in range(seq_len)))
        id = f"seq_{idx}"
        yield Seq(id, seq)


def reverse_complement(seq: Seq, id_modifier: Callable | None = None) -> Seq:
    if id_modifier is None:
        id = f"{seq.id}_rev"
    else:
        id = id_modifier(seq.id)

    seq = Seq(
        id=id, seq="".join((REVERSE_COMPLEMENT[nucl.upper()] for nucl in seq.seq[::-1]))
    )
    return seq


def add_seqs(seq1: Seq, seq2: Seq, id_modifier: Callable | None = None):
    if id_modifier is None:
        id = f"{seq1.id}+{seq2.id}"
    else:
        id = id_modifier(seq1.id, seq2.id)
    seq = Seq(id=id, seq=seq1.seq + seq2.seq)
    return seq
