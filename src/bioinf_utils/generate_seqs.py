from collections.abc import Callable
from typing import Generator
from collections import namedtuple
import random
import math

Seq = namedtuple("Seq", ["id", "seq"])

NUCLEOTIDES = "ATCG"


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
