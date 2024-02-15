import random
from abc import ABC

import Bio.Data.CodonTable
import numpy as np
from Bio import SeqIO, Phylo, AlignIO
from Bio.Align.Applications import ClustalwCommandline
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


# ======================================================================================================================
# Utils
# ======================================================================================================================


def concat_strs(strs: list[str]) -> str:
    return "".join(strs)


def sep(*args, **kwargs):
    print("\n", "=" * 120, "\n", sep="")
    print(*args, **kwargs)
    print("\n", "=" * 120, "\n", sep="")


def invert_map(m: dict[str, str]) -> dict[str, list[str]]:
    return {v: [kp for kp in m if m[kp] == v] for k, v in m.items()}


def decode_amino_to_dna(amino: str, mapper: dict[str, list[str]]) -> str:
    return random.choice(mapper[amino])


# ======================================================================================================================
# Markov
# ======================================================================================================================


class ChainType:
    def __init__(self):
        pass

    def distance(self, s1: str, s2: str) -> int:
        raise NotImplementedError

    def get_matrix(self) -> dict[str, dict[str, float]]:
        raise NotImplementedError


class JC69(ChainType, ABC):
    def __init__(self, miu: float):
        super().__init__()
        self.miu = miu

    def get_matrix(self) -> dict[str, dict[str, float]]:
        acgt = "ACGT"
        xy = 1 / len(acgt) * self.miu
        xx = 1 - (len(acgt) - 1) * xy
        return {c: {c2: (xx if c == c2 else xy) for c2 in acgt} for c in acgt}


class MarkovChain:
    @staticmethod
    def new(start: str, model: ChainType):
        return MarkovChain(start, model.get_matrix())

    def __init__(self, start: str, model: dict[str, dict[str, float]]):
        self.state = start
        self.chain = model

    def _map_single(self, c: str) -> str:
        return np.random.choice(
            list(self.chain[c].keys()), p=list(self.chain[c].values())
        )

    def step(self):
        self.state = "".join(self._map_single(c) for c in self.state)

    def simulate(self, n: int, callback: callable = None):
        for _ in range(n):
            self.step()
            if callback is not None:
                callback(self.state)

    def clone(self):
        return MarkovChain(self.state, self.chain)


# ======================================================================================================================
# Solution
# ======================================================================================================================


# For
# ```
# d1 = {'a': 1, 'b': 2}
# d2 = {'c': 3, 'd': 4}
# ```
# `join([d1, d2])` returns `{'a': 1, 'b': 2, 'c': 3, 'd': 4}`
def join(seqs: list[dict]) -> dict:
    return {k: v for d in seqs for k, v in d.items()}


def traverse_tree(node, miu: float, mc: MarkovChain) -> dict[str, str]:
    n = (node.branch_length / 100) / miu
    mc.simulate(int(n))
    children_seqs = [traverse_tree(child, miu, mc.clone()) for child in node.clades]
    if children_seqs:
        return join(children_seqs)
    else:
        return {node.name: mc.state}


def solution():
    # Zadanie 1.
    origin_seq = next(SeqIO.parse("PAH.fa", "fasta"))
    sep("Origin sequence:", str(origin_seq), str(origin_seq.seq), sep="\n")
    amino_to_dna = invert_map(Bio.Data.CodonTable.standard_dna_table.forward_table)

    dna = concat_strs([decode_amino_to_dna(el, amino_to_dna) for el in origin_seq.seq])
    dna_seq = Seq(dna)
    assert dna_seq.translate() == origin_seq.seq  # sanity check

    # Zadanie 2.
    expected_phylo_tree = Phylo.read("tree", "newick")
    sep("Expected tree:", expected_phylo_tree)

    # Zadanie 3.
    miu = 1 / 10000
    jc69 = JC69(miu)
    mc = MarkovChain.new(dna, jc69)
    generated_seqs = traverse_tree(expected_phylo_tree.clade, jc69.miu, mc.clone())

    # Zadanie 4.
    translated_seqs = {k: Seq(v).translate() for k, v in generated_seqs.items()}
    written = SeqIO.write(
        [SeqRecord(Seq(v), id=k, description=k) for k, v in translated_seqs.items()],
        "bialka.fa",
        "fasta",
    )
    assert written == len(translated_seqs)  # sanity check

    # Zadanie 5.
    ClustalwCommandline(cmd="clustalw", infile=f"bialka.fa", matrix="BLOSUM")()

    # Zadanie 6.
    alignment = AlignIO.read("bialka.aln", "clustal")
    calculator = DistanceCalculator("identity")
    distance_matrix = calculator.get_distance(alignment)
    tree_constructor = DistanceTreeConstructor()
    generated_phylo_tree_upgma = tree_constructor.upgma(distance_matrix)
    generated_phylo_tree_nj = tree_constructor.nj(distance_matrix)

    Phylo.write(generated_phylo_tree_upgma, "tree_upgma", "newick")
    Phylo.write(generated_phylo_tree_nj, "tree_nj", "newick")


if __name__ == "__main__":
    solution()
