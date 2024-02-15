#!/usr/bin/env python3

import random

import Bio.motifs
import Bio.Seq


class MotifData:
    def __init__(self, motif, value, idx):
        self.motif = motif
        self.value = value
        self.id = idx

    def __gt__(self, other):
        return self.value > other.value

    @staticmethod
    def of(motif, value):
        return MotifData(motif, value, None)


class BestKeeper:
    def __init__(self, size):
        self.size = size
        self.data = []

    def normalize(self):
        self.data.sort(reverse=True)
        self.data = self.data[: self.size]

    def add(self, motif_data):
        self.add_all([motif_data])

    def add_all(self, motif_datas: []):
        self.data.extend(motif_datas)
        self.normalize()

    def __iadd__(self, other: []):
        self.add_all(other)
        return self

    def into_list(self):
        return self.data


def sliding_window(elements, window_size):
    for i in range(len(elements) - window_size + 1):
        yield elements[i : i + window_size]


def without_ith_element(elements, i):
    return elements[:i] + elements[i + 1 :]


def consensus_with_unified_motif(sequences, motif_length, selected_motif):
    best = MotifData.of(selected_motif, 0)

    for seq_id, seq in enumerate(sequences):
        for window in sliding_window(seq, motif_length):
            motif = Bio.motifs.create([window] + selected_motif.instances)
            pwm = motif.counts.normalize(
                pseudocounts={"A": 0.6, "C": 0.4, "G": 0.4, "T": 0.6}
            )
            value = pwm.log_odds().std({"A": 0.3, "C": 0.2, "G": 0.2, "T": 0.3})
            if (current := MotifData(motif, value, seq_id)) > best:
                best = current

    return (
        consensus_with_unified_motif(
            without_ith_element(sequences, best.id), motif_length, best.motif
        )
        if len(sequences) > 1
        else best
    )


def consensus_single(sequences, motif_length, keep=5):
    motifs = BestKeeper(keep)

    first, *rest = sequences
    for window in sliding_window(first, motif_length):
        motif = Bio.motifs.create([window])
        motifs += [consensus_with_unified_motif(rest, motif_length, motif)]

    return motifs.into_list()


def consensus(sequences, motif_length=7, iterations=10, keep=5):
    best = BestKeeper(keep)

    for _ in range(iterations):
        permuted_sequences = list(sequences)
        random.shuffle(permuted_sequences)
        best += consensus_single(permuted_sequences, motif_length, keep=keep)

    return best.into_list()


def test1():
    # some, more of less, random sequences
    data = [
        Bio.Seq.Seq("ATGAGGTCTATATCGTACGATCATCGCG"),
        Bio.Seq.Seq("ACTGATCGATCGATCATCGATCATCATC"),
        Bio.Seq.Seq("ATCGATCGATCGATCATCGATGCTACGC"),
        Bio.Seq.Seq("ATCGATTCGATCATCGTACGTACTACAT"),
        Bio.Seq.Seq("ATCGATCGATCGATCGATCGATCGATCG"),
        Bio.Seq.Seq("CATCGATCGTACGATCTACGATCGTGAC"),
        Bio.Seq.Seq("CAGTGAGCTACGTACGATGCTACGTACG"),
        Bio.Seq.Seq("ATCGATGCATCGATCGATCGATCGATCT"),
    ]
    motifs = consensus(data, 7)
    print("Top motifs:", "=" * 80, sep="\n")
    for motif in motifs:
        print("Motif (pfm):", motif.motif.format("pfm"))
        print("Value:", motif.value)
        print("=" * 80)


if __name__ == "__main__":
    test1()
