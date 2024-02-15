#!/usr/bin/env python3
# ============================================================
# | This is a solution to the problem                        |
# | to find the best local alignment                         |
# | between two sequences ignoring                           |
# | stop codons.                                             |
# ============================================================
# | Author: Aleksander Tudruj                                |
# | Index:  at429630                                         |
# | MIMUW WBO 2023                                           |
# ============================================================
# | Usage:                                                   |
# |   ./best_local_alignment_at429630.py -f <filename>       |
# |   ./best_local_alignment_at429630.py -s                  |
# |                                                          |
# |  -f <filename> - read sequences from fasta file          |
# |  -s            - read sequences from stdin               |
# |  -h           - show help                                |
# |                                                          |
# | OR                                                       |
# |  Use function `best_local_alignment_without_stop_codons` |
# |  to get the best local alignment between two sequences   |
# |  ignoring stop codons.                                   |
# ============================================================

import argparse
from sys import stderr
from Bio import SeqIO, Seq
from Bio.Align import substitution_matrices, PairwiseAligner, Alignment
from Bio.Data import CodonTable


def _flatten(iterables):
    return [item for sublist in iterables for item in sublist]


def _chunks(lst, n):
    for i in range(0, len(lst), n):
        yield lst[i : i + n]


class AlignmentResult:
    def __init__(self, score, seq1, seq2):
        self.score = score
        self.seq1 = seq1
        self.seq2 = seq2

    def __repr__(self) -> str:
        return f"AlignmentResult({self.score}, {self.seq1}, {self.seq2})"

    def result(self) -> (Seq, Seq):
        return self.seq1, self.seq2


def _get_best_alignments(alignments: [AlignmentResult]) -> [AlignmentResult]:
    if len(alignments) < 1:
        return []
    best_score = max(alignments, key=lambda x: x.score)
    return list(filter(lambda x: x.score == best_score.score, alignments))


def _amino_alignment_to_dna(alignment: Alignment, seq: Seq) -> Seq:
    result = Seq.MutableSeq("")
    for kodon in alignment:
        if kodon is None:
            result += "---"
        else:
            result += seq[:3]
            seq = seq[3:]
    return Seq.Seq(result)


def _align_locally(seq1: Seq, seq2: Seq, codon_table) -> [AlignmentResult]:
    seq1_translated = seq1.translate(table=codon_table)
    seq2_translated = seq2.translate(table=codon_table)
    matrix = substitution_matrices.load("BLOSUM62")
    aligner = PairwiseAligner()
    aligner.substitution_matrix = matrix
    aligner.mode = "local"
    alignments = aligner.align(seq1_translated, seq2_translated)
    return [
        AlignmentResult(
            alignments.score,
            _amino_alignment_to_dna(alignment[0], seq1),
            _amino_alignment_to_dna(alignment[1], seq2),
        )
        for alignment in alignments
    ]


def _split_by_stop_codon(seq1, codon_table):
    if (offset := len(seq1) % 3) != 0:
        seq1 = seq1[:-offset]

    parted = [Seq.MutableSeq("")]
    for chunk in _chunks(seq1, 3):
        if chunk not in codon_table.stop_codons:
            parted[-1] += chunk
        else:
            parted.append(Seq.MutableSeq(""))
    parted = [part for part in parted if len(part) > 0]
    return [Seq.Seq(part) for part in parted]


def _opt_local_alignment(seq1: Seq, seq2: Seq, codon_table) -> [AlignmentResult]:
    seq1_split = _split_by_stop_codon(seq1, codon_table)
    seq2_split = _split_by_stop_codon(seq2, codon_table)

    alignments = [
        _align_locally(seq1_part, seq2_part, codon_table)
        for seq1_part in seq1_split
        for seq2_part in seq2_split
    ]

    return _get_best_alignments(_flatten(alignments))


def best_local_alignment_without_stop_codons(
    seq1: Seq, seq2: Seq, codon_table=CodonTable.standard_dna_table
) -> (Seq, Seq):
    """
    Find the best local alignment between two sequences ignoring stop codons.
    :param seq1: the first sequence
    :param seq2: the second sequence
    :param codon_table: the codon table to use for translation (default: standard_dna_table)
    :return: the best local alignment between seq1 and seq2 ignoring stop codons
    """
    offsets = [0, 1, 2]

    alignments = [
        _opt_local_alignment(seq1[pad1:], seq2[pad2:], codon_table)
        for pad1 in offsets
        for pad2 in offsets
    ]

    return _get_best_alignments(_flatten(alignments))[0].result()


# =======================================
# |             IO Utils                |
# =======================================


class SequenceProvider:
    def get_sequences(self):
        raise NotImplementedError()


class FileSequenceProvider(SequenceProvider):
    def __init__(self, filename):
        self.filename = filename

    def get_sequences(self):
        all_seqs = SeqIO.parse(self.filename, "fasta")
        return [seq.seq for seq in all_seqs][:2]


class StdioSequenceProvider(SequenceProvider):
    def get_sequences(self):
        return [
            input(f"Enter {ordinal} sequence\n").upper()
            for ordinal in ["first", "second"]
        ]


def run(sequence_provider):
    seq1, seq2 = sequence_provider.get_sequences()
    seq1_local, seq2_local = best_local_alignment_without_stop_codons(seq1, seq2)
    return str(seq1_local), str(seq2_local)


def print_result(result):
    seq1, seq2 = result
    print("Sequences:")
    print(seq1)
    print(seq2)


def main():
    parser = argparse.ArgumentParser(
        prog="aligner",
        description="Find the best alignment between two translated sequences.",
        epilog="MIMUW@2023 Aleksander Tudruj at429630",
    )
    parser.add_argument(
        "-f",
        "--file",
        type=str,
        help="File with sequences. The first two are taken. Rest is ignored.",
    )
    parser.add_argument(
        "-s", "--stdio", action="store_true", help="Read sequences from stdin"
    )
    args = parser.parse_args()

    predicates = [
        (args.file and args.stdio, "You must provide only one of -f or -s argument"),
        (not args.file and not args.stdio, "You must provide either -f or -s argument"),
    ]

    for predicate, message in predicates:
        if predicate:
            print(message, file=stderr)
            parser.print_help()
            exit(1)

    if args.file:
        provider = FileSequenceProvider(args.file)
    elif args.stdio:
        provider = StdioSequenceProvider()
    else:
        raise ValueError("Invalid arguments")

    result = run(provider)
    print(result)


if __name__ == "__main__":
    main()
