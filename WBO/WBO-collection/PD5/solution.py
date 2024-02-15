import os

from Bio import Seq, SeqIO, AlignIO, Phylo
from Bio.Align import substitution_matrices
from Bio.Align.Applications import ClustalwCommandline

import numpy as np


# 6. Praca domowa (2 pkt). Zaimplementuj prosty algorytm uliniowienia progresywnego. Jako dane, Twoja funkcja powinna
# przyjmować listę sekwencji i drzewo (np. pochodzące z metody UPGMA albo nj) oraz macierz substytucji i wartość kary
# za przerwy. Funkcja powinna liczyć koszt substytucji dla profili wg miary sumy par i zwracać multiulioniowienie.


def simple_progressive_alignment(seqs_list, tree, matrix, gap_open, gap_extend):
    result, names = align_node(tree.root, seqs_list, matrix, (gap_open, gap_extend))
    mapped_result = {name: seq for name, seq in zip(names, result)}
    sorted_result = [(seq.id, mapped_result[seq.id]) for seq in seqs_list]
    return sorted_result


def get_len(seqs):
    if len(set(map(len, seqs))) != 1:
        raise ValueError("Sequences must be the same length")
    return len(seqs[0])


MATCH = 1
GAP1 = 2
GAP2 = 4
GAPB = 6


def align_children(seqs1, seqs2, gap_cost, matrix):
    n = get_len(seqs1)
    m = get_len(seqs2)

    (gap_open, gap_extend) = gap_cost

    dp_matrix = np.zeros((n + 1, m + 1), dtype=float)
    dp_source = np.zeros((n + 1, m + 1), dtype=int)
    seqs1_res = ["" for _ in range(len(seqs1))]
    seqs2_res = ["" for _ in range(len(seqs2))]

    def get_gap_cost(id_i, id_j, kind):
        if kind == "i" and (dp_source[id_i, id_j] & GAP1 > 0):
            return gap_extend
        if kind == "j" and (dp_source[id_i, id_j] & GAP2 > 0):
            return gap_extend
        return gap_open

    def get_cost(amino1, amino2, id_i, id_j):
        if amino1 == "-" or amino2 == "-":
            return get_gap_cost(id_i, id_j, GAP1 if amino1 == "-" else GAP2)
        return matrix[amino1, amino2]

    dp_matrix[0, :] = gap_extend
    dp_matrix[:, 0] = gap_extend
    dp_matrix[0, 0] = 0
    dp_source[0, :] = GAP1
    dp_source[:, 0] = GAP2
    dp_source[0, 0] = GAPB

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            match = dp_matrix[i - 1, j - 1] + sum(
                [
                    get_cost(seq1[i - 1], seq2[j - 1], i - 1, j - 1)
                    for seq1 in seqs1
                    for seq2 in seqs2
                ]
            )
            gap_seq2 = dp_matrix[i - 1, j] + get_gap_cost(i - 1, j, "j") * m
            gap_seq1 = dp_matrix[i, j - 1] + get_gap_cost(i, j - 1, "i") * n

            dp_matrix[i, j], dp_source[i, j] = max(
                zip([match, gap_seq1, gap_seq2], [MATCH, GAP1, GAP2]),
                key=lambda x: x[0],
            )

    i_it = n
    j_it = m
    while i_it > 0 or j_it > 0:
        if dp_source[i_it, j_it] == MATCH and i_it > 0 and j_it > 0:
            for seq_it in range(len(seqs1)):
                seqs1_res[seq_it] = seqs1[seq_it][i_it - 1] + seqs1_res[seq_it]
            for seq_it in range(len(seqs2)):
                seqs2_res[seq_it] = seqs2[seq_it][j_it - 1] + seqs2_res[seq_it]
            i_it -= 1
            j_it -= 1
        elif dp_source[i_it, j_it] & GAP2 > 0 and i_it > 0:
            for seq_it in range(len(seqs1)):
                seqs1_res[seq_it] = seqs1[seq_it][i_it - 1] + seqs1_res[seq_it]
            for seq_it in range(len(seqs2)):
                seqs2_res[seq_it] = "-" + seqs2_res[seq_it]
            i_it -= 1
        elif dp_source[i_it, j_it] & GAP1 > 0 and j_it > 0:
            for seq_it in range(len(seqs1)):
                seqs1_res[seq_it] = "-" + seqs1_res[seq_it]
            for seq_it in range(len(seqs2)):
                seqs2_res[seq_it] = seqs2[seq_it][j_it - 1] + seqs2_res[seq_it]
            j_it -= 1
        else:
            raise ValueError("Something went wrong")

    return seqs1_res + seqs2_res


def align_node(node, seqs, matrix, gap_cost):
    if len(node.clades) == 0:
        return [str(seq.seq) for seq in seqs if seq.id == node.name], [
            seq.id for seq in seqs if seq.id == node.name
        ]

    children = [align_node(child, seqs, matrix, gap_cost) for child in node.clades]

    while len(children) > 1:
        (a_seq, a_names), (b_seq, b_names) = children.pop(), children.pop()
        res = align_children(a_seq, b_seq, gap_cost, matrix), a_names + b_names
        children.append(res)
    return children[0]


def read_alignment(path: str, source: str):
    def translate_seq(seq: Seq) -> Seq:
        seq.seq = seq.seq.translate(stop_symbol="")
        return seq

    sequences = list(SeqIO.parse(f"{path}/{source}.fa", "fasta"))
    translated_sequences = list(map(translate_seq, sequences))
    SeqIO.write(translated_sequences, f"{path}/{source}_translated.fa", "fasta")

    ClustalwCommandline(
        cmd="clustalw",
        infile=f"{path}/{source}_translated.fa",
        outfile=f"{path}/{source}_translated_clustalw.aln",
    )()

    os.rename(  # rename files to match clustalw output
        f"{path}/{source}_translated.dnd", f"{path}/{source}_translated_clustalw.dnd"
    )

    alignment = AlignIO.read(
        f"{path}/{source}_translated_clustalw.aln", format="clustal"
    )
    phylo_tree = Phylo.read(f"{path}/{source}_translated_clustalw.dnd", format="newick")
    return alignment, phylo_tree, sequences


def test():
    result_alignment, phylo_tree, sequences = read_alignment("data", "histones")
    matrix = substitution_matrices.load("PAM250")
    sorted_result = simple_progressive_alignment(
        sequences, phylo_tree, matrix, -1.0, -0.1
    )

    print(matrix)

    print("\n\n", "RESULT:", "\n\n")
    for seq in sorted_result:
        print(seq[0], ": ", seq[1])

    print("\n\n", "EXPECTED:", "\n\n")
    for seq in result_alignment:
        print(seq.id, ": ", seq.seq)

    print("\n\n", "ORIGINAL:", "\n\n")
    for seq in sequences:
        print(seq.seq)


if __name__ == "__main__":
    test()
