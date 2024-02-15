from Bio.Align import PairwiseAligner, substitution_matrices
from Bio import SeqIO
from numpy import mean

# Zadania na laboratorium:

# 1. Zapoznaj się z metodami uliniowienia par sekwencji przy użyciu obiektu Bio.Align.PairwiseAligner, obiektami typu
# Alignment. (jeśli ktoś używa starszej wersji biopythona, obsługa macierzy substytucji może być jeszcze w module
# Bio.SubsMat.MatrixInfo. Warto wtedy zajrzeć do dokumentacji tej wersji. Opis wykorzystania macierzy substytucji
# jest na stronie 97. Nasz serwer jupyter już został uaktualniony do wersji 1.81)

# DONE

# 2. Wczytaj sekwencje DNA histonów histones.fa i czynników bZIP bzips.fa do pamięci, najlepiej przy pomocy modułu
# Bio.SeqIO. Pliki są też dostępne na naszym serwerze jupyter w folderze /.


def seq_to_obj(seqs, mapper):
    return [mapper(seq.seq) for seq in seqs]


def read_dna_histones_and_bzips(mapper=lambda x: str(x)):
    histones = SeqIO.parse("histones.fa", "fasta")
    bzips = SeqIO.parse("bzips.fa", "fasta")
    return seq_to_obj(histones, mapper), seq_to_obj(bzips, mapper)


# 3. Dokonaj porównań pomiędzy sekwencjami DNA białek histonowych i bzip – dla każdej pary policz oceny dla
# najlepszych globalnych i lokalnych uliniowień z prostą macierzą substytucji (-1 za mismatch, 1 za match) i
# afiniczną funkcją kary za przerwy (-1.0 za otwarcie przerwy i -0,5 za jej rozszerzenie). Wylicz średnią ocenę w
# ramach grupy bzip, w ramach grupy histonów i pomiędzy grupami.


def score(seq1, seq2, mode, match, mismatch, gap_open=None, gap_extend=None):
    aligner = PairwiseAligner()
    aligner.mode = mode
    aligner.match = match
    aligner.mismatch = mismatch
    if gap_open is not None and gap_extend is not None:
        aligner.open_gap_score = gap_open
        aligner.extend_gap_score = gap_extend

    alignments = aligner.align(seq1, seq2)
    return alignments[0].score


def simple_and_affine_global_and_local(seq1, seq2):
    return [
        score(
            seq1,
            seq2,
            mode=mode,
            match=1,
            mismatch=-1,
            gap_open=gap_open,
            gap_extend=gap_extend,
        )
        for gap_open, gap_extend in [(None, None), (-1.0, -0.5)]
        for mode in ["global", "local"]
    ]


def avg_for_seqs(seqs1, seqs2):
    all_scores = zip(
        *[
            simple_and_affine_global_and_local(seq1, seq2)
            for seq1 in seqs1
            for seq2 in seqs2
        ]
    )

    labels = [
        f"Mode: {mode}, Affine: {(gap_open, gap_extend) != (None, None)}"
        for gap_open, gap_extend in [(None, None), (-1.0, -0.5)]
        for mode in ["global", "local"]
    ]

    return list(zip(labels, [mean(scores) for scores in all_scores]))


def task3():
    histones, bzips = read_dna_histones_and_bzips()

    print("=" * 120)
    print("DNA")

    print("Histones vs Histones")
    print(avg_for_seqs(histones, histones))

    print("Bzips vs Bzips")
    print(avg_for_seqs(bzips, bzips))

    print("Histones vs Bzips")
    print(avg_for_seqs(histones, bzips))


# 4. Dokonaj tłumaczenia sekwencji DNA na białka(używając metody translate dla DNA), powtórz obliczenia z punktu 3,
# ale na sekwencjach białkowych używając macierzy substutucji BLOSUM80


def score_blosom80(seq1, seq2, mode, match, mismatch):
    aligner = PairwiseAligner()
    aligner.mode = mode
    aligner.match = match
    aligner.mismatch = mismatch
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM80")
    alignments = aligner.align(seq1, seq2)
    return alignments[0].score


def global_and_local_blosom80(seq1, seq2):
    return [
        score_blosom80(seq1, seq2, mode=mode, match=1, mismatch=-1)
        for mode in ["global", "local"]
    ]


def avg_for_seqs_blosom80(seqs1, seqs2):
    all_scores = zip(
        *[global_and_local_blosom80(seq1, seq2) for seq1 in seqs1 for seq2 in seqs2]
    )

    labels = [f"Mode: {mode}" for mode in ["global", "local"]]

    return list(zip(labels, [mean(scores) for scores in all_scores]))


def task4():
    histones, bzips = read_dna_histones_and_bzips(mapper=lambda seq: seq.translate())

    print("=" * 120)
    print("Translated")

    print("Histones vs Histones")
    print(avg_for_seqs_blosom80(histones, histones))

    print("Bzips vs Bzips")
    print(avg_for_seqs_blosom80(bzips, bzips))

    print("Histones vs Bzips")
    print(avg_for_seqs_blosom80(histones, bzips))


if __name__ == "__main__":
    task3()
    task4()
