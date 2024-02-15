import hmmlearn.hmm as hmm
import numpy as np
from Bio import SeqIO

# Dziś na wykładzie omówimy ukryte modele Markowa i algorytmy rekonstrukcji parametrów z danych.
#
# Zadania na dziś:

# 0. Zapoznaj się z modułem hmmlearn. Będziemy go wykorzystywać do uczenia modeli z emisjami

# 1. Wyspy CpG znajdują się często w genomach, w szczególności genomie ludzkim. Spróbuj zdefiniować ukryty model
# Markowa, który ma 2 stany i spróbuj nauczyć go na sekwncji zawartej w pliku cpg.fa. Zrób to zarówno dla sekwencji
# liter (ACGT), jak i dla sekwencji dwunukleotydów (AA,AC,AG,AT, itp…) Czy możesz zinterpretować macierze emisji i
# przypisać jeden ze stanów do wysp CpG? Wykonaj kilkakrotnie proces uczenia (Baum-Welch) i zobacz czy wyniki są
# podobne. Jak interpretujesz prawdopodobieństwa w macierzy przejść. Czy coś możesz powiedzieć o średniej długości
# wysp CpG?
#
# Warto przyjrzeć się przykładowi użycia klasy CategoricalHMM


def dna_to_ordinal(dna: str) -> int:
    return {"A": 0, "C": 1, "G": 2, "T": 3}[dna]


full_training_seq = SeqIO.read("data/cpg.fa", "fasta").seq
full_training_seq = [dna_to_ordinal(dna) for dna in full_training_seq]

# Weźmy 8 stanów ukrytych i 4 emisje (A, C, G, T)
# Stany podzielimy na dwie kategorie - nie-CpG i CpG - po cztery stany w każdej kategorii
# Pierwszy kategoria, nie-CpG, będzie miała stany 0, 1, 2, 3
# Druga kategoria, CpG, będzie miała stany 4, 5, 6, 7
# Emisje będą odpowiadały literom A, C, G, T

p_state_switch = 0.01
p4 = (1 - p_state_switch) / 4
ps4 = p_state_switch / 4
p_cg = 0.4
ps4_cg = p_state_switch / 4
p4_cg = (1 - p_state_switch - p_cg) / 3

states_matrix = np.array(
    [
        # (nCpG)  A   C   G   T  | (CpG)  A    C    G    T
        [p4, p4, p4, p4, ps4, ps4, ps4, ps4],
        [p4, p4, p4, p4, ps4, ps4, ps4, ps4],
        [p4, p4, p4, p4, ps4, ps4, ps4, ps4],
        [p4, p4, p4, p4, ps4, ps4, ps4, ps4],
        [ps4, ps4, ps4, ps4, p4, p4, p4, p4],
        [ps4_cg, ps4_cg, ps4_cg, ps4_cg, p4_cg, p4_cg, p_cg, p4_cg],
        [ps4, ps4, ps4, ps4, p4, p4, p4, p4],
        [ps4, ps4, ps4, ps4, p4, p4, p4, p4],
    ]
)

emissions_matrix = np.array(
    [
        # A  C  G  T
        [1, 0, 0, 0],
        [0, 1, 0, 0],
        [0, 0, 1, 0],
        [0, 0, 0, 1],
        [1, 0, 0, 0],
        [0, 1, 0, 0],
        [0, 0, 1, 0],
        [0, 0, 0, 1],
    ]
)

cpg_dna_hmm = hmm.CategoricalHMM(n_components=8, init_params="")
cpg_dna_hmm.n_features_ = 4
cpg_dna_hmm.startprob_ = np.array([0.125] * 8)
cpg_dna_hmm.transmat_ = states_matrix
cpg_dna_hmm.emissionprob_ = emissions_matrix


training_reshaped = np.array(full_training_seq).reshape(1, -1)
cpg_dna_hmm.fit(training_reshaped)
# We need to adjust start probabilities, because the training set will favor the starting state of the training set,
# and as the training set have only one sample (the whole sequence), the start probabilities will be 1.0 for the
# starting state of the training set, and 0.0 for the rest
cpg_dna_hmm.startprob_ = np.array([0.125] * 8)

print(cpg_dna_hmm.transmat_)
print(cpg_dna_hmm.emissionprob_)
print(cpg_dna_hmm.startprob_)

# praca domowa (2 pkt):
#
# Wykorzystaj model nauczony na danych o CpG i przetestuj które z 30 sekwencji w pliku cpg_test.fa są naprawdę
# wyspami CpG. Jako wynik proszę przysłać program w pythonie i wynik w pliku tekstowym (w każdej linii możemy podać
# parawdopodobieństwo a'posteriori tego, że dana sekwencja z wejścia pochodzi z modelu CpG).

request_seqs = SeqIO.parse("data/cpg_test.fa", "fasta")


def map_to_cpg_state(state: int) -> int:
    if state in [0, 1, 2, 3]:
        return 0
    elif state in [4, 5, 6, 7]:
        return 1
    else:
        raise ValueError(f"Unknown state: {state}")


def count_cpg(states: [int]) -> int:
    assert all(map(lambda s: s in [0, 1], states))
    return sum(states)


for req_seq in request_seqs:
    req_name = req_seq.description
    req_seq = req_seq.seq
    req_seq = [dna_to_ordinal(dna) for dna in req_seq]
    req_seq = np.array(req_seq).reshape(1, -1)
    result = cpg_dna_hmm.decode(req_seq)
    res_len = len(result[1])
    cpg_states = [map_to_cpg_state(state) for state in result[1]]
    cpg_count = count_cpg(cpg_states)
    cpg_ratio = cpg_count / res_len
    is_cpg = cpg_ratio > 0.5
    print(
        f"Sequence: {req_name} "
        f"| Ratio: {cpg_ratio:.4f} "
        f"| {'CpG island' if cpg_ratio > 0.5 else 'Not CpG island'}"
    )
