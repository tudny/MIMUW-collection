Dziś na wykładzie omawiamy zagadnienie mapowania krótkich odczytów do długiego genomu. Interesuje nas głównie
zrozumienie działania tranformaty Burrows'a-Wheeler'a i sposobu wyszukiwania w tekstach przy jej pomocy:

1. Napisz funkcję `compute_BWT(txt)`, zwracającą transformatę Burrows'a-Wheeler'a zadanego tekstu (jako napis), a także
   macierze `C(x)` i `OCC(i,x)` jako słowniki (dla alfabetu `ACTG`) - Ważne, aby podczas sortowania uniknąć używania
   kwadratowego rozmiaru pamięci, co dzieje się gdy np. wygenerujemy wszystkie sufiksy naszego tekstu...

2. Zaaplikuj tę funkcję do poniższej sekwencji:
   `CGAGCCGCTTTCCATATCTATTAACGCATAAAAAACTCTGCTGGCATTCACAAATGCGCAGGGGTAAAACGTTTCCTGTAGCACCGTGAGTTATACTTTGT`

3. zaimplementuj mapowanie Last-to-First i korzystającą z niego metodę znajdującą wystąpienia krótkich podsłów w
   tekście na podstawie `BWT(T)`. Wyszukaj słowa "AAAC" w powyższym tekście.

4. Korzystając z metody z wykładu dzielenia zapytania (query) na pół, napisz funkcję, która pozwoli Ci znaleźć
   wystąpienie słowa `AAAACTCCGCTGGCATTCACAAAT` w powyższym tekście z co najwyżej 1 błędem.

praca domowa (bonus) Wykorzystaj funkcje z dzisiejszego labu, aby wyszukać w genomie E. coli (
Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.dna.chromosome.Chromosome.fa) pozycje wystąpień promotorów (
ecoli_proms.fa ) UWAGA: powinniśmy unikać przechowywania w pamięcie struktur o rozmiarze kwadratowym względem rozmiaru
wejścia (np. nie zmieści nam się w pamięci tablica wszystkich sufiksów)