# Projekt zaliczeniowy 2 (20 pkt)

## Treść

Tym razem nasze zadanie polegać będzie na napisaniu programu, który dla zadanej listy fragmentów białek uzupełnia je do pełnych sekwencji używając programu BLAST, a następnie wyszukuje w tych sekwencjach domen przy pomocy algorytmu HMMER, oraz wyszukuje motywów DNA w promotorach tych sekwencji (w grupach A i B). 

1. (6 pkt) Napisz program extend.py, który dla zadanej listy fragmentów białek w formacie fasta (protein_fragments.fa), znajdzie w bazie lokalnej stworzonej z podanego większego pliku fasta (genes_e_coli.fa) najbliższe białko i zwróci nowy plik fasta z wynikowymi sekwencjami genów już przetłumaczonych na białka. (jak uruchamiać blast lokalnie). Proszę zwrócić uwagę, że mamy fragmenty białek i sekwencje DNA genów.

2. (7 pkt) Napisz program scan_pfam.py, który na podstawie pliku z białkami wykonuje zapytanie do serwera hmmscan (serwer online tu, opis api) i pobiera pliki wynikowe w formacie tsv. Korzystając z plików pobranych z serwera hmmer poda nam w wyniku plik csv, w którym będziemy mieli w wierszach kolejne identyfikatory białek z pliku FASTA, zaś w kolumnach będzie miał kolejne identyfikatory domen białkowych PFAM. Na przecięciu wiersza i kolumny stawiamy 0, jeśli dana domena nie została znaleziona w danym białku, a 1 w przeciwnym wypadku. Niestety moduły do parsowania wyjścia z HMMera w biopythonie często sprawiają problemy, dlatego polecam ręcznie wczytywać pliki tsv.

3. (7 pkt ) Napisz program, który dla genów znalezionych na podstawie fragmentów z grup A i B (w oryginalnym pliku fasta z fragmentami, są one opisane) wyszuka po 5 najelpszych motywów długości 10 (parametr Państwa programu) w promotorach. (Promotory genów e-coli są w pliku proms_e_coli_fixed.fa), Można wykorzystać własną metodę consensus lub program MEME.

