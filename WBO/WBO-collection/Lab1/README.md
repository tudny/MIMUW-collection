# Laboratorium 1

Zasadniczo pracujemy w języku python. Jeśli ktoś ma kłopoty z interpreterem w labie, to może korzystać z serwera
jupyter, używając hasła: "studenci-wbo-2023". Jeśli korzystają Państwo z serwera, proszę o założenie własnego,
podpisanego, pliku w unikalnym katalogu.

W laboratorium zajmiemy się następującymi kwestiami:

1. Spróbujemy wczytać zestaw sekwencji z pliku FASTA przy pomocy biblioteki Bio.Seq z pakietu Biopython. Mamy do
   dyspozycji plik test_fasta Na początek warto go wczytać i wypisać do innego pliku zestaw sekwencji odwrotnie
   komplementarnych (np. metodą .reverse_complement ) do sekwencji wejściowych. Warto zapoznać się z tutorialem

2. Napiszemy funkcje kmers(s,k), która na podstawie sekwencji DNA daje nam wszystkie k-mery (podsekwencje długości k)
   składające się na widmo sekwencji s. Rozważymy dwa warianty – gdy widmo zawiera tylko podsekwencje s, i drugą, gdzie
   do widma zaliczamy też fragmenty sekwencji komplementarnej.

3. Napiszemy funkcje euler(kmers) i hamilton(kmers), które tworzą graf na podstawie widma k-merów i znajdują w nim
   wszystkie ścieżki (ta funkcja może być naiwna, jeśli chodzi o złożoność). Graf można reprezentować np. jako słownik,
   w którym każdemu k-merowi odpowiadają te, do których prowadzi krawędź grafu.

4. Spróbujemy rozwiązać problem znajdowania zestawu sond komplementarnych z ostatniego slajdu wykładu.

## Praca domowa (2pkt.)

Napisz program, który znajduje minimalne k, dla którego istnieje zbiór sond długości k komplementarnych do zestawu
sekwencji (każda sonda komplementarna do dokładnie jednej sekwencji z zestawu) zadanego w pliku .fasta. Rozwiązanie
polega na przysłaniu programu (z komentarzami), który znajduje minimalne k oraz wyniku jego uruchomienia dla wybranych
genów drożdży Saccharomyces cerevisiae yeast.fa.

Np. dla pliku przykładowego test_fasta, rozwiązaniem jest k=3 (jednym ze zbiorów P jest ["CGT", "GCA", "TTT", "GGG"]

Prace domowe trzeba przesłać przez moodle przed następnymi zajęciami w laboratorium. Należy przesłać kod rozwiązania i
wynikowe k dla podanego zbioru danych "drożdżowych". (np. w komentarzu).

4. Zastanów się nad problemem zestawu sond komplementarnych, gdzie zamiast równej długości wszystkich sond założymy, że
powinny mieć tę samą temperaturę topnienia.




