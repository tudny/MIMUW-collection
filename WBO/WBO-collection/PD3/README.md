# Praca domowa (2pkt):

Napisz program, który wylicza optymalne lokalne uliniowienie dla sekwencji DNA dla różnych możliwych tłumaczeń na
białka (tzn zakładając standardową tablicę kodonów oraz macierz substytucji BLOSUM62)  przy uwzględnieniu trzech różnych
ramek odczytu w obu sekwencjach. Dla uproszczenia załóżmy, że insercje i delecje powinny występować tylko “trójkami”
nukleotydów. Kodony stopu nie mogą być częścią roziwązania (lokalnego uliniowienia), ale lokalne uliniowienie może
pochodzić z ramki odczytu z kodonem stopu, jeśli całe mieści się przed lub po kodonie stop.(Prace domowe jak zwykle
przesyłamy do następnych zajęć przez moodle – proszę załączyć rozwiązanie jako pliki .py lub linki do działających
notebooków, nie jako pliki .ipynb)

Np. dla sekwencji "ATTGGCGAGCACA" i "GATAGATAGCTAGCTAGTA" optymalne uliniowienie może wyglądać tak ("TTGGCGAGC","
CTAGCTAGT")

(*to już nie na punkty) Zastanów się, jakby wyglądał algorytm programowania dynamicznego w tym zadaniu gdybyśmy
rozważali ogólną postać insercji i delecji, a nie tylko wielokrotności 3.