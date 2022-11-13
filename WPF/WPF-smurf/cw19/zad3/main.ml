
(* 

Jaka jest minimalna liczba monet potrzebna do wydania n zł reszty, przyjmując, że monety mogą być przekazywane w obie strony. (Na przykład, żeby otrzymać 9 zł mogę dać 1 zł i otrzymać 10 zł, razem 2). W obrocie są dostępne monety o zadanych (w postaci tablicy) nominałach.
przeszukujemy graf, algorytmem BFS

wierzchołki = ile mamy pieniędzy do wydania jeszcze
krawędź z i do j, jeżeli istnieje nominał |i-j|
szukamy najkrótszej ścieżki z kwoty 0 do kwoty n
ograniczamy się do wierzchołków o numerach od min(0, n - max nominał) do n + max nominał
n = 99
nominały = 100, 500, 501 [ 0 -> 100 -> 600 -> 99 ]
n = 1 nominały = 2, 3 [ 0 -> 3 -> 1 ]
złożoność czasowa: O( liczba wierzchołków + liczba krawędzi w grafie ) = O(liczba wierzchołków * liczba nominałów) = O(liczba nominałów) * O(n+max nominał)
złożoność pamięciowa = O(n + max nominał)

*)
