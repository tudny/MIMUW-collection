# Java-C-performance-comparison

Program dla podanego grafu odpowiada wykonuje algorytm Floyda-Warsahlla, a następnie odpowiada offline na zapytania o najkrótsze ścieżki między parą wierzchołków.
Dla grafu `G = (V, E)` złożoność wynosi `O(|V|3+|E|)`.

* Porównywany język: `C`
* Środowisko: `x86_64 Linux 4.19.128-microsoft-standard` (WSL)
* Wersja Javy: `javac 11.0.10`
* Wersja C: `gcc version 9.3.0 (Ubuntu 9.3.0-17ubuntu1~20.04)`

Wyniki dla testu wygenerowanego generatorem (seed = 42 dla podanego środowiska i kompilatora): `|V| = 1463`, `|E| = 1615`

* Java (Just In Time): `0m9.545s`
* Java (-Xint): `0m9.773s`
* C(-O2): `0m2.931s`
* C: `0m12.933s`
