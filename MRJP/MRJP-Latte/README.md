# Latte

## Zmiany po rozmowie
- Naprawiono błędy wypisane przez Panią w mailu
- Napisałem LCSE dla wyrażeń arytmetycznych dla intów (+, -, *, /, %)
- Niestety nie dałem rady napisać nic więcej :((

## Algorytm Phi
[https://c9x.me/compile/bib/braun13cc.pdf](https://c9x.me/compile/bib/braun13cc.pdf)

## Uruchomienie

Zbudowanie kompilatora:

```bash
make
```

Uruchomienie kompilatora:

```bash
./latc_llvm <plik.lat>
```

Projekt używa Cabala do zarządzania zależnościami.
`Makefile` uruchamia dodatkowo generowanie plików zależnych od gramatyki oraz buduje dodatkowe zależności w katalogu `lib`.

## Type Checker

Mniej oczywiste przypadki, które akceptuję:

- Pozwalam na stworzenie zmiennej o takiej samej nazwie jak zmienna pętli.

```java
int main() {
  for (int x : new int[1]) {
    int x;
  }
  return 0;
}
```

- Pozwalam na napisanie wartości zmiennej pętli. Zmiana wartości nie zmienia przebiegu pętli.

```java
int main() {
  for (int x : new int[1]) {
    x = 1;
  }
  return 0;
}
```
