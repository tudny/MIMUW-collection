# JPP | Zadanie 2 | interpreter

### Nazwa języka
`Jabba the Lang`
### Rozszerzenie plików
`.jbb`

## Specyfikacja języka

Typ języka: imperatywny

Język nie posiada szczególnych udziwnień w porównaniu do Latte.

### Typy, zmienne, stałe, literały, przekazywanie argumentów

Języka trzy typy: `Integer`, `Boolean`, `String`. Każda może być zadeklarowana jako `var` lub `val`. Pierwsza z nich oznacza zmienną możliwą do modyfikacji, druga oznacza stałą. Wszystkie zmienne są inicjalizowane wartością domyślną. `Integer` = `0`, `Boolean` = `false`, `String` = `""`. Wszystkie zmienne są typowane statycznie. W przypadku próby przypisania zmiennej typu `Integer` wartości typu `String` zostanie zgłoszony wyjątek.

Każdy z trzech typów posiada literały. `Integer` - liczby całkowite, `Boolean` - `true` i `false`, `String` - ciągi znaków w cudzysłowach.

Na liczbach można wykonywać operacje arytmetyczne: `+`, `-`, `*`, `/`, `%`. Na `Boolean` można wykonywać operacje logiczne: `&&`, `||`, `!`. Na `String` można wykonywać operacje konkatenacji: `+`. Na `Integer` można wykonywać operacje porównania: `==`, `!=`, `<`, `>`, `<=`, `>=`.

Standardowe przekazania zmiennej do funkcji odbywa się przez referencję (nie kopię referencji). 
Przykład:
```kotlin
fun foo(var x: Integer) : Unit {
    x = 10;
}
var y: Integer = 5;
foo(y);
// x = 10, ponieważ x zostało przekazane przez referencję
```
Przy czym gdyby `y` było typu `val`, to wtedy program zgłosiłby błąd na etapie sprawdzania typów.
Aby przekazać kopię należy dopisać słowo kluczowe `new` przed typem.
Przykład:
```kotlin
fun foo(var x: new Integer) : Unit {
    x = 10;
}
var y: Integer = 5;
foo(y);
// x = 5, ponieważ x zostało przekazane przez wartość
```

Wymusza to na użytkowniku zastanowienia się przed skopiowaniem zmiennej, ponieważ jest to kosztowna operacja (może nie w przypadku naszych prostych typów, ale w ogólności).

W przypadku zwykłych operacji na zmiennych zawsze dochodzi do kopiowania. Nie pozwalamy na referencje poza wywołaniami funkcji. Przykład:
```kotlin
var x: Integer = 5;
var y: Integer = x;
y = 10;
// x = 5, ponieważ x zostało przekazane przez wartość
```

## Funkcje standardowe (w tym print)

```kotlin
// wypisuje `str` na stdout
writeStr(str: String) : Unit

// wypisuje `num` na stdout
writeInt(num: Integer) : Unit

// zamienia zmienną typu `Integer` na `String`
toString(Integer) : String

// zamienia zmienną typu `String` na `Integer`
toInt(String) : Integer

// sprawdza warunek i w przypadku jego nie spełnienia wywłaszcza program
assert(Boolean) : Unit

// zamyka program z kodem błędu podanym w argumencie
exit(Integer) : Unit
```

### Instrukcje warunkowe, pętle, break, continue, while-finally
Język posiada instrukcje warunkowe `if` (`if else`) oraz pętle `while` oraz `for`. Każda z instrukcji wymaga bloku kodu `{}`,
ale nie wymaga nawiasów okrągłych `()`. Przykład:
```kotlin
if true {
    writeInt(5);
}
// wyświetli 5

if true {
    writeInt(5);
} else {
    writeInt(10);
}
// wyświetli 5

while x < 10 {
    writeInt(x);
    x = x + 1;
} // wyświetli 5, 6, 7, 8, 9
```
W bloku pętli `while` oraz `for` można dodać instrukcję `break` lub `continue`.
W przypadku pętli `while` można dodać blok `finally` na końcu, który wykona się po wyjściu
z pętli, jeżeli nie wystąpiło wywołanie instrukcji `break`.
 Przykład:
```kotlin
while x < 10 {
    if x == 5 {
        break;
    }
    writeInt(x);
    x = x + 1;
} // wyświetli 0, 1, 2, 3, 4

x = 0;
while x < 2 {
    if x == 1 {
        break;
    }
    x = x + 1;
} finally {
    writeStr("finally");
} // nic nie wyświetli

x = 0;
while x < 2 {
    if x == 10 {
        break;
    }
    x = x + 1;
} finally {
    writeStr("finally");
} // wyświetli "finally"
```

Pętle `for` są dwóch rodzajów. Range loop oraz generator loop. Range loop wygląda następująco:
```kotlin
for x in 0..10 {
    writeInt(x);
} // wyświetli 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10
```
Pętla `for x in e1..e2` jest równoważna pętli `while e1 <= e2 { x = e1; e1 = e1 + 1; }`. Stała `x` jest typu `Integer` i jest widoczna tylko wewnątrz pętli. 

Generator loop wygląda następująco:
```kotlin
for x in gen() {
    writeInt(x);
} // wyświetli 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10
```
Pętla `for x in gen()` jest równoważna pętli `val g: Gen[T] = gen(); while hasNext(g) { x = next(g); }`. Stała `x` jest typu `T` i jest widoczna tylko wewnątrz pętli, gdzie T to typ zwracany przez generator g.

### Funkcje

Funkcje można deklarować na dowolnym poziomie zagnieżdżenia. Funkcje mogą zwracać wartość lub nie. Funkcje mogą przyjmować dowolną liczbę argumentów. Przykład:
```kotlin
fun foo(val x: Integer, val y: Integer) : Integer {
    return x + y;
}

fun bar(val x: Integer, val y: Integer) : Unit {
    writeInt(x + y);
}
```
Widoczność zmiennych dla funkcji jest taka sama jak dla bloków kodu `{}`. Przykład:
```kotlin
fun foo(val x: Integer, val y: Integer) : Integer {
    var z: Integer = 5;
    return x + y + z;
}
// zmienna z jest widoczna tylko wewnątrz funkcji foo
```
```kotlin
fun foo(val x: Integer, val y: Integer) : Integer {
    var z: Integer = 5;
    fun bar(val x: Integer, val y: Integer) : Integer {
        return x + y + z; // z widoczny jest z funkcji foo; x, y przysłonięte 
    }
    return bar(x + 1, y + 1) + z;
}
foo(1, 2); // zwróci 15
```

### Błędy składni, typów i wykonania

Błędna składnia (parser, lexer) zostaną sprawdzone jako pierwsze. 
W przypadku błędu składni program nie zostanie uruchomiony.
Typy zostaną sprawdzone przed wykonaniem programu. Błędy typów to na przykład:
- przypisanie wartości złego typu do zmiennej (w tym literałów)
  - `var x: Integer = "abc";` albo `var x: Integer = 10; x = "abc";`
- próba zmiany stałej
  - `val c: Integer = 10; c = 20;`
- przekazanie błędnego typu do funkcji
  - `fun foo(val x: Integer) : Integer { return x; } foo("abc");`
- przekazanie złej liczby argumentów do funkcji
  - `fun foo(val x: Integer) : Integer { return x; } foo(1, 2);`
- przekazanie stałej w miejsce zmiennej do funkcji
  - `fun foo(var x: Integer) : Integer { return x; } val y: Integer = 10; foo(y);`
- sprawdzenie zwracanego typu przez funkcję
  - `fun foo(val x: Integer) : String { return x; }`
- sprawdzenie zwracanych typów przez generator
  - `fun foo(val x: Integer) : Gen[String] { yield x; }`
- użycie niezadeklarowanej zmiennej
  - `var x: Integer = 10; x = x + y;`

W trakcie wykonania programu mogą wystąpić błędy wykonania. Błędy wykonania to na przykład:
- dzielenie przez 0 (przy czym dzielenie przez literał 0 jest sprawdzane podczas typowania)

### ~~Generatory~~

Zdecydowałem się nie pisać generatorów. Zamiast tego powstały lambdy.

### Lambdy i przekazywanie funkcji jako argumentów

Lambdy są anonimowymi funkcjami. Przykład:
```kotlin
var g: (var String) -> String = |var x: String| "World";

g("Hello"); // zwróci "World"
```

Funkcje mogą być przekazywane jako argumenty do funkcji. Przykład:
```kotlin
fun foo(val x: new Integer, val f: (val$Integer) -> Integer) : Integer {
    return f(x);
}

fun bar(val x: new Integer) : Integer {
    return x + 1;
}

val res = foo(10, bar); // zwróci 11
assert(res == 11);
```

Lambdy mogę nie przyjmować argumentów. Przykład:
```kotlin
var g = || 10;

for i = 0..100 {
    assert(g() == 10);
    writeInt(g());
    writeStr("\n");
}
```

Syntax dla lambd:
```kotlin
|<argumenty>| <ciało>
```

Argumenty lambdy mogą być zadeklarowane jako `val` lub `var`. Przykład:
```kotlin
var g = |var x: Integer| x + 1;
```

Tak samo jak w przypadku funkcji argument lambdy może być zadeklarowany jako `new`. Przykład:
```kotlin
var g = |new x: Integer| x + 1;
```

Typowanie zmiennych funkcyjnych jest dosyć proste. Przykład:
```kotlin
var g: (var Integer) -> Integer = |var x: Integer| x + 1;

var x = 10;
g(x); // zwróci 11
```
Funkcja przyjmuje zmienną `Integer` i zwraca `Integer`.
Zmienna przyjęta w funkcji `g` jest zmienną lokalną dla lambdy.
Jest ona modyfikowalną referencją do zmiennej przekazanej do lambdy.

Tutaj lambda przyjmuje dwa argumenty, z czego jeden jest kopią 
(`$` jest odpowiednikiem `new` w typach lambd).
```kotlin
fun foo(
    val f: new (var Integer, var$Integer) -> Unit, 
    var x: Integer, 
    var y: new Integer
): Unit {make
    f(x, y);
}

var x = 10;
var y = 20;
foo(|var x: Integer, var y: new Integer| {
    val temp = x;
    x = y;
    y = temp;
}, x, y);

assert(x == 20);
assert(y == 20); // y jest kopią, więc nie zmieni się
```

Jak widać na przykładzie powyżej, lambdy mogą mieć dwa rodzaje ciała:
- wyrażeniowe
- blokowe

Wyrażeniowe lambdy zwracają wartość ostatniego wyrażenia. Przykład:
```kotlin
var g = |var x: Integer| x + 1;
```

Blokowe lambdy zwracają wartość wskazaną przez instrukcję `return`. Przykład:
```kotlin
var g = |var x: Integer| {
    return x + 1;
};
```


### Dołączone pliki
Dodatkowe przykłady pozostają w katalogu `Lang/Examples/`.
Reguły składni są w pliku `Lang/jabba.cf`.

Na koniec przykład programu w języku Jabba:

```kotlin

fun foo_123(val str1: String, var num: new Integer) : Integer {
    writeStr(str1);

    // zwiększy kopię x o 1 i go zwróci
    num = num + 1;
    return num;
}

fun bar_abc(var num: Integer) : Unit {

    // zwiększy x dwukrotnie i zmieni oryginał
    num = num * 2;
}

fun main() : Unit {
    // x = -1
    var x: Integer = 10 - 11;

    // na stdout wpisze "abc", zwróci 0, zmienna x będzie miała wartość -1
    foo_123("abc", x);

    // po wykonaniu funkcji x będzie równa -2
    bar_abc(x);
}

// Wywołanie funkcji main (jak w Pythonie)
main();

```

### Testowanie

Wszystkie testy można uruchomić za pomocą skryptu `./check_examples.sh`, który:
- kompiluje Interpreter poleceniem `make Interpreter`
- Uruchamia Interpreter na wszystkich plikach w katalogu `Lang/Examples/`
- Uruchamia Interpreter na wszystkich plikach w katalogu `good/`
- Uruchamia Interpreter na wszystkich plikach w katalogu `bad/` i sprawdza czy zwrócił błąd
- Kompiluje testy type checker'a poleceniem `make TypeCheckerTest`
- Uruchamia testy type checker'a poleceniem `./TypeCheckerTest`


### Wymagania i kompilacja

Wersja kompilatora:
```
➜ ghc --version
The Glorious Glasgow Haskell Compilation System, version 8.10.7
```

Kompilacja:
```
make
```
Zostanie wygenerowany plik `Interpreter` oraz `interpreter`.

Poniższe polecenie uruchamia `check_examples.sh`:
```
make test
```



## Tabelka deklaracji
```txt
  Na 15 punktów
  01 (trzy typy)
  02 (literały, arytmetyka, porównania)
  03 (zmienne, przypisanie)
  04 (print)
  05 (while, if)
  06 (funkcje lub procedury, rekurencja)
  07 (przez zmienną / przez wartość)
  08 (zmienne read-only i pętla for)
  Na 20 punktów
  09 (przesłanianie i statyczne wiązanie)
  10 (obsługa błędów wykonania)
  11 (funkcje zwracające wartość)
  Na 30 punktów
  12 (4) (statyczne typowanie)
  13 (2) (funkcje zagnieżdżone ze statycznym wiązaniem)
  16 (1) (break, continue)
  17 (4) (funkcje wyższego rzędu, anonimowe, domknięcia)
  ~~18 (3) (generatory)~~

Razem: min(30, \sigma = 31) = 30
```