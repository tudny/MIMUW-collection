# Pytania do zadania 2 (maptel)

## Numery telefonów

Według treści zadania
```
...jako ciąg maksymalnie 22 cyfr zapisanych w ASCII i jest zakończony
znakiem o wartości zero. ...
```
oraz
```
Należy sprawdzać przynajmniej, czy identyfikator słownika jest
poprawny, [...] czy przekazany napis nie
jest za długi, [...] i czy kończy się znakiem o wartości
zero.
```
powinniśmy sprawdzać poprawność argumentów za pomocą asercji. W przypadku kod korzystającego z naszej biblioteki
```c
#include "maptel.h"
...
int foo() {
    unsigned long id = ... 
    ...
    char *str = malloc(sizeof(char) * 15);
    str[0] = '1';
    ...
    str[14] = '1';
    // Cały ciąg to znane nam znaki różne od '\0', lecz nie kończy się on znakiem '\0'.
    
    maptel_insert(id, str, ...);
...
}
```

Jak w takiej sytuacji powinien zachować się nasz program przy:
* 1.1) kompilacji w wersji DEBUG? Nie możemy odczytać długości tablicy za alokowanej dynamicznie za pomocą wskaźnika, a miejsce w pamięci
  `str + 15` może nie być znakiem o wartości zero, co więcej może spowodować naruszenie bezpieczeństwa przy odczycie.
* 1.2) kompilacji w wersji RELEASE? Tutaj spodziewalibyśmy się, że kod bez włączonych asercji może po prostu zakończyć się komentarzem 'core dumped'.

Pytanie jest skierowane do każdej z funkcji przyjmującej ciąg znaków jako argument.

Oczywiście, gdy ciąg znaków będzie dłuższy niż TEL_NUM_MAX_LEN to nie ma problemu. Gorzej w napisanym wyżej przykładzie.

## Przypisania wyniku

W wersji DEBUG, jak powinna zachować się funkcja `maptel_transform` w przypadku, gdy rzeczywista długość tablicy `tel_dst` będzie różna od zmiennej `len`?


