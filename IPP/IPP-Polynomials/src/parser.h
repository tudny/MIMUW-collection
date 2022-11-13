/** @file
 * Interfejs modułu parsowania linii.
 *
 * @author Aleksander Tudruj
 * @date 17.05.2021
 * */

#ifndef PARSER_PARSER_H
#define PARSER_PARSER_H

#include <stdbool.h>
#include <limits.h>
#include "poly.h"

/** Minimalna wartość współczynnika wieloimianu */
static const long long MIN_COEFF = LONG_MIN;
/** Maksymalna wartość współczynnika wieloimianu */
static const long long MAX_COEFF = LONG_MAX;

/** Minimalna wartość wykładnika w jednomianie */
static const unsigned long long MIN_EXP = 0;
/** Maksymalna wartość wykładnika w jednomianie */
static const unsigned long long MAX_EXP = 2147483647;

/** Minimalna wartość stopnia wieloimianu */
static const unsigned long long MIN_DEG = 0;
/** Maksymalna wartość stopnia wieloimianu */
static const unsigned long long MAX_DEG = 18446744073709551615ULL;

/** Minimalna wartość argumentu COMPOSE */
static const unsigned long long MIN_COMP = 0;
/** Maksymalna wartość argumentu COMPOSE */
static const unsigned long long MAX_COMP = 18446744073709551615ULL;

/** Typ enumeracyjny przydatny przy wczytywaniu */
typedef enum {
    LONG,
    ULONG
} NumberType;

/**
 * Sprawdzenie czy ciąg znaków może być wielomianem.
 * Sprawdzenie dotyczny całej linii.
 * @param[in] str : sprawdzany ciąg znaków
 * @param[out] p : wskaźnik na miejsce, do którego zostanie wpisany wczytany
 * wielomian, o ile takowy będzie na wyjściu. W przypadku braku wielomianu
 * w str będą tam śmieci, zatem zawartość nie jest zdefiniowana.
 * @return czy str może być wielomianem
 * */
bool CanBePoly(char *str, Poly *p);

/**
 * Sprawdzenie czy znak jest w zakresie.
 * Sprawdzenie czy znak c jest w zakresie [begin, end] z porządku ASCII.
 * @param[in] begin : początek przedziału
 * @param[in] end : koniec przedziału
 * @param[in] c : sprawdzany znak
 * */
bool isInRange(char begin, char end, char c);

/**
 * Sprawdzenie czy w str znajduje się coś mogącego być stopniem wielomianu.
 * W przypadku poprawnego wczytania czegoś będącego stopniem wielomianu
 * do zmiennej deg zostanie wpisana ta wartość, a wskaźnik endPtr zostanie
 * ustawiony na pierwszy znak nie będący stopniem wielomianu.
 * @param[in] str : sprawdzany ciąg znaków
 * @param[out] deg : wynik wczytania stopnia wielomianu
 * @param[out] endPtr : pierwszy znak za stopniem wielomianu
 * @return czy w str znajduje się stopień wielomianu
 * */
bool canBeDeg(char *str, size_t *deg, char **endPtr);

/**
 * Sprawdzenie czy w str znajduje się coś mogącego być współczynnikiem
 * wielomianu.
 * W przypadku poprawnego wczytania czegoś będącego współczynnikiem wielomianu
 * do zmiennej number zostanie wpisana ta wartość, a wskaźnik endPtr zostanie
 * ustawiony na pierwszy znak nie będący współczynnikiem wielomianu.
 * @param[in] str : sprawdzany ciąg znaków
 * @param[out] number : wynik wczytania współczynnika wielomianu
 * @param[out] endPtr : pierwszy znak za współczynnikiem wielomianu
 * @return czy w str znajduje się współczynnik wielomianu
 * */
bool canBeCoeff(char *str, poly_coeff_t *number, char **endPtr);

/**
 * Sprawdzenie czy w str znajduje się coś mogącego być arggumentem
 * polecenie COMPOSE.
 * W przypadku poprawnego wczytania czegoś będącego argumentem
 * do zmiennej number zostanie wpisana ta wartość, a wskaźnik endPtr zostanie
 * ustawiony na pierwszy znak nie będący współczynnikiem wielomianu.
 * @param[in] str : sprawdzany ciąg znaków
 * @param[out] comp : wynik wczytania argumentu
 * @param[out] endPtr : pierwszy znak za argumentem
 * @return czy w str znajduje się argument COMPOSE
 * */
bool canBeComp(char *str, size_t *comp, char **endPtr);

#endif //PARSER_PARSER_H
