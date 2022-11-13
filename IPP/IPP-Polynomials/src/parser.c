/** @file
 * Implementacja modułu parsowania linii.
 *
 * @author Aleksander Tudruj
 * @date 17.05.2021
 * */

#include <errno.h>
#include <stdlib.h>
#include <string.h>
#include "parser.h"
#include "memory.h"

#define INIT_MONO_SIZE 4 ///< podstawowy rozmiar tablicy jednomianów

/**
 * Sprawdzenie czy znak c jest cyfrą.
 * @param[in] c : sprawdzany znak
 * @return czy c jest cyfrą dziesiętną
 * */
static bool isDigit(char c) {
    return isInRange('0', '9', c);
}

/**
 * Sprawdzenie czy znak pod wskaźnikiem c jest znakiem pattern.
 * @param[in] c : wskaźnik na sprawdzany znak
 * @param[in] pattern : oczekiwany wskaźnik
 * @return *c == pattern
 * */
static bool is(const char * const c, char pattern) {
    return *c == pattern;
}

/**
 * Sprawdzenie czy ciąg zawiera poprawne nawiasowanie.
 * Nawiasowanie sprawdzane jest dla nawiasów '(' oraz ')'.
 * @param[in] str : sprawdzany ciąg znaków
 * @return czy str zwiera poprawne nawiasowanie
 * */
static bool hasPropperBrackets(char *str) {
    long bracketValue = 0;
    for  (; *str != '\0' && bracketValue >= 0; ++str) {
        if (*str == '(')
            bracketValue++;
        else if (*str == ')')
            bracketValue--;
    }

    return (bracketValue == 0);
}

/**
 * Sprawdzenie czy linia nie paru błędów.
 * Sprawdzenie czy linia nie zawiera oczywistych błędów, które dyskwalifikują
 * ją z bycia poprawną.
 * @param str : sprawdzana linia
 * @return czy linia zawiera oczywiste błędy
 * */
static bool hasTypos(char *str) {
    size_t len = strlen(str);

    for (size_t i = 1; i + 1 < len; ++i) {
        if (str[i] == '+' && (str[i - 1] != ')' || str[i + 1] != '('))
            return true;
    }

    for (size_t i = 1; i < len; ++i) {
        if (str[i - 1] == '-' && str[i] == '-')
            return true;
    }

    for (size_t i = 0; i < len; ++i)
        if (!isDigit(str[i]))
            if (str[i] != ')' && str[i] != '(')
                if (str[i] != '+' && str[i] != '-' && str[i] != ',')
                    return true;

    return false;
}

/**
 * Sprawdzenie czy ciąg znaków może być liczbą.
 * Obsługiwane są dwa typy liczb: 64-bitowe ze znakiem i bez.
 * @param[in] str : sprawdzany ciąg znaków
 * @param[out] number : wskaźnik na miejsce, do którego zostanie wpisana liczba
 * @param[out] endPtr : wskaźnik na miejsce, w którym zostanie zakończone
 * wpisywanie, o ile będzie ono poprawne; *endPtr = str jeżeli nie ma liczby
 * @param[in] numberType : wczytywany typ
 * @return czy w ciągu znaków str znajduje się liczba danego typu
 * */
static bool canBeNumber(char *str,
                        void *number,
                        char **endPtr,
                        NumberType numberType) {
    *endPtr = str;

    if (*str == '\0') {
        return false;
    }

    char *strPtr = str;
    *((long long *) number) = 0; // LL i LLU zachowają się tak samo
    errno = 0;

    if (is(strPtr, '-')) {
        if (numberType == ULONG || !isDigit(*(strPtr + 1))) {
            return false;
        }
    }

    if (isDigit(*strPtr) || is(strPtr, '-')) {
        switch (numberType) {
            case LONG:
                *((long long *) number) =
                        strtoll(strPtr, &strPtr, 10);
                break;
            case ULONG:
                *((unsigned long long *) number) =
                        strtoull(strPtr, &strPtr, 10);
                break;
            default:
                break;
        }

        if (strPtr != str && errno == 0) {
            *endPtr = strPtr;
            return true;
        }
    }

    return false;
}

/**
 * Sprawdzenie czy w str znajduje się coś mogącego być wykładnikiem jednomianu.
 * W przypadku poprawnego wczytania czegoś będącego wykładnikiem jednomianu
 * do zmiennej deg zostanie wpisana ta wartość, a wskaźnik endPtr zostanie
 * ustawiony na pierwszy znak nie będący wykładnikiem jednomianu.
 * @param[in] str : sprawdzany ciąg znaków
 * @param[out] number : wynik wczytania wykładnika jednomianu
 * @param[out] endPtr : pierwszy znak za wykładnikiem jednomianu
 * @return czy w str znajduje się wykłądnik jednomianu
 * */
static bool canBeExp(char *str, poly_exp_t *number, char **endPtr) {
    unsigned long long tempNumber;
    bool toReturn = canBeNumber(str, &tempNumber, endPtr, ULONG);
    if (!toReturn)
        return toReturn;

    *number = (poly_exp_t) tempNumber;
    if (!(MIN_EXP <= tempNumber && tempNumber <= MAX_EXP)) {
        *endPtr = str;
        toReturn = false;
    }

    return toReturn;
}

static bool canBeMono(char *str, Mono *m, char **endPtr);

/**
 * Dodanie jednomianu do tablicy jednomianów.
 * Ewentualne rozszerzenie tablicy jeżeli to konieczne.
 * @param[in] m : dodawany jednomian
 * @param[in,out] tab : wskaźnik na tablicę, do której dodajemy
 * @param[in,out] elems : wskaźnik na liczbę elementów w talbicy
 * @param[in,out] memSize : wskaźnik na długość tablicy
 * */
static void addSinleExtend(Mono m, Mono **tab, size_t *elems, size_t *memSize) {
    if (*tab == NULL) {
        *tab = safeCalloc(INIT_MONO_SIZE, sizeof(Mono));
        *elems = 0;
        *memSize = INIT_MONO_SIZE;
    }

    if (*elems == *memSize) {
        *memSize <<= 1;
        *tab = safeRealloc(*tab, *memSize * sizeof(Mono));
    }

    (*tab)[(*elems)++] = m;
}

/**
 * Sprawdzenie czy w str znajduje się coś mogącego być wielomianem.
 * W przypadku poprawnego wczytania czegoś będącego wielomianem
 * do zmiennej p zostanie wpisana ta wartość, a wskaźnik endPtr zostanie
 * ustawiony na pierwszy znak nie będący wielomianem.
 * @param[in] str : sprawdzany ciąg znaków
 * @param[out] p : wynik wczytania wielomianu
 * @param[out] endPtr : pierwszy znak za wielomianem
 * @return czy w str znajduje się wielomian
 * */
static bool canBePoly(char *str, Poly *p, char **endPtr) {
    // sprawdzenie czy wielomian jest stały
    poly_coeff_t number;
    char *end;

    if (canBeCoeff(str, &number, &end)) {
        *endPtr = end;
        *p = PolyFromCoeff(number);
        return true;
    }

    // nie stały
    Mono *monos = NULL;
    char *strPtr = str;
    size_t monosCnt = 0;
    size_t memorySize = 0;
    Mono tempM;
    bool lastMonoCreated = true;

    if (canBeMono(strPtr, &tempM, &end)) {
        strPtr = end;
        addSinleExtend(tempM, &monos, &monosCnt, &memorySize);

        while (*strPtr != '\0' && is(strPtr, '+') &&
               (lastMonoCreated = canBeMono(strPtr + 1, &tempM, &end))) {
            addSinleExtend(tempM, &monos, &monosCnt, &memorySize);
            strPtr = end;
        }

        if (lastMonoCreated) {
            *p = PolyAddMonos(monosCnt, monos);
            safeFree((void **) &monos);
            *endPtr = end;
            return true;
        }
    }

    for (size_t i = 0; i < monosCnt; ++i)
        MonoDestroy(&monos[i]);

    safeFree((void **) &monos);
    *endPtr = str;
    return false;
}

/**
 * Sprawdzenie czy w str znajduje się coś mogącego być jednomianem.
 * W przypadku poprawnego wczytania czegoś będącego jednomianem
 * do zmiennej m zostanie wpisana ta wartość, a wskaźnik endPtr zostanie
 * ustawiony na pierwszy znak nie będący jednomianem.
 * @param[in] str : sprawdzany ciąg znaków
 * @param[out] m : wynik wczytania jednomianu
 * @param[out] endPtr : pierwszy znak za jednomianem
 * @return czy w str znajduje się jednomian
 * */
static bool canBeMono(char *str, Mono *m, char **endPtr) {
    char *strPtr = str;
    poly_exp_t number;
    Poly tempP = PolyZero();

    if (is(strPtr, '(')) {
        if (canBePoly(str + 1, &tempP, &strPtr)) {
            if (is(strPtr, ',')) {
                if (canBeExp(strPtr + 1, &number, &strPtr)) {
                    if (is(strPtr, ')')) {
                        if (PolyIsZero(&tempP)) number = 0;
                        *m = MonoFromPoly(&tempP, number);
                        *endPtr = strPtr + 1;
                        return true;
                    }
                }
            }
        }
    }

    *endPtr = str;
    PolyDestroy(&tempP);
    return false;
}

bool canBeDeg(char *str, size_t *deg, char **endPtr) {
    unsigned long long tempNumber;
    bool toReturn = canBeNumber(str, &tempNumber, endPtr, ULONG);
    if (!toReturn)
        return toReturn;

    *deg = (size_t) tempNumber;
    if (!(MIN_DEG <= tempNumber && tempNumber <= MAX_DEG)) {
        *endPtr = str;
        toReturn = false;
    }

    return toReturn;
}

bool canBeComp(char *str, size_t *comp, char **endPtr) {
    return canBeDeg(str, comp, endPtr);
}

bool canBeCoeff(char *str, poly_coeff_t *number, char **endPtr) {
    long long tempNumber;
    bool toReturn = canBeNumber(str, &tempNumber, endPtr, LONG);
    if (!toReturn)
        return toReturn;

    *number = (poly_coeff_t) tempNumber;
    if (!(MIN_COEFF <= tempNumber && tempNumber <= MAX_COEFF)) {
        *endPtr = str;
        toReturn = false;
    }

    return toReturn;
}

bool CanBePoly(char *str, Poly *p) {

    if (!hasPropperBrackets(str))
        return false;

    if (hasTypos(str))
        return false;

    char *end;
    bool canBe = canBePoly(str, p, &end);

    if (canBe && *end != '\0') {
        PolyDestroy(p);
        return false;
    }

    return canBe && *end == '\0';
}

bool isInRange(const char begin, const char end, const char c) {
    return begin <= c && c <= end;
}
