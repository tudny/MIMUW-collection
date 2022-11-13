/** @file
 * Implementacja modułu wczytywania wejścia linia po linii.
 *
 * @author Aleksander Tudruj
 * @date 17.05.2021
 * */

#include <stdio.h>
#include <stdbool.h>
#include "reader.h"
#include "memory.h"

// na końcu jest '\0'
static char *whiteChars = "\t\n\v\f\r"; ///< Tablica białych znaków.
static size_t whiteCharsNumber = 6; ///< Liczba białych znaków.

/**
 * Sprawdzenie czy znak jest znakiem białym niebędącym spacją.
 * @param[in] c : sprawdzany znak
 * @return czy c jest znakiem białym
 * */
static bool isWhiteCharacter(char c) {
    for (size_t i = 0; i < whiteCharsNumber; ++i)
        if (c == whiteChars[i])
            return true;

    return false;
}

/**
 * Dwukrotne rozszerzenie stringa [string] o długości [actualSize].
 * @param[in,out] str : rozszerzany ciąg znaków
 * @param[in,out] actualSize : rozmiar ciągu znaków
 * */
static void extendString(char **str, size_t *actualSize) {
    size_t newSize = *actualSize << 1;

    *str = safeRealloc(*str, sizeof(char) * newSize);
    *actualSize = newSize;
}
/**
 * Dodatnie znaku [c] do stringa pod wskaźnikiem [str]. Zmiana rozmiaru
 * o ile to konieczne. Pominięcie jeżeli pierwszy znak to '#', gdyż wówczas
 * linia jest komentarzem.
 * @param[in,out] str : string do którego dodajemy
 * @param[in,out] id : indeks do którego wstawiamy
 * @param[in,out] size : wskaźnik rozmiar stringa
 * @param[in,out] c : wstawiany znak
 * @param[in,out] force : czy powinniśmy wstawić niezależnie od komentarza
 * */
static void addChar(char **str, size_t *id, size_t *size, char c, bool force) {
    if (*id + 1 == *size)
        extendString(str, size);

    if ((*str)[0] != '#' || force)
        (*str)[(*id)++] = c;
}

bool readLine(char **line, size_t *len, size_t *lastLineLength) {
    char c;
    bool isNotEOF;
    size_t id = 0, size = *lastLineLength;

    (*line)[0] = 0;
    while ((isNotEOF = (scanf("%c", &c) != EOF)) && c != EOL) {
        if (isWhiteCharacter(c)) c = '$'; // wstawienie błędnego znaku
        addChar(line, &id, &size, c, false);
    }

    if (isNotEOF || id != 0) {
        addChar(line, &id, &size, '\0', true);
        *len = id - 1;
        *lastLineLength = size;
        return true;
    }
    else {
        return false;
    }
}
