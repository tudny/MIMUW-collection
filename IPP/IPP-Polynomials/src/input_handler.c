/** @file
 * Implementacja modułu odpowiedzialnego za sprawdzanie wejścia.
 *
 * @author Aleksander Tudruj
 * @date 17.05.2021
 * */

#include <string.h>
#include <stdio.h>
#include "input_handler.h"
#include "parser.h"
#include "command_handler.h"

#define SIZE(x) ((size_t) (sizeof (x) / sizeof (x)[0])) ///< rozmiar tablicy X

/**
 * Struktura opisująca komendę i jej handler.
 * */
typedef struct {
    char *const name; ///< nazwa komendy
    void (*command)(char *const, size_t, Stack *); ///< obsługa komendy
    bool hasArgument; ///< czy polecenie potrzebuje argumentu
} Command;

/** Lista komend */
Command commands[] = {
        {"ZERO", handleZero, false},
        {"IS_COEFF", handleIsCoeff, false},
        {"IS_ZERO", handleIsZero, false},
        {"CLONE", handleClone, false},
        {"ADD", handleAdd, false},
        {"MUL", handleMul, false},
        {"NEG", handleNeg, false},
        {"SUB", handleSub, false},
        {"IS_EQ", handleIsEq, false},
        {"DEG_BY", handleDegBy, true},
        {"DEG", handleDeg, false},
        {"AT", handleAt, true},
        {"PRINT", handlePrint, false},
        {"POP", handlePop, false},
        {"COMPOSE", handleCompose, true}
};

/**
 * Sprawdzenie czy ciąg znaków zaczyna się danym ciągiem.
 * @param[in] start : oczekiwany początek
 * @param[in] str : sprawdzany ciąg znaków
 * @param[in] canBeMore : flaga dająca mażliwość czy ciąg znaków może mieć coś
 * więcej niż tylko zadany początek
 * @return czy str zaczyna się na start
 * */
static bool startsWithOrIs(char *const start, char *const str, bool canBeMore) {
    size_t startSize = strlen(start);
    size_t strLen = strlen(str);

    if (startSize > strLen)
        return false;

    for (size_t i = 0; i < startSize; ++i) {
        if (start[i] != str[i])
            return false;
    }

    return canBeMore || strLen == startSize;
}

/**
 * Wypisanie błędu o błędym wielomianie.
 * @param[in] lineNumber : numer linii
 * */
static void wrongPoly(size_t lineNumber) {
    printError(lineNumber, "WRONG POLY");
}

bool isComment(const char *const str) {
    return *str == '#';
}

bool isEmpty(const char *const str) {
    return *str == '\0';
}

void printError(size_t lineNumber, const char *const error) {
    fprintf(stderr, "ERROR %zu %s\n", lineNumber, error);
}

void wrongCommand(size_t lineNumber) {
    printError(lineNumber, "WRONG COMMAND");
}

void stackUnderflow(size_t lineNumber) {
    printError(lineNumber, "STACK UNDERFLOW");
}

bool pretendsToBeCommand(const char *const str) {
    assert(str != NULL);

    return isInRange('A', 'Z', *str) ||
           isInRange('a', 'z', *str);
}

void handleCommand(char *str, size_t lineNumber, Stack *stack) {
    for (size_t i = 0; i < SIZE(commands); ++i) {
        if (startsWithOrIs(commands[i].name, str, commands[i].hasArgument)) {
            (*commands[i].command)(str, lineNumber, stack);
            return;
        }
    }

    wrongCommand(lineNumber);
}

void handlePoly(char *str, size_t lineNumber, Stack *stack) {
    Poly p;
    if (CanBePoly(str, &p))
        pushStack(stack, p);
    else
        wrongPoly(lineNumber);
}
