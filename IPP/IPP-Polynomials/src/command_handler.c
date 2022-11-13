/** @file
 * Implementacja modułu odpowiedzialnego za obsługę komend.
 *
 * @author Aleksander Tudruj
 * @date 17.05.2021
 * */

#include <string.h>
#include <stdio.h>
#include "command_handler.h"
#include "input_handler.h"
#include "parser.h"

/**
 * Sprawdzenie czy stos zawiera odpowiednią liczbę elementów.
 * Jeżeli stos zawiera mniej niż x elementów wypisywany jest błąd.
 * @param[in] lineNumber : numer linii do wypisania błędu
 * @param[in] stack : sprawdzany stos
 * @param[in] x : oczekiwana liczba elementów
 * @return czy stos zawiera co najmniej x elementów
 * */
static bool stackHasXPolys(size_t lineNumber, Stack *stack, size_t x) {
    if (sizeStack(stack) < x) {
        stackUnderflow(lineNumber);
        return false;
    }

    return true;
}


void handleZero(__attribute__((unused)) char *const str,
                size_t lineNumber,
                Stack *stack) {
    if (!stackHasXPolys(lineNumber, stack, 0))
        return;

    pushStack(stack, PolyZero());
}

void handleIsCoeff(__attribute__((unused)) char *const str,
                   size_t lineNumber,
                   Stack *stack) {
    if (!stackHasXPolys(lineNumber, stack, 1))
        return;

    Poly top = topStack(stack);
    printf("%d\n", PolyIsCoeff(&top) ? 1 : 0);
}

void handleIsZero(__attribute__((unused)) char *const str,
                  size_t lineNumber,
                  Stack *stack) {
    if (!stackHasXPolys(lineNumber, stack, 1))
        return;

    Poly top = topStack(stack);
    printf("%d\n", PolyIsZero(&top) ? 1 : 0);
}

void handleClone(__attribute__((unused)) char *const str,
                 size_t lineNumber,
                 Stack *stack) {
    if (!stackHasXPolys(lineNumber, stack, 1))
        return;

    Poly top = topStack(stack);
    Poly cloned = PolyClone(&top);
    pushStack(stack, cloned);
}

void handleAdd(__attribute__((unused)) char *const str,
               size_t lineNumber,
               Stack *stack) {
    if (!stackHasXPolys(lineNumber, stack, 2))
        return;

    Poly a = takeStack(stack);
    Poly b = takeStack(stack);
    Poly sum = PolyAdd(&a, &b);

    PolyDestroy(&a);
    PolyDestroy(&b);

    pushStack(stack, sum);
}

void handleMul(__attribute__((unused)) char *const str,
               size_t lineNumber,
               Stack *stack) {
    if (!stackHasXPolys(lineNumber, stack, 2))
        return;

    Poly a = takeStack(stack);
    Poly b = takeStack(stack);
    Poly sum = PolyMul(&a, &b);

    PolyDestroy(&a);
    PolyDestroy(&b);

    pushStack(stack, sum);
}

void handleNeg(__attribute__((unused)) char *const str,
               size_t lineNumber,
               Stack *stack) {
    if (!stackHasXPolys(lineNumber, stack, 1))
        return;

    Poly top = takeStack(stack);
    Poly negated = PolyNegProperty(&top);

    pushStack(stack, negated);
}

void handleSub(__attribute__((unused)) char *const str,
               size_t lineNumber,
               Stack *stack) {
    if (!stackHasXPolys(lineNumber, stack, 2))
        return;

    Poly a = takeStack(stack);
    Poly b = takeStack(stack);
    Poly sum = PolySubProperty(&a, &b);

    pushStack(stack, sum);
}

void handleIsEq(__attribute__((unused)) char *const str,
                size_t lineNumber,
                Stack *stack) {
    if (!stackHasXPolys(lineNumber, stack, 2))
        return;

    Poly a = topStack(stack);
    Poly b = secondTopStack(stack);

    printf("%d\n", PolyIsEq(&a, &b) ? 1 : 0);
}

void handleDeg(__attribute__((unused)) char *const str,
               size_t lineNumber,
               Stack *stack) {
    if (!stackHasXPolys(lineNumber, stack, 1))
        return;

    Poly a = topStack(stack);

    printf("%d\n", PolyDeg(&a));
}

void handleDegBy(char *const str,
                 size_t lineNumber,
                 Stack *stack) {
    char *name = "DEG_BY";
    char *spaceAndArgument = (char *) str + strlen(name);
    size_t argument;
    char *endPtr;

    if (*spaceAndArgument != ' ' ||
            !canBeDeg(spaceAndArgument + 1, &argument, &endPtr) ||
            *endPtr != '\0') {
        printError(lineNumber, "DEG BY WRONG VARIABLE");
        return;
    }

    if (!stackHasXPolys(lineNumber, stack, 1))
        return;

    Poly a = topStack(stack);
    poly_exp_t degBy = PolyDegBy(&a, argument);
    printf("%d\n", degBy);
}

void handleAt(char *const str,
              size_t lineNumber,
              Stack *stack) {
    char *name = "AT";
    char *spaceAndArgument = ((char *) str) + strlen(name);
    poly_coeff_t argument;
    char *endPtr;

    if (*spaceAndArgument != ' ' || *(spaceAndArgument + 1) == '\0' ||
            !canBeCoeff(spaceAndArgument + 1, &argument, &endPtr) ||
            *endPtr != '\0'){
        printError(lineNumber, "AT WRONG VALUE");
        return;
    }

    if (!stackHasXPolys(lineNumber, stack, 1))
        return;

    Poly a = takeStack(stack);
    Poly at = PolyAt(&a, argument);

    PolyDestroy(&a);

    pushStack(stack, at);
}

void handlePrint(__attribute__((unused)) char *const str,
                 size_t lineNumber,
                 Stack *stack) {
    if (!stackHasXPolys(lineNumber, stack, 1))
        return;

    Poly top = topStack(stack);
    PrintPolyNormalized(&top);
    printf("\n");
}

void handlePop(__attribute__((unused)) char *const str,
               size_t lineNumber,
               Stack *stack) {
    if (!stackHasXPolys(lineNumber, stack, 1))
        return;

    popStack(stack);
}

void handleCompose(char *const str,
                   size_t lineNumber,
                   Stack *stack) {
    char *name = "COMPOSE";
    char *spaceAndArgument = (char*) str + strlen(name);
    size_t argument;
    char *endPtr;

    if (*spaceAndArgument != ' ' ||
        !canBeComp(spaceAndArgument + 1, &argument, &endPtr) ||
        *endPtr != '\0') {
        printError(lineNumber, "COMPOSE WRONG PARAMETER");
        return;
    }

    if (!stackHasXPolys(lineNumber, stack, argument + 1)
        || !stackHasXPolys(lineNumber, stack, argument))
        return;

    Poly base = takeStack(stack);
    Poly *substitutes = safeCalloc(argument, sizeof(Poly));

    for (size_t i = argument; i > 0; --i) {
        substitutes[i - 1] = takeStack(stack);
    }

    Poly res = PolyCompose(&base, argument, substitutes);
    pushStack(stack, res);

    PolyDestroy(&base);

    for (size_t i = 0; i < argument; i++) {
        PolyDestroy(&substitutes[i]);
    }

    safeFree((void **) &substitutes);
}


