/** @file
 * Implementacja struktury stosu przechowującej wielomiany
 *
 * @author Aleksander Tudruj <at429630@students.mimuw.edu.pl>
 * @date 12.05.2021
 * */

#include <stdio.h>
#include "stack.h"

/**
 * Nowy element stosu.
 * Utworzenie nowego pustego elementu stosu. W przypadku braku pamięci
 * funkcja zwróci kod błędu wedle działania modułu memory.
 * @param[in] p : wielomian @f$p@f$
 * @return element stosu z wielomianem @f$p@f$
 * */
static SElem *createSElem(Poly *p) {
    SElem *sElem = safeMalloc(sizeof(SElem));
    sElem->p = *p;

    return sElem;
}

/**
 * Usunięcie z pamięci elementu stosu.
 * Usunięcie elementu może wystąpić w dwóch wariantach: głębokim i płytkim.
 * Usunięcie głębokie usuwa również wielomian w wierzchołu.
 * Usuniecie płytkie usuwa sam element.
 * @param sElem : usywany element
 * @param deep : czy usunięcie ma być głębokie
 * */
static void destorySingle(SElem *sElem, bool deep) {
    if (deep) PolyDestroy(&sElem->p);
    safeFree((void **) &sElem);
}

/**
 * Wypisanie elementu stosu na standardowe wyjście
 * Wypisanie elementu stosu na standardowe wyjście wraz zawartym wielomianem
 * oraz wywołanie rekurencyjne funkcji na kolejnych elementów.
 * @param sElem : wypisywany element stosu
 * */
static void printSElem(SElem *sElem) {
    if (sElem == NULL) {
        printf("TAIL{null}");
        return;
    }

    printf("Node{");
    PrintPolyNormalized(&sElem->p);
    printf("}-> ");

    printSElem(sElem->next);
}

/**
 * Ściągnięcie wierzchołka stosu.
 * Ściągnięcie wierzchołka stosu występuje w dwóch wariantach. Z usunięciem
 * wielomianu oraz bez.
 * @param stack : rozpatrywany stos
 * @param deep : czy ściągnięcie ma być głębokie
 * */
static void popStackDeep(Stack *stack, int deep) {
    assert(!isEmptyStack(stack));

    SElem *sElem = stack->head;
    stack->head = stack->head->next;
    --stack->size;

    destorySingle(sElem, deep);
}

Stack createEmptyStack() {
    Stack stack = {.head = NULL, .size = 0};

    return stack;
}

void destoryStack(Stack *stack) {
    SElem *elem = stack->head;
    while (elem != NULL) {
        SElem *temp = elem;
        elem = elem->next;
        destorySingle(temp, true);
    }
}

void pushStack(Stack *stack, Poly p) {
    SElem *sElem = createSElem(&p);

    sElem->next = stack->head;
    stack->head = sElem;
    ++stack->size;
}

bool isEmptyStack(Stack *stack) {
    return sizeStack(stack) == 0;
}

size_t sizeStack(Stack *stack) {
    return stack->size;
}

void popStack(Stack *stack) {
    popStackDeep(stack, true);
}

Poly topStack(Stack *stack) {
    assert(!isEmptyStack(stack));

    return stack->head->p;
}

Poly secondTopStack(Stack *stack) {
    assert(sizeStack(stack) >= 2);

    return stack->head->next->p;
}

Poly takeStack(Stack *stack) {
    Poly p = topStack(stack);
    popStackDeep(stack, false);
    return p;
}

void printStack(Stack *stack) {
    printSElem(stack->head);
    printf("\n");
}
