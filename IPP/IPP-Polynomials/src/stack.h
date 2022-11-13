/** @file
 * Interfejs struktury stosu przechowującej wielomiany
 *
 * @author Aleksander Tudruj <at429630@students.mimuw.edu.pl>
 * @date 12.05.2021
 * */

#ifndef STACK_H
#define STACK_H

#include "poly.h"
#include "memory.h"

/** To jest pojedynczy element stosu. */
typedef struct SElem {
   struct SElem *next; ///< następny element stosu
   Poly p; ///< wielomian w wierzchołku
} SElem;

/** To jest struktura stosu */
typedef struct Stack {
    SElem *head; ///< góra stosu
    size_t size; ///< rozmiar stosu
} Stack;

/**
 * Funkcja tworzy pusty stos
 * Funkcja tworzy pusty stos i ustawia początkowe wartości. W przypadku
 * braku pamięci program kończy wedle moduły memory.
 * @return wskaźnik na nowy pusty stos
 * */
Stack createEmptyStack();

/**
 * Funkcja usuwa stos z pamięci.
 * Funkcja usuwa stos z pamięci wraz z całą zawartością.
 * @param[in] stack : usuwany stos
 * */
void destoryStack(Stack *stack);

/**
 * Dodanie wielomianu do stosu.
 * Dodanie wielomianu na szczyt stosu.
 * @param[in,out] stack : stos, do którego dodajemy
 * @param[in] p : dodawany wielomian
 * */
void pushStack(Stack *stack, Poly p);

/**
 * Sprawdzenie czy stos jest pusty.
 * Sprawdzenie czy stos nie zawiera żadnych elementów.
 * @param[in] stack : sprawdzany stos
 * @return czy stos jest pusty
 * */
bool isEmptyStack(Stack *stack);

/**
 * Rozmiar stosu.
 * Funkcja zwraca aktualny rozmiar stosu.
 * @param[in] stack : stos, którego rozmiar sprawdzamy
 * */
size_t sizeStack(Stack *stack);

/**
 * Usunięcie wierzchołka.
 * Usunięcie wierzchołka ze stosu. W przypadku stosu pustego program
 * zakończy się fałszywą asercją.
 * @param[in,out] stack : stos, z którego usuwamy wierzchołek
 * */
void popStack(Stack *stack);

/**
 * Zwrócenie wierzchołka.
 * Zwrócenie wielomianu z wierzchołka ze stosu. W przypadku stosu pustego
 * program zakończy się fałszywą asercją.
 * @param[in] stack : rozpatrywany stos
 * @return wierzchołek rozpatrywanego stosu
 * */
Poly topStack(Stack *stack);

/**
 * Zwrócenie elementu pod wierzchołkiem wierzchołka.
 * Zwrócenie wielomianu spod wierzchołka stosu. W przypadku stosu pustego
 * lub o jednym elemencie program zakończy się fałszywą asercją.
 * @param[in] stack : rozpatrywany stos
 * @return wielomian spod wierzchołka rozpatrywanego stosu
 * */
Poly secondTopStack(Stack *stack);

/**
 * Zabranie wierzchołka.
 * Zabranie wierzchołka ze stosu t.j. usunięcie do ze stosu
 * oraz zwrócenie go. Efekt identyczny jak po wywołaniu topStack(),
 * a potem popStack(), lecz w jednej funkcji.
 * @param[in,out] stack : rozpatrywany stos
 * @return wierzchołek rozpatrywanego stosu
 * */
Poly takeStack(Stack *stack);

/**
 * Wypisanie stosu.
 * Wypisanie stosu na standardowe wyjście wraz z wielomianami.
 * @param[in] stack : wypisywany stos
 * */
void printStack(Stack *stack);

#endif //STACK_H
