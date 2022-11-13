/** @file
 * Interfejs modułu odpowiedzialnego za sprawdzanie wejścia.
 *
 * @author Aleksander Tudruj
 * @date 17.05.2021
 * */

#ifndef INPUT_HANDLER_H
#define INPUT_HANDLER_H

#include <stdbool.h>
#include "stack.h"

/**
 * Sprawdzenie czy linia jest komentarzem.
 * @param[in] str : sprawdzana linia
 * @return czy linia jest komentarzem
 * */
bool isComment(const char *str);

/**
 * Sprawdzenie czy linia jest pusta.
 * @param[in] str : sprawdzana linia
 * @return czy linia jest pusta
 * */
bool isEmpty(const char *str);

/**
 * Wypisanie błędu.
 * Wypisanie błędu w linii lineNumber o treści error.
 * @param[in] lineNumber : numer linii
 * @param[in] error : treść błędu
 * */
void printError(size_t lineNumber, const char *error);

/**
 * Wypisanie błędu komendy.
 * Wypisanie informacji o wprowadzeniu błędnej komendy.
 * @param[in] lineNumber : numer linii
 * */
void wrongCommand(size_t lineNumber);

/**
 * Wypisanie błędu stosu.
 * Wypisanie informacji o braku wystarczającej liczby elementów na stosie.
 * @param[in] lineNumber : numer linii
 * */
void stackUnderflow(size_t lineNumber);

/**
 * Sprawdzenie linii pod kątem bycia komendą.
 * Sprawdzenie czy linia jest kandydatem do bycia poprawną komendą.
 * @param[in] str : sprawdzany ciąg znaków
 * @return czy ciąg może być komendą
 * */
bool pretendsToBeCommand(const char *str);

/**
 * Obsługa linii kandydującej do bycia komendą.
 * Wyszukany jest odpowiedni handler komendy oraz ewentualnie wypisany
 * jest błąd.
 * @param[in] str : sprawdzany ciąg znaków
 * @param[in] lineNumber : numer linii
 * @param[in,out] stack : stos kalkulatora
 * */
void handleCommand(char *str, size_t lineNumber, Stack *stack);

/**
 * Obsługa linii kandydującej do bycia wielomianem.
 * Ewentualnie wypisany odpowiedni jest błąd.
 * @param[in] str : sprawdzany ciąg znaków
 * @param[in] lineNumber : numer linii
 * @param[in,out] stack : stos kalkulatora
 * */
void handlePoly(char *str, size_t lineNumber, Stack *stack);

#endif //INPUT_HANDLER_H
