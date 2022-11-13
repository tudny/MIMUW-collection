/** @file
 * Interfejs modułu zarządzania pamięcią.
 *
 * @author Aleksander Tudruj
 * @date 17.05.2021
 * */

#ifndef MEMORY_H
#define MEMORY_H

/**
 * Bezpieczne zamknięcie programu w przypadku wskaźnika na NULL.
 * Sprawdzenie czy wskaźnik wskazuje na NULL. Jeżeli tak, program kończy
 * się kodem błędu EXIT_FAILURE (1).
 * @param[in] ptr : sprawdzany wskaźnik
 * */
void *requireNotNull(void *ptr);

/**
 * Bezpieczna alternatywa malloc'a.
 * Funkcja działa podobnie do standardowej funkcji malloc, lecz w przypadku
 * braku pamięci program zakończony zostaje kodem błędu EXIT_FAILURE (1).
 * @param[in] memoryBlockSize : rezerwowany rozmiar w pamięci w bajtach
 * */
void *safeMalloc(size_t memoryBlockSize);

/**
 * Bezpieczna alternatywa calloc'a.
 * Funkcja działa podobnie do standardowej funkcji calloc, lecz w przypadku
 * braku pamięci program zakończony zostaje kodem błędu EXIT_FAILURE (1).
 * @param[in] numberOfElements : liczba rezerwowanych elementów
 * @param[in] sizeOfElement : rozmiar elementu w bajtach
 * */
void *safeCalloc(size_t numberOfElements, size_t sizeOfElement);

/**
 * Bezpieczna alternatywa realloc'a.
 * Funkcja działa podobnie do standardowej funkcji realloc, lecz w przypadku
 * braku pamięci program zakończony zostaje kodem błędu EXIT_FAILURE (1).
 * @param[in] ptr : wskaźnik na realloc'owaną pamięć
 * @param[in] newSize : nowy rozmiar pamięci w bajtach
 * */
void *safeRealloc(void *ptr, size_t newSize);

/**
 * Bezpieczna alternatywa free'a.
 * Funkcja działa podobnie do standardowej funkcji free, lecz po zwolnieniu
 * pamięci wskaźnik zostanie ustawiony na NULL.
 * @param[in,out] ptr : wskaźnik na czyszczony wskaźnik
 * */
void safeFree(void **ptr);

#endif //MEMORY_H
