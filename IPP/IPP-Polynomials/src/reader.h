/** @file
 * Interfejs modułu wczytywania wejścia linia po linii.
 *
 * @author Aleksander Tudruj
 * @date 17.05.2021
 * */

#ifndef READER_H
#define READER_H

#include <stdbool.h>

#define EOL '\n' ///< znak końca linii

/**
 * Wczytanie linii z stdin do zmiennej [line].
 * W przypadku pustej linii na końcu pliku lub końca pliku zwracane
 * jest false. W przeciwnym przypadku zwracane jest true.
 * @param[out] line :  wskaźnik na linię do wczytania danych
 * @param[out] len : wskaźnik na zmienną z długością linii
 * @param[in,out] lastLineLength : wielkość tablicy na linię
 * @return czy została wczytana nowa linia
 * */
bool readLine(char **line, size_t *len, size_t *lastLineLength);

#endif // READER_H