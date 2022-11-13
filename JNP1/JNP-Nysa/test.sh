#!/bin/bash
# Skrypt testujący do zadania Similar Line - małego zadania z IPP.
# @author Aleksander "tudny" Tudruj

# Ustawienie RUN_MEMORY_TEST na 1 spowoduje dodatkowe sprawdzenie
# programu za pomocą valgrinda.
RUN_MEMORY_TEST=0

CORE='\033['
RED=${CORE}'0;31m'
GREEN=${CORE}'0;32m'
CLEAR=${CORE}'0m'

# Sprawdzenie czy liczba argumentów jest odpowiednia.
# Uruchamiający powinien zapewnić dokładnie dwa argumenty.
# Pierwszy z nich to plik wykonywalny.
# Drugi to ścieżka do folderu z testami.
if (($# != 2))
then
    echo -e "${RED}Error${CLEAR}\nWrong number of arguments. Should be 2."
    exit 1
fi

# Odczywanie argumentów użytkownika do zmiennych.
prog=$1
dir=$2

# Jeżeli na końcu nazwy folderu z testami nie ma separatora
# jest on dopisywany.
if [ "${dir:(-1)}" != "/" ]; then
  dir="${dir}/"
fi

# Jeżeli ścieżka do programu jest absolutna, to wywołanie odbywa sie bez "./"
if [[ "${prog}" == /* ]]; then
  exec="${prog}"
else
  exec="./${prog}"
fi

# Sprawdzenie czy plik wykonywalny istnieje.
if [ ! -f "${prog}" ]; then
  echo -e "${RED}File ${prog} doesn't exists!${CLEAR}"
  exit 1
fi

# Sprawdzenie czy folder z testami istnieje.
if [ ! -d "${dir}" ]; then
  echo -e "${RED}Dir ${dir} doesn't exists!${CLEAR}"
  exit 1
fi

# Odpowiednio licznik testów poprawnych i błędnych.
correct=0
error=0

for test in "${dir}"*.in;
do
    IN_TEST="${test}"
    OUT_TEST="${test%.in}".out
    ERR_TEST="${test%.in}".err

    OUT_MY="${test%.in}"_temp.out
    ERR_MY="${test%.in}"_temp.err

    "${exec}" < "${IN_TEST}" > "${OUT_MY}" 2> "${ERR_MY}"

    if diff -Z "${OUT_TEST}" "${OUT_MY}" > /dev/null && \
        diff -Z "${ERR_TEST}" "${ERR_MY}" > /dev/null
    then
        echo -e "${GREEN}${IN_TEST}: \tOK!${CLEAR}"
        ((correct = correct + 1))

        # VALGRIND TEST | ONLY IF: RUN_MEMORY_TEST == 1
        if [[ "${RUN_MEMORY_TEST}" == 1 ]]; then
          valgrind --error-exitcode=123 --leak-check=full \
	            --show-leak-kinds=all --errors-for-leak-kinds=all \
              "${exec}" < "${IN_TEST}" > /dev/null 2> /dev/null

          if [[ "$?" == 123 ]]; then
            echo -e "${RED}Memory error | Valgrid test failed! ${CLEAR}"
          else
            echo -e "${GREEN}Memory test passed. ${CLEAR}"
          fi
        fi
        # VALGRIND TEST END

    else
        echo -e "${RED}${IN_TEST}: \tERROR!${CLEAR}"
        ((error = error + 1))
    fi

    rm "${OUT_MY}" "${ERR_MY}"
done

sum=$((correct + error))

echo -e "===========================================";
echo -e "Number of tests: \t$((sum))"
echo -e "${GREEN}Tests passed: \t\t${correct} \t($((correct * 100 / sum))%)${CLEAR}"
echo -e "${RED}Tests not passed: \t${error} \t($((error * 100 / sum))%)${CLEAR}"
echo -e "===========================================";

exit $error
