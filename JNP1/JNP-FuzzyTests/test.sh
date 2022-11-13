#!/bin/bash

COMPILER=g++-10
if [ "$HOSTNAME" = "students" ]; then
    echo -e "Detected running script on students. Settings g++ as default compiler"
    COMPILER=g++
else
    echo -e "students machine not detected. Settings g++-10 as default compiler"
fi

CORE='\033['
RED=${CORE}'0;31m'
GREEN=${CORE}'0;32m'
CLEAR=${CORE}'0m'
BLUE=${CORE}'0;34m'

errors=0
correct=0

function test_compilable {
    directory=$1

    if ${COMPILER} -Wall -Wextra -O2 -std=c++20 -lm "${directory}"/*.cc fuzzy.cc; then
        echo -e "${GREEN}Compiled${CLEAR}"

        if ./a.out; then
            echo -e "${GREEN}Finished successfully${CLEAR}"
            (( correct = correct + 1 ))
        else
            echo -e "${RED}Finished with exit code $?${CLEAR}"
            (( errors = errors + 1 ))
        fi

        rm a.out
    else
        echo -e "${RED}Test didn't compile.${CLEAR}"
        (( errors = errors + 1 ))
    fi
}

function test_uncompilable {
    directory=$1

    if ${COMPILER} -Wall -Wextra -O2 -std=c++20 -lm "${directory}"/*.cc fuzzy.cc 2> /dev/null; then
        echo -e "${RED}Compiled and shouldn't.${CLEAR}"
        rm a.out
        (( errors = errors + 1 ))
    else
        echo -e "${GREEN}Test didn't compile and shouldn't.${CLEAR}"
        (( correct = correct + 1 ))
    fi
}

# Running the tester

function test_single {
    directory="$1"
    echo -e "=========================================="
    echo -e "${BLUE}Testing $directory${CLEAR}"

    if [[ $directory == "bad_"* ]]; then
        test_uncompilable "${directory}"
    else
        test_compilable "${directory}"
    fi

    echo -e "${BLUE}End${CLEAR}"
}

function test_all {

    echo -e "+============================+"
    echo -e "| Testing Fuzzy Numbers Task |"
    echo -e "+============================+"

    for directory in */ ; do
        test_single "$directory"
    done

    (( ratio = correct * 100 / (correct + errors) ))

    echo -e "+============================+"
    echo -e "| Raport:"
    echo -e "| ${RED}errors  = ${errors}${CLEAR}"
    echo -e "| ${GREEN}correct = ${correct}${CLEAR}"
    echo -e "| ratio   = ${ratio}%"
    echo -e "+============================+"
}

# Checking args

if [ -z $1 ]; then
    echo -e "No path to fuzzy.h (fuzzy.cc) directory"
    exit 1
fi

# Linking task solutions

if [ ! -f fuzzy.h ]; then
    ln -s $1/fuzzy.h .
fi

if [ ! -f fuzzy.cc ]; then
    ln -s $1/fuzzy.cc .
fi

# Testing all

if [ -z $2 ]; then
    echo -e "No test selected. Testing all."
    test_all
else
    test_single "$2"
fi
