#!/usr/bin/env bash

mkdir build

gcc -Wall -Wextra -O2 -c -o build/test.o src/test.c
g++ -Wall -Wextra -O2 -c -o build/maptel.o src/maptel.cc

g++ -Wall -Wextra -O2 -o build/test build/test.o build/maptel.o src/maptel.h
