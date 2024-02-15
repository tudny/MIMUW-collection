#!/bin/bash

# This is a test script.

make Interpreter > /dev/null 2>&1

for file in Lang/Examples/*.jbb; do
    echo "Testing $file"
    ./Interpreter $file > /dev/null 2>&1
    if [ $? -ne 0 ]; then
        echo "Failed to eval $file"
        exit 1
    fi
done

for file in good/*.jbb; do
    echo "Testing $file"
    ./Interpreter $file > /dev/null 2>&1
    if [ $? -ne 0 ]; then
        echo "Failed to eval $file"
        exit 1
    fi
done

for file in bad/*.jbb; do
    echo "Testing $file"
    ./Interpreter $file > /dev/null 2>&1
    if [ $? -eq 0 ]; then
        echo "Evaled $file, but should not have"
        exit 1
    fi
done

make TypeCheckerTest > /dev/null 2>&1

./TypeCheckerTest
if [ $? -ne 0 ]; then
    echo "TypeCheckerTest failed"
    exit 1
fi

echo -e "\e[32mAll tests passed\e[0m"
