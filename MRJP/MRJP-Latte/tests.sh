#!/bin/bash

dirs=(
    lattests/extensions/arrays1
    lattests/extensions/objects1
    lattests/extensions/objects2
    lattests/extensions/struct
    lattests/good
    lattests/bad
    ../mrjp-tests/good
    # ../mrjp-tests/bad
)

total=0
passed=0
not_passed=0

for dir in ${dirs[@]}; do
    for file in $dir/*.lat; do
        if ./latc_llvm $file > dist-newstyle/stdout.out 2>dist-newstyle/stdout.err; then

            file_wo="${file%.*}"

            input_file=$file_wo.input
            if [ -f $input_file ]; then
                lli $file_wo.bc > dist-newstyle/llvm.out 2>dist-newstyle/llvm.err < $input_file
                exit_code=$?
            else
                lli $file_wo.bc > dist-newstyle/llvm.out 2>dist-newstyle/llvm.err
                exit_code=$?
            fi

            expected=$file_wo.output

            if diff dist-newstyle/llvm.out $expected; then
                echo "Passed: $file"
                passed=$((passed+1))
            else
                echo "! Failed: $file"
                not_passed=$((not_passed+1))
            fi
        else
            dir_name=$(dirname $file)
            if [ $dir_name == "lattests/bad" ] || [ $dir_name == "../mrjp-tests/bad" ]; then
                echo "Passed: $file"
                passed=$((passed+1))
            else
                echo "! Failed: $file"
                not_passed=$((not_passed+1))
            fi
        fi
        total=$((total+1))
    done
done

echo "====================================================================="
echo "Total: $total"
echo "Passed: $passed"
echo "Not passed: $not_passed"
