#!/bin/bash

if [ -n "$1" ]; then
    BUILD_TYPE=$1
else
    BUILD_TYPE=Debug
fi

THIS_SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
EXEC_PATH="cmake-build-${BUILD_TYPE,,}/matmul"
# ./gen-mult-ex.py --n 7 --m 7 --k 7 --dens 0.3 --seed 42 --output output-mult-sq --name "{}.csr"
ARGS="-a ${THIS_SCRIPT_DIR}/data/output-mult-sq-100/A.csr -b ${THIS_SCRIPT_DIR}/data/output-mult-sq-100/B.csr -t 2D -g 100 -v"

if [ ! -f "$EXEC_PATH" ]; then
    echo "Executable not found: $EXEC_PATH"
    exit 1
fi

echo "Running $EXEC_PATH with $ARGS"

# shellcheck disable=SC2086
mpiexec -n 4 --hostfile hostfile "$EXEC_PATH" $ARGS
