#!/bin/bash

set -e

# Go to script location
cd "$(dirname "$0")"

if [ "$#" -ne 4 ]; then
    echo "Usage: $0 <N> <M> <K> <OUTPUT_DIR>"
    exit 1
fi

N=$1
M=$2
K=$3
OUTPUT_DIR=$4

python3 gen-example-matrices.py --output "$OUTPUT_DIR" --rows "$N" --cols "$M" --dens 0.3 --num 1 --name a{}.csr
python3 gen-example-matrices.py --output "$OUTPUT_DIR" --rows "$M" --cols "$K" --dens 0.3 --num 1 --name b{}.csr
