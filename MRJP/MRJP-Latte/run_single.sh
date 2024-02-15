
# set -e

lat_file="$1".lat
ll_file="$1".ll
bc_file="$1".bc

echo "\n========================== latc_llvm =========================="

./latc_llvm $lat_file

echo "\n========================== clang =========================="

clang -o $1 $bc_file -g

echo "\n========================== valgrind =========================="

valgrind --error-exitcode=123 --leak-check=full \
		--show-leak-kinds=all --errors-for-leak-kinds=all \
		--track-origins=yes \
        ./$1

echo "\n========================== raw =========================="

./$1
