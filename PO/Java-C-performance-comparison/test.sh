# Script testing speed of a few programming language.
# author: Aleksander 'tudny' Tudruj at UW 
# This case Java vs C
# Usage: ./path/to/script/test.sh <temporary_test_name> <test_seed>  
# First of all it generates test with name <temporary_test_name> with
# seed <test_seed>, than compiles and runs each language one by one
# providing generated test.

# Bash Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
CLEAR='\033[0m'

# Checking if terminal supports colors.
if [ -x /usr/bin/tput ] && tput setaf 1 >&/dev/null; then
    color_prompt=yes
else
    color_prompt=
fi

# Checking number of agruments given to script.
if (( $# != 2 )); then
    if [ "$color_prompt" = yes ]; then
        echo -e "${RED}Wrong number of arguments.${CLEAR}"
    else
        echo -e "Wrong number of arguments."
    fi
    echo -e "It should be $0 <temporary_test_name> <test_seed>."
    exit 1
fi

test=$1
seed=$2

# Compiling all solutions and generator.
make

# Generating test.
./gen/gen $seed > "${test}".in

# Testing Java
echo -e "Java time"
time make runJava < "${test}".in > java_"${test}".out
echo -e ""

# Testing Java with -Xint flag
echo -e "Java time -Xint"
time make runJavaXint < "${test}".in > javaXint_"${test}".out
echo -e ""

# Testing C
echo -e "C time"
time make runC < "${test}".in > c_"${test}".out
echo -e ""

# Testing C O2
echo -e "C O2 time"
time make runC_O < "${test}".in > c_o_"${test}".out
echo -e ""

# Removing outputs and input
rm java_"${test}".out javaXint_"${test}".out c_"${test}".out c_o_"${test}".out
rm "${test}".in

# Removing compilation files
make clean

# END
