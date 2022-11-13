# Source files for algorithm implementation
SRC = src
JAVA_FILE = com/github/tudny/algo/ShortestPath
C_FILE = $(SRC)/c/shortest_path

# Source of test generator file
GEN_FILE = gen/gen

.PHONY : all runJava runJavaXint runC runOcaml gen

all : $(SRC)/java/$(JAVA_FILE).class $(C_FILE) $(GEN_FILE) Makefile

# Compiling Java solution
$(SRC)/java/$(JAVA_FILE).class : $(SRC)/java/$(JAVA_FILE).java
	( cd $(SRC)/java && javac $(JAVA_FILE).java )

# Running solution in Java
runJava : $(SRC)/java/$(JAVA_FILE).class
	( cd $(SRC)/java && java $(JAVA_FILE) )

# Running solution in Java with -Xint flag
runJavaXint: $(SRC)/java/$(JAVA_FILE).class
	( cd $(SRC)/java && java $(JAVA_FILE) -Xint )

# Compiling C solution
$(C_FILE) : $(C_FILE).c
	gcc -Wall -Wextra -o $(C_FILE) $(C_FILE).c

# Running solution in C 
runC : $(C_FILE)
	./$(C_FILE)

# Compiling C O2 solution
$(C_FILE)_O : $(C_FILE).c
	gcc -Wall -Wextra -O2 -o $(C_FILE)_O $(C_FILE).c

# Running solution in C 
runC_O : $(C_FILE)_O
	./$(C_FILE)_O

# Cleaning whole directory from all Java, C and OCaml compilation files.
clean :
	@find . -name "*.class" -type f -delete
	@rm -f $(C_FILE) $(C_FILE)_O
	@rm -f gen/gen

# Compiling test generator
$(GEN_FILE) : $(GEN_FILE).c
	gcc -Wextra -Wall -o $(GEN_FILE) $(GEN_FILE).c

# Running test generator
runGen : $(GEN_FILE)
	./$(GEN_FILE)
