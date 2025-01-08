
.PHONY: all clean

hello: hello.c
	mpicc -o hello hello.c

run: hello
	mpirun -np 10 hello

all: hello

clean:
	rm -f hello


benchmark: benchmark.c
	mpicc -o benchmark benchmark.c

run-benchmark: benchmark
	mpiexec -n 2 benchmark
