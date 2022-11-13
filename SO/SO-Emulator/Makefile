PROJECT = so_emulator

export N=4

NASM = nasm
NASMFLAGS = -DCORES=$N -f elf64 -w+all -w+error -g
OBJS = so_emulator.o so_emulator_example.o
CC = gcc
CFLAGS = -DCORES=$N -c -Wall -Wextra -std=c17 -O2 -g
LDFLAGS = -pthread

.PHONY : all clean valgrind

all : $(PROJECT)

$(PROJECT) : $(OBJS)

so_emulator.o : so_emulator.asm 
	$(NASM) $(NASMFLAGS) $< -o $@

so_emulator_example.o : so_emulator_example.c

clean :
	rm -f $(PROJECT) $(OBJS)

valgrind : $(PROJECT)
	valgrind --error-exitcode=123 --leak-check=full \
		--show-leak-kinds=all --errors-for-leak-kinds= all \\
		--track-origins=yes ./$(PROJECT)
