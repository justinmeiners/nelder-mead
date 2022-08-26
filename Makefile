
FLAGS=-std=c99 -Wall -Wextra -pedantic
CFLAGS=-O3 -ffast-math -fno-common
BIN=nm_test
SRC=test.c

$(BIN): $(SRC)
	gcc $(FLAGS) $(CFLAGS) -o $@ $< -lm

test: $(BIN)
	time ./$(BIN)
