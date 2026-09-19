#	BT - Fast algorithms for matrix multiplication
#	Makefile
#	Author: xpetkob00 Branislav Peťko

CC=gcc
CFLAGS=-std=c17 -pedantic -Wall -Wextra
NAME = mtx-multiplication
SRC = $(wildcard src/*.c)
OBJ = $(SRC:.c=.o)


.PHONY: main debug test clean

main: $(NAME)

$(NAME): $(OBJ)
	$(CC) $(CFLAGS) $^ -o $@

debug: clean
	$(CC) $(CFLAGS) -DDEBUG $(SRC) -o $(NAME)

test: main

clean:
	rm -f src/*.o $(NAME)
