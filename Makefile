#	BT - Fast algorithms for matrix multiplication
#	Makefile
#	Author: xpetkob00 Branislav Peťko

CC=gcc
CFLAGS=-std=c17 -pedantic -Wall -Wextra -O3

TARGETS = imtx-multiplication fmtx-multiplication dmtx-multiplication
SRC = $(wildcard src/*.c)

OBJ_INT = $(SRC:.c=-i.o)
OBJ_FLOAT = $(SRC:.c=-f.o)
OBJ_DOUBLE = $(SRC:.c=-d.o)


.PHONY: main debug test clean

main: $(TARGETS)

imtx-multiplication: $(OBJ_INT)
	$(CC) $(CFLAGS) -DTYPE_INT $^ -o $@
fmtx-multiplication: $(OBJ_FLOAT)
	$(CC) $(CFLAGS) -DTYPE_FLOAT $^ -o $@
dmtx-multiplication: $(OBJ_DOUBLE)
	$(CC) $(CFLAGS) -DTYPE_DOUBLE $^ -o $@

%-i.o: %.c
	$(CC) $(CFLAGS) -DTYPE_INT -c $< -o $@
%-f.o: %.c
	$(CC) $(CFLAGS) -DTYPE_FLOAT -c $< -o $@
%-d.o: %.c
	$(CC) $(CFLAGS) -DTYPE_DOUBLE -c $< -o $@


debug: CFLAGS += -DDEBUG
debug: clean main


test: main


clean:
	rm -f src/*.o $(TARGETS)
