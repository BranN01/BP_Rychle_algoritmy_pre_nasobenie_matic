#	BT - Fast algorithms for matrix multiplication
#	Makefile
#	Author: xpetkob00 Branislav Peťko

CC=gcc
CFLAGS=-std=c17 -pedantic -Wall -Wextra -O3 -march=native
N = 8



.PHONY: main debug test clean


main: 
	$(CC) $(CFLAGS) -DN=$(N) src/naive.c -o naive
	$(CC) $(CFLAGS) -DN=$(N) src/naive_opt.c -o naive_opt
	$(CC) $(CFLAGS) -DN=$(N) src/divide_conquer.c -o divide_conquer
	$(CC) $(CFLAGS) -DN=$(N) src/strassen.c -o strassen


debug: CFLAGS += -DDEBUG
debug: clean main


test:
	@n=2; \
	for idx in $$(seq 1 12); do \
		for algorithm in naive divide_conquer; do \
			$(CC) $(CFLAGS) -DN=$$n src/$$algorithm.c -o $$algorithm; \
			printf "%-16s N=%-5s  " $$algorithm $$n >&2; \
			for jdx in $$(seq 1 10); do \
				./$$algorithm >> out.txt; \
				printf "▰" >&2; \
			done; \
			printf " ✔\n" >&2; \
		done; \
		n=$$((n * 2)); \
	done
	$(MAKE) clean


clean:
	rm -f naive naive_opt divide_conquer strassen
