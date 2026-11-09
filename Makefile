#	BT - Fast algorithms for matrix multiplication
#	Makefile
#	Author: xpetkob00 Branislav Peťko

CC=gcc
CFLAGS=-std=c17 -pedantic -Wall -Wextra -O3 -march=native
N = 4096



.PHONY: main debug test clean gzip


main: 
	$(CC) $(CFLAGS) -DN=$(N) src/naive.c -o naive
	$(CC) $(CFLAGS) -DN=$(N) src/block.c -o block
	$(CC) $(CFLAGS) -DN=$(N) src/strassen.c -o strassen


debug: CFLAGS += -DDEBUG
debug: clean main


test:
	@n=2; \
	for idx in $$(seq 1 12); do \
		for algorithm in naive block strassen; do \
			$(CC) $(CFLAGS) -DN=$$n src/$$algorithm.c -o $$algorithm; \
			printf "%-16s N=%-5s  " $$algorithm $$n >&2; \
			./$$algorithm >> out.txt; \
		done; \
		n=$$((n * 2)); \
	done
	$(MAKE) clean


clean:
	rm -f naive naive_opt divide_conquer strassen graphs.zip


gzip: clean
	zip -r graphs.zip graphs/
