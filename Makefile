#	BT - Fast algorithms for matrix multiplication
#	Makefile
#	Author: xpetkob00 Branislav Peťko

CC=gcc
CFLAGS=-std=c17 -pedantic -Wall -Wextra -O3 -march=native
N = 512
NAIVE_START = 1



.PHONY: main debug test test2 clean gzip


main: 
	$(CC) $(CFLAGS) -DN=$(N) -DNAIVE_START=$(NAIVE_START) src/naive.c -o naive
	$(CC) $(CFLAGS) -DN=$(N) -DNAIVE_START=$(NAIVE_START) src/block.c -o block
	$(CC) $(CFLAGS) -DN=$(N) -DNAIVE_START=$(NAIVE_START) src/strassen.c -o strassen


debug: CFLAGS += -DDEBUG
debug: clean main


test:
	@for algorithm in naive block strassen; do 														\
		n=2; 																						\
		for idx in $$(seq 1 12); do 																\
			$(CC) $(CFLAGS) -DN=$$n -DNAIVE_START=$(NAIVE_START) src/$$algorithm.c -o $$algorithm; 	\
			printf "%-16s N=%-5s  " $$algorithm $$n >&2; 											\
			./$$algorithm >> out.txt; 																\
			n=$$((n * 2)); 																			\
		done; 																						\
	done
	$(MAKE) clean


cutoff_test:
	@for algorithm in block strassen; do 															\
		naive_start=2; 																				\
		for idx in $$(seq 1 10); do 																\
			$(CC) $(CFLAGS) -DN=$(N) -DNAIVE_START=$$naive_start src/$$algorithm.c -o $$algorithm; 	\
			printf "%-16s CUTOFF=%-5s  " $$algorithm $$naive_start >&2; 							\
			./$$algorithm >> out.txt; 																\
			naive_start=$$((naive_start * 2)); 														\
		done; 																						\
	done
	$(MAKE) clean


clean:
	rm -f naive block strassen graphs.zip


gzip: clean
	zip -r graphs.zip graphs/
