/** 
 * BT - Fast algorithms for matrix multiplication
 * 
 * @file naive_opt.c
 * @details implementation of optimized naive algorithm for matrix multiplication
 * 
 * @author xpetkob00 Branislav Peťko
*/

#define _POSIX_C_SOURCE 199309L

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#ifndef N
#define N 4096
#endif

#ifdef DEBUG
    #define debug_print(matrix, n)                          \
        do{                                                 \
            printf("----------< DEBUG >----------\n");      \
            printf("%s:\n", #matrix);                       \
            for(size_t idx = 0; idx < n; idx++) {         \
                for(size_t jdx = 0; jdx < n; jdx++) {     \
                    printf("%f, ", matrix[idx * n + jdx]);  \
                }                                           \
                printf("\n");                               \
            }                                               \
            printf("----------< DEBUG >----------\n");      \
            printf("\n\n");                                 \
        }while(0)
#else
    #define debug_print(matrix, n)
#endif

typedef double mtx_t;



static void naive_opt(mtx_t *A, mtx_t *B, mtx_t *C) {

    for(size_t idx = 0; idx < N; idx++) {
        for(size_t kdx = 0; kdx < N; kdx++) {
            mtx_t reg = A[idx * N + kdx];
            for(size_t jdx = 0; jdx < N; jdx++) {
                C[idx * N + jdx] += reg * B[kdx * N + jdx];
            }
        }
    }
}


int main() {

    fprintf(stdout, "\n[INFO] \t Matrix %ldx%ld -> Optimized naive\n", (size_t)N, (size_t)N);

    // allocating space for matrices
    mtx_t *A = calloc((size_t)N * N, sizeof(mtx_t));
    mtx_t *B = calloc((size_t)N * N, sizeof(mtx_t));
    mtx_t *C = calloc((size_t)N * N, sizeof(mtx_t));
    if(A == NULL || B == NULL || C == NULL) {
        fprintf(stderr, "Memory allocation failed!");
        return 1;
    }

    srand(42);
    // filling matrices with random numbers
    for(size_t idx = 0; idx < N * N; idx++) {
        A[idx] = rand() / (mtx_t)RAND_MAX - 0.5;
        B[idx] = rand() / (mtx_t)RAND_MAX - 0.5;
    }
    debug_print(A, N);
    debug_print(B, N);


    struct timespec start, end;
    clock_gettime(CLOCK_MONOTONIC, &start);

    naive_opt(A, B, C);

    clock_gettime(CLOCK_MONOTONIC, &end);
    unsigned long ns = (end.tv_sec - start.tv_sec) * 1000000000.0 + (end.tv_nsec - start.tv_nsec);
    double sec = ns / 1000000000.0;
    unsigned long us = (ns + 500) / 1000;   // microseconds

    printf("Total time: \t %fs  |  microsecs/element: %f  |  microsecs/mul: %f\n\n",
           sec, us / ((double)N * N), us / ((double)N * N * N));
    debug_print(C, N);


    free(A);
    free(B);
    free(C);

    return 0;
}