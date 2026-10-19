/** 
 * BT - Fast algorithms for matrix multiplication
 * 
 * @file divide_conquer.c
 * @details implementation of divide and conquer algorithm for matrix multiplication
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



static void divide_conquer(mtx_t *A, mtx_t *B, mtx_t *C, size_t stride, size_t n) {

    // threshold to stop recursion
    if(n == 1) {
        C[0] += A[0] * B[0];
        return;
    }

    // dividing matrix by 2
    size_t half = n / 2;

    // setting pointers to the matrix parts after divide
    mtx_t *A11 = A, *A12 = A + half, *A21 = A + half * stride, *A22 = A + half * stride + half;
    mtx_t *B11 = B, *B12 = B + half, *B21 = B + half * stride, *B22 = B + half * stride + half;
    mtx_t *C11 = C, *C12 = C + half, *C21 = C + half * stride, *C22 = C + half * stride + half;

    // C11 = A11 * B11 + A12 * B21
    divide_conquer(A11, B11, C11, stride, half);
    divide_conquer(A12, B21, C11, stride, half);

    // C12 = A11 * B12 + A12 * B22
    divide_conquer(A11, B12, C12, stride, half);
    divide_conquer(A12, B22, C12, stride, half);

    // C21 = A21 * B11 + A22 * B21
    divide_conquer(A21, B11, C21, stride, half);
    divide_conquer(A22, B21, C21, stride, half);

    // C22 = A21 * B12 + A22 * B22
    divide_conquer(A21, B12, C22, stride, half);
    divide_conquer(A22, B22, C22, stride, half);
}


int main() {

    fprintf(stdout, "\n[INFO] \t Matrix %ldx%ld -> Divide and conquer\n", (size_t)N, (size_t)N);

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

    divide_conquer(A, B, C, N, N);

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