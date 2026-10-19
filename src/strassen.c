/** 
 * BT - Fast algorithms for matrix multiplication
 * 
 * @file strassen.c
 * @details implementation of Strassen algorithm for matrix multiplication
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

#define NAIVE_START 32

typedef double mtx_t;



/*static void strassen(mtx_t *A, size_t stride_A, mtx_t *B, size_t stride_B,
                     mtx_t *C, size_t stride_C, size_t stride, size_t n) {

    // 
    if(n <= NAIVE_START) {
        for(size_t idx = 0; idx < n; idx++) {
            for(size_t kdx = 0; kdx < n; kdx++) {
                mtx_t reg = A[idx * stride_A + kdx];
                for(size_t jdx = 0; jdx < n; jdx++) {
                    C[idx * stride_C + jdx] += reg * B[kdx * stride_B + jdx];
                }
            }
        }
        return;
    }


    size_t half = n / 2;

    mtx_t *A11 = A, *A12 = A + half, *A21 = A + half * stride_A, *A22 = A + half * stride_A + half;
    mtx_t *B11 = B, *B12 = B + half, *B21 = B + half * stride_B, *B22 = B + half * stride_B + half;
    

    // multiplication TODO


    
    mtx_t *C11 = C, *C12 = C + half, *C21 = C + half * stride_C, *C22 = C + half * stride_C + half;

    for(size_t idx = 0; idx < half; idx++) {
        for(size_t jdx = 0; jdx < half; jdx++) {
            size_t pos = idx * half + jdx;

            C11[idx * stride_C + jdx] += I[pos] + IV[pos] - V[pos] + VII[pos];
            C12[idx * stride_C + jdx] += III[pos] + V[pos];
            C21[idx * stride_C + jdx] += II[pos] + IV[pos];
            C22[idx * stride_C + jdx] += I[pos] + III[pos] - II[pos] + VI[pos];
        }
    }

}

void add(mtx_t *X, size_t stride_X, mtx_t *Y, size_t stride_Y,
         mtx_t *Z, size_t stride_Z, size_t n) {
    
    for(size_t idx = 0; idx < n; idx++) {
        for(size_t jdx = 0; jdx < n; jdx++) {
            Z[idx * stride_Z + jdx] = X[idx * stride_X + jdx] + Y[idx * stride_Y + jdx];
        }
    }
}

void sub(mtx_t *X, size_t stride_X, mtx_t *Y, size_t stride_Y,
         mtx_t *Z, size_t stride_Z, size_t n) {
    
    for(size_t idx = 0; idx < n; idx++) {
        for(size_t jdx = 0; jdx < n; jdx++) {
            Z[idx * stride_Z + jdx] = X[idx * stride_X + jdx] - Y[idx * stride_Y + jdx];
        }
    }
}
*/


int main() {

    fprintf(stdout, "\n[INFO] \t Matrix %ldx%ld -> Strassen\n", (size_t)N, (size_t)N);

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

    //strassen(A, B, C, N, N);

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