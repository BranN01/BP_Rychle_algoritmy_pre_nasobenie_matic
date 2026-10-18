/** 
 * BT - Fast algorithms for matrix multiplication
 * 
 * @file naive.h
 * @details implementation of naive algorithm for matrix multiplication
 * 
 * @author xpetkob00 Branislav Peťko
*/

#include "naive.h"


void naive(size_t n) {

    printf("Algorithm: \t Naive\n");

    struct timespec start, end;
    clock_gettime(CLOCK_MONOTONIC, &start);

    // naive algorithm
    for(size_t idx = 0; idx < n; idx++) {
        for(size_t jdx = 0; jdx < n; jdx++) {
            for(size_t kdx = 0; kdx < n; kdx++) {
                g.mtx_C[idx * n + jdx] += g.mtx_A[idx * n + kdx] * g.mtx_B[kdx * n + jdx];
            }
        }
    }

    clock_gettime(CLOCK_MONOTONIC, &end);
    double sec = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1000000000.0;

    printf("Total time: \t %fs  |  microsecs/element: %f  |  microsecs/mul: %f\n\n",
           sec, sec * 1000000.0 / ((double)n * n), sec * 1000000.0 / ((double)n * n * n));
    debug_print(g.mtx_C, n);
    fprintf(stderr, "========================================================================\n");

    // clearing C matrix
    memset(g.mtx_C, 0, n * n * sizeof(mtx_t));
}

void naive_optimized(size_t n) {

    printf("Algorithm: \t Optimized Naive\n");

    struct timespec start, end;
    clock_gettime(CLOCK_MONOTONIC, &start);

    // optimized naive algorithm
    for(size_t idx = 0; idx < n; idx++) {
        for(size_t kdx = 0; kdx < n; kdx++) {
            mtx_t reg = g.mtx_A[idx * n + kdx];
            for(size_t jdx = 0; jdx < n; jdx++) {
                g.mtx_C[idx * n + jdx] += reg * g.mtx_B[kdx * n + jdx];
            }
        }
    }

    clock_gettime(CLOCK_MONOTONIC, &end);
    double sec = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1000000000.0;

    printf("Total time: \t %fs  |  microsecs/element: %f  |  microsecs/mul: %f\n\n",
           sec, sec * 1000000.0 / ((double)n * n), sec * 1000000.0 / ((double)n * n * n));
    debug_print(g.mtx_C, n);
    fprintf(stderr, "========================================================================\n");

    // clearing C matrix
    memset(g.mtx_C, 0, n * n * sizeof(mtx_t));
}
