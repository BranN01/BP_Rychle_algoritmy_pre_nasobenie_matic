/** 
 * BT - Fast algorithms for matrix multiplication
 * 
 * @file basic.h
 * @details implementation of basic algorithm for matrix multiplication
 * 
 * @author xpetkob00 Branislav Peťko
*/

#include "basic.h"


void basic(unsigned n) {

    printf("Algorithm: \t\t Basic\n");

    struct timespec start, end;
    clock_gettime(CLOCK_MONOTONIC, &start);

    for(unsigned idx = 0; idx < n; idx++) {
        for(unsigned jdx = 0; jdx < n; jdx++) {
            g.mtx_C[idx * n + jdx] = 0;
            for(unsigned kdx = 0; kdx < n; kdx++) {
                g.mtx_C[idx * n + jdx] += g.mtx_A[idx * n + kdx] * g.mtx_B[kdx * n + jdx];
            }
        }
    }

    clock_gettime(CLOCK_MONOTONIC, &end);
    double sec = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1000000000.0;

    printf("Time to compute: \t %fs\n\n", sec);
    debug_print(g.mtx_C, n);
    fprintf(stderr, "========================================================================\n");

    // clearing C matrix
    memset(g.mtx_C, 0, n * n * sizeof(float));
}


void basic_optimized(unsigned n) {

    printf("Algorithm: \t\t Optimized Basic\n");

    struct timespec start, end;
    clock_gettime(CLOCK_MONOTONIC, &start);

    for(unsigned idx = 0; idx < n; idx++) {
        for(unsigned kdx = 0; kdx < n; kdx++) {
            float reg = g.mtx_A[idx * n + kdx];
            for(unsigned jdx = 0; jdx < n; jdx++) {
                g.mtx_C[idx * n + jdx] += reg * g.mtx_B[kdx * n + jdx];
            }
        }
    }

    clock_gettime(CLOCK_MONOTONIC, &end);
    double sec = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1000000000.0;

    printf("Time to compute: \t %fs\n\n", sec);
    debug_print(g.mtx_C, n);
    fprintf(stderr, "========================================================================\n");

    // clearing C matrix
    memset(g.mtx_C, 0, n * n * sizeof(float));
}
