/** 
 * BT - Fast algorithms for matrix multiplication
 * 
 * @file strassen.h
 * @details implementation of Strassen algorithm for matrix multiplication
 * 
 * @author xpetkob00 Branislav Peťko
*/

#include "strassen.h"


void strassen(unsigned n) {

    printf("Algorithm: \t\t Strassen\n");

    struct timespec start, end;
    clock_gettime(CLOCK_MONOTONIC, &start);

    // ALG

    clock_gettime(CLOCK_MONOTONIC, &end);
    double sec = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1000000000.0;

    printf("Time to compute: \t %fs\n\n", sec);
    debug_print(g.mtx_C, n);
    fprintf(stderr, "========================================================================\n");

    // clearing C matrix
    memset(g.mtx_C, 0, n * n * sizeof(float));
}
