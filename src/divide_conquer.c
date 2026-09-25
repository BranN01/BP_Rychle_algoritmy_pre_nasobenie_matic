/** 
 * BT - Fast algorithms for matrix multiplication
 * 
 * @file divide_conquer.h
 * @details implementation of divide and conquer algorithm for matrix multiplication
 * 
 * @author xpetkob00 Branislav Peťko
*/

#include "divide_conquer.h"


void divide_conquer(unsigned n) {

    printf("Algorithm: \t\t Divide and conquer\n");

    struct timespec start, end;
    clock_gettime(CLOCK_MONOTONIC, &start);

    // main algorithm to compute the multiplication
    compute_divide_conquer(g.mtx_A, g.mtx_B, g.mtx_C, n, n);

    clock_gettime(CLOCK_MONOTONIC, &end);
    double sec = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1000000000.0;

    printf("Time to compute: \t %fs\n\n", sec);
    debug_print(g.mtx_C, n);
    fprintf(stderr, "========================================================================\n");

    // clearing C matrix
    memset(g.mtx_C, 0, n * n * sizeof(mtx_t));
}

void compute_divide_conquer(mtx_t *A, mtx_t *B, mtx_t *C, unsigned stride, unsigned n) {

    // threshold to stop recursion
    if(n == 1) {
        C[0] += A[0] * B[0];
        return;
    }

    // dividing matrix by 2
    unsigned half = n / 2;

    // setting pointers to the matrix parts after divide
    mtx_t *A11 = A, *A12 = A + half, *A21 = A + half * stride, *A22 = A + half * stride + half;
    mtx_t *B11 = B, *B12 = B + half, *B21 = B + half * stride, *B22 = B + half * stride + half;
    mtx_t *C11 = C, *C12 = C + half, *C21 = C + half * stride, *C22 = C + half * stride + half;

    // C11 = A11 * B11 + A12 * B21
    compute_divide_conquer(A11, B11, C11, stride, half);
    compute_divide_conquer(A12, B21, C11, stride, half);

    // C12 = A11 * B12 + A12 * B22
    compute_divide_conquer(A11, B12, C12, stride, half);
    compute_divide_conquer(A12, B22, C12, stride, half);

    // C21 = A21 * B11 + A22 * B21
    compute_divide_conquer(A21, B11, C21, stride, half);
    compute_divide_conquer(A22, B21, C21, stride, half);

    // C22 = A21 * B12 + A22 * B22
    compute_divide_conquer(A21, B12, C22, stride, half);
    compute_divide_conquer(A22, B22, C22, stride, half);
}
