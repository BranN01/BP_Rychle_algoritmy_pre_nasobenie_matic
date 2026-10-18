/** 
 * BT - Fast algorithms for matrix multiplication
 * 
 * @file strassen.h
 * @details implementation of Strassen algorithm for matrix multiplication
 * 
 * @author xpetkob00 Branislav Peťko
*/

#include "strassen.h"


void strassen(size_t n) {

    printf("Algorithm: \t Strassen\n");

    struct timespec start, end;
    clock_gettime(CLOCK_MONOTONIC, &start);

    // main algorithm to compute the multiplication
    //compute_strassen(g.mtx_A, g.mtx_B, n);

    clock_gettime(CLOCK_MONOTONIC, &end);
    double sec = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1000000000.0;

    printf("Total time: \t %fs  |  microsecs/element: %f  |  microsecs/mul: %f\n\n",
           sec, sec * 1000000.0 / ((double)n * n), sec * 1000000.0 / ((double)n * n * n));
    debug_print(g.mtx_C, n);
    fprintf(stderr, "========================================================================\n");

    // clearing C matrix
    memset(g.mtx_C, 0, n * n * sizeof(mtx_t));
}

/*void compute_strassen(mtx_t *A, size_t stride_A, mtx_t *B, size_t stride_B,
                        mtx_t *C, size_t stride_C, size_t n) {

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
}*/
