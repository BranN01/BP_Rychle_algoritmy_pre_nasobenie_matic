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

    //compute_strassen(g.mtx_A, g.mtx_B, n);

    clock_gettime(CLOCK_MONOTONIC, &end);
    double sec = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1000000000.0;

    printf("Time to compute: \t %fs\n\n", sec);
    debug_print(g.mtx_C, n);
    fprintf(stderr, "========================================================================\n");

    // clearing C matrix
    memset(g.mtx_C, 0, n * n * sizeof(mtx_t));
}

/*void compute_strassen(mtx_t *A, unsigned stride_A, mtx_t *B, unsigned stride_B,
                      mtx_t *C, unsigned stride_C, unsigned n) {

    // 
    if(n <= BNAIVE_START) {
        for(unsigned idx = 0; idx < n; idx++) {
            for(unsigned kdx = 0; kdx < n; kdx++) {
                mtx_t reg = A[idx * stride_A + kdx];
                for(unsigned jdx = 0; jdx < n; jdx++) {
                    C[idx * stride_C + jdx] += reg * B[kdx * stride_B + jdx];
                }
            }
        }
        return;
    }


    unsigned half = n / 2;

    mtx_t *A11 = A, *A12 = A + half, *A21 = A + half * stride_A, *A22 = A + half * stride_A + half;
    mtx_t *B11 = B, *B12 = B + half, *B21 = B + half * stride_B, *B22 = B + half * stride_B + half;
    

    // multiplication TODO


    
    mtx_t *C11 = C, *C12 = C + half, *C21 = C + half * stride_C, *C22 = C + half * stride_C + half;

    for(unsigned idx = 0; idx < half; idx++) {
        for(unsigned jdx = 0; jdx < half; jdx++) {
            unsigned pos = idx * half + jdx;

            C11[idx * stride_C + jdx] += I[pos] + IV[pos] - V[pos] + VII[pos];
            C12[idx * stride_C + jdx] += III[pos] + V[pos];
            C21[idx * stride_C + jdx] += II[pos] + IV[pos];
            C22[idx * stride_C + jdx] += I[pos] + III[pos] - II[pos] + VI[pos];
        }
    }
}

void add(mtx_t *X, unsigned stride_X, mtx_t *Y, unsigned stride_Y,
         mtx_t *Z, unsigned stride_Z, unsigned n) {
    
    for(unsigned idx = 0; idx < n; idx++) {
        for(unsigned jdx = 0; jdx < n; jdx++) {
            Z[idx * stride_Z + jdx] = X[idx * stride_X + jdx] + Y[idx * stride_Y + jdx];
        }
    }
}

void sub(mtx_t *X, unsigned stride_X, mtx_t *Y, unsigned stride_Y,
         mtx_t *Z, unsigned stride_Z, unsigned n) {
    
    for(unsigned idx = 0; idx < n; idx++) {
        for(unsigned jdx = 0; jdx < n; jdx++) {
            Z[idx * stride_Z + jdx] = X[idx * stride_X + jdx] - Y[idx * stride_Y + jdx];
        }
    }
}*/
