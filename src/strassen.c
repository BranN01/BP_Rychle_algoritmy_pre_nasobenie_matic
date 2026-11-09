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
#include <string.h>

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


static void strassen(mtx_t *A, size_t stride_A, mtx_t *B, size_t stride_B,
                     mtx_t *C, size_t stride_C, size_t n, mtx_t *buffer) {

    // threshold to stop recursion
    if(n <= NAIVE_START) {
        for(size_t idx = 0; idx < n; idx++) {
            for(size_t jdx = 0; jdx < n; jdx++) {
                for(size_t kdx = 0; kdx < n; kdx++) {
                    C[idx * stride_C + jdx] += A[idx * stride_A + kdx] * B[kdx * stride_B + jdx];
                }
            }
        }
        return;
    }


    size_t half = n / 2;
    mtx_t *A11 = A, *A12 = A + half, *A21 = A + half * stride_A, *A22 = A + half * stride_A + half;
    mtx_t *B11 = B, *B12 = B + half, *B21 = B + half * stride_B, *B22 = B + half * stride_B + half;
    mtx_t *I = buffer, *II = I + half * half, *III = II + half * half, *IV = III + half * half,
          *V = IV + half * half, *VI = V + half * half, *VII = VI + half * half;
    mtx_t *temp1 = VII + half * half, *temp2 = temp1 + half * half, *next_free = temp2 + half * half;

    memset(I, 0, 7 * half * half * sizeof(mtx_t));


    //        temp1         temp2
    // I = (A11 + A22) * (B11 + B22)
    add(A11, stride_A, A22, stride_A, temp1, half, half);
    add(B11, stride_B, B22, stride_B, temp2, half, half);
    strassen(temp1, half, temp2, half, I, half, half, next_free);

    //        temp1
    // I = (A21 + A22) * B11
    add(A21, stride_A, A22, stride_A, temp1, half, half);
    strassen(temp1, half, B11, stride_B, II, half, half, next_free);

    //                temp2
    // III = A11 * (B12 - B22)
    sub(B12, stride_B, B22, stride_B, temp2, half, half);
    strassen(A11, stride_A, temp2, half, III, half, half, next_free);
 
    //               temp2
    // IV = A22 * (B21 - B11)
    sub(B21, stride_B, B11, stride_B, temp2, half, half);
    strassen(A22, stride_A, temp2, half, IV, half, half, next_free);
 
    //        temp1
    // V = (A11 + A12) * B22
    add(A11, stride_A, A12, stride_A, temp1, half, half);
    strassen(temp1, half, B22, stride_B, V, half, half, next_free);
 
    //         temp1         temp2
    // VI = (A21 - A11) * (B11 + B12)
    sub(A21, stride_A, A11, stride_A, temp1, half, half);
    add(B11, stride_B, B12, stride_B, temp2, half, half);
    strassen(temp1, half, temp2, half, VI, half, half, next_free);
 
    //          temp1         temp2
    // VII = (A12 - A22) * (B21 + B22)
    sub(A12, stride_A, A22, stride_A, temp1, half, half);
    add(B21, stride_B, B22, stride_B, temp2, half, half);
    strassen(temp1, half, temp2, half, VII, half, half, next_free);

    
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


void sort(double *array, int n) {

    for(int idx = 1; idx < n; idx++) {
        double key = array[idx];
        int jdx = idx - 1;

        while(jdx >= 0 && array[jdx] > key) {
            array[jdx + 1] = array[jdx];
            jdx--;
        }
        array[jdx + 1] = key;
    }
}


int main() {
    
    // allocating space for matrices
    mtx_t *A = calloc((size_t)N * N, sizeof(mtx_t));
    mtx_t *B = calloc((size_t)N * N, sizeof(mtx_t));
    mtx_t *C = calloc((size_t)N * N, sizeof(mtx_t));
    mtx_t *buffer = calloc((size_t)N * N * 3, sizeof(mtx_t)); // buffer for Strassen matrices
    if(A == NULL || B == NULL || C == NULL || buffer == NULL) {
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

    int run_count;
    if(N >= 2 && N <= 2048) {
        run_count = 1000;
    }
    else {
        run_count = 100;
    }
    double total_time[run_count], us_element[run_count], us_mul[run_count];
    int step = run_count / 10;
    fprintf(stdout, "\n[INFO] \t Matrix %ldx%ld -> Strassen -> Run count: %d\n", (size_t)N, (size_t)N, run_count);


    fprintf(stderr, "Progress: ");
    for(int idx = 0; idx < run_count; idx++) {
        memset(C, 0, (size_t)N * N * sizeof(mtx_t));

        struct timespec start, end;
        clock_gettime(CLOCK_MONOTONIC, &start);

        strassen(A, N, B, N, C, N, N, buffer);

        clock_gettime(CLOCK_MONOTONIC, &end);
        unsigned long ns = (end.tv_sec - start.tv_sec) * 1000000000.0 + (end.tv_nsec - start.tv_nsec);
        double sec = ns / 1000000000.0;
        unsigned long us = (ns + 500) / 1000;   // microseconds
        total_time[idx] = sec;
        us_element[idx] = us / ((double)N * N);
        us_mul[idx] = us / ((double)N * N * N);

        if((idx + 1) % step == 0) {
            fprintf(stderr, "▰");
        }
    }
    fprintf(stderr, "\n\n");
    
    debug_print(C, N);
    sort(total_time, run_count);
    sort(us_element, run_count);
    sort(us_mul, run_count);
    int median = run_count / 2;
    printf("Total time: \t %fs  |  microsecs/element: %f  |  microsecs/mul: %f\n\n",
           total_time[median], us_element[median], us_mul[median]); 


    free(A);
    free(B);
    free(C);
    free(buffer);

    return 0;
}