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
    #define debug_mtx_print(matrix, n)                          \
        do{                                                     \
            printf("----------< DEBUG >----------\n");          \
            printf("%s:\n", #matrix);                           \
            for(int idx = 0; idx < n; idx++) {                  \
                for(int jdx = 0; jdx < n; jdx++) {              \
                    printf("%f, ", matrix[idx * n + jdx]);      \
                }                                               \
                printf("\n");                                   \
            }                                                   \
            printf("\n\n");                                     \
        }while(0)

    #define debug_array_print(array, n)                         \
        do{                                                     \
            printf("microsecs/element:\n");                     \
            for(int idx = 0; idx < n; idx++) {                  \
                printf("%f, ", array[idx]  / ((double)N * N));  \
            }                                                   \
            printf("\n----------< DEBUG >----------\n");        \
            printf("\n\n");                                     \
        }while(0)
#else
    #define debug_mtx_print(matrix, n)
    #define debug_array_print(array, n)
#endif

#define NAIVE_START 1

typedef double mtx_t;



void add(mtx_t *X, size_t stride_X, mtx_t *Y, size_t stride_Y, mtx_t *Z, size_t n) {
    
    for(size_t idx = 0; idx < n; idx++) {
        for(size_t jdx = 0; jdx < n; jdx++) {
            Z[idx * n + jdx] = X[idx * stride_X + jdx] + Y[idx * stride_Y + jdx];
        }
    }
}


void sub(mtx_t *X, size_t stride_X, mtx_t *Y, size_t stride_Y, mtx_t *Z, size_t n) {
    
    for(size_t idx = 0; idx < n; idx++) {
        for(size_t jdx = 0; jdx < n; jdx++) {
            Z[idx * n + jdx] = X[idx * stride_X + jdx] - Y[idx * stride_Y + jdx];
        }
    }
}


static void strassen(mtx_t *A, size_t stride_A, mtx_t *B, size_t stride_B, mtx_t *C, size_t n) {

    // threshold to stop recursion
    if(n <= NAIVE_START) {
        // naive
        for(size_t idx = 0; idx < n; idx++) {
            for(size_t jdx = 0; jdx < n; jdx++) {
                for(size_t kdx = 0; kdx < n; kdx++) {
                    C[idx * n + jdx] += A[idx * stride_A + kdx] * B[kdx * stride_B + jdx];
                }
            }
        }
        return;
    }


    // preparations - matrix block parts pointers
    size_t half = n / 2;
    mtx_t *A11 = A, *A12 = A + half, *A21 = A + half * stride_A, *A22 = A + half * stride_A + half;
    mtx_t *B11 = B, *B12 = B + half, *B21 = B + half * stride_B, *B22 = B + half * stride_B + half;
    
    // 7 results + 2 temp matrices
    mtx_t *I = calloc(half * half, sizeof(mtx_t));
    mtx_t *II = calloc(half * half, sizeof(mtx_t));
    mtx_t *III = calloc(half * half, sizeof(mtx_t));
    mtx_t *IV = calloc(half * half, sizeof(mtx_t));
    mtx_t *V = calloc(half * half, sizeof(mtx_t));
    mtx_t *VI = calloc(half * half, sizeof(mtx_t));
    mtx_t *VII = calloc(half * half, sizeof(mtx_t));
    mtx_t *temp1 = calloc(half * half, sizeof(mtx_t));
    mtx_t *temp2 = calloc(half * half, sizeof(mtx_t));
    if(I == NULL || II == NULL || III == NULL || IV == NULL || V == NULL || VI == NULL ||
       VII == NULL || temp1 == NULL || temp2 == NULL) {
        fprintf(stderr, "Memory allocation failed!");
        exit(1);
    }


    //        temp1         temp2
    // I = (A11 + A22) * (B11 + B22)
    add(A11, stride_A, A22, stride_A, temp1, half);
    add(B11, stride_B, B22, stride_B, temp2, half);
    strassen(temp1, half, temp2, half, I, half);

    //        temp1
    // I = (A21 + A22) * B11
    add(A21, stride_A, A22, stride_A, temp1, half);
    strassen(temp1, half, B11, stride_B, II, half);

    //                temp2
    // III = A11 * (B12 - B22)
    sub(B12, stride_B, B22, stride_B, temp2, half);
    strassen(A11, stride_A, temp2, half, III, half);
 
    //               temp2
    // IV = A22 * (B21 - B11)
    sub(B21, stride_B, B11, stride_B, temp2, half);
    strassen(A22, stride_A, temp2, half, IV, half);
 
    //        temp1
    // V = (A11 + A12) * B22
    add(A11, stride_A, A12, stride_A, temp1, half);
    strassen(temp1, half, B22, stride_B, V, half);
 
    //         temp1         temp2
    // VI = (A21 - A11) * (B11 + B12)
    sub(A21, stride_A, A11, stride_A, temp1, half);
    add(B11, stride_B, B12, stride_B, temp2, half);
    strassen(temp1, half, temp2, half, VI, half);
 
    //          temp1         temp2
    // VII = (A12 - A22) * (B21 + B22)
    sub(A12, stride_A, A22, stride_A, temp1, half);
    add(B21, stride_B, B22, stride_B, temp2, half);
    strassen(temp1, half, temp2, half, VII, half);

    
    mtx_t *C11 = C, *C12 = C + half, *C21 = C + half * n, *C22 = C + half * n + half;

    // results into C
    for(size_t idx = 0; idx < half; idx++) {
        for(size_t jdx = 0; jdx < half; jdx++) {
            size_t pos = idx * half + jdx;

            C11[idx * n + jdx] += I[pos] + IV[pos] - V[pos] + VII[pos];
            C12[idx * n + jdx] += III[pos] + V[pos];
            C21[idx * n + jdx] += II[pos] + IV[pos];
            C22[idx * n + jdx] += I[pos] + III[pos] - II[pos] + VI[pos];
        }
    }

    // freeing
    free(I);
    free(II);
    free(III);
    free(IV);
    free(V);
    free(VI);
    free(VII);
    free(temp1);
    free(temp2);
}


// insertion sort
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
    debug_mtx_print(A, N);
    debug_mtx_print(B, N);

    // how many times will the algorithm run
    int run_count;
    switch(N) {
        case 2: case 4: case 8: case 16: case 32: case 64: case 128:
            run_count = 10000; break;
        case 256: case 512:
            run_count = 1000; break;
        case 1024:
            run_count = 100; break;
        case 2048:
            run_count = 10; break;
        default:
            run_count = 1;
    }

    double total_time[run_count];
    int step = run_count / 10;
    fprintf(stdout, "\n[INFO] \t Matrix %ldx%ld -> Strassen -> Run count: %d\n", (size_t)N, (size_t)N, run_count);


    // algorithm measurement
    fprintf(stderr, "Progress: ");
    for(int idx = 0; idx < run_count; idx++) {
        memset(C, 0, (size_t)N * N * sizeof(mtx_t));

        struct timespec start, end;
        clock_gettime(CLOCK_MONOTONIC, &start);

        strassen(A, N, B, N, C, N);

        clock_gettime(CLOCK_MONOTONIC, &end);
        unsigned long ns = (end.tv_sec - start.tv_sec) * 1000000000.0 + (end.tv_nsec - start.tv_nsec);
        unsigned long us = (ns + 500) / 1000;   // microseconds
        total_time[idx] = (double)us;

        if(step == 0) {
            continue;
        }
        if((idx + 1) % step == 0) {
            fprintf(stderr, "▰");
        }
    }
    fprintf(stderr, "\n\n");
    
    // getting median and final result
    debug_mtx_print(C, N);
    sort(total_time, run_count);
    int median = run_count / 2;
    debug_array_print(total_time, run_count);


    printf("Total time: \t %fs  |  microsecs/element: %f  |  microsecs/mul: %f\n\n",
           total_time[median] / 1000000.0, total_time[median] / ((double)N * N),
           total_time[median] / ((double)N * N * N)); 


    free(A);
    free(B);
    free(C);

    return 0;
}
