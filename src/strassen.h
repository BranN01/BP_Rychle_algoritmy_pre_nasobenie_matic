/** 
 * BT - Fast algorithms for matrix multiplication
 * 
 * @file strassen.h
 * @details Strassen algorithm for matrix multiplication
 * 
 * @author xpetkob00 Branislav Peťko
*/

#ifndef STRASSEN_H
#define STRASSEN_H

#define _POSIX_C_SOURCE 199309L
#define NAIVE_START 16

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "main.h"


/**
 * @brief function to run Strassen algorithm for matrix multiplication
 * 
 * @param n matrix size (n x n)
 * 
 */
void strassen(size_t n);


/**
 * @brief main Strassen algorithm implementation
 * 
 * @param A
 * @param stride_A
 * @param B
 * @param stride_B
 * @param C
 * @param stride_C
 * @param n
 * 
 */
void compute_strassen(mtx_t *A, unsigned stride_A, mtx_t *B, size_t stride_B,
                      mtx_t *C, unsigned stride_C, size_t n);


/**
 * @brief
 * 
 * @param X
 * @param stride_X
 * @param Y
 * @param stride_Y
 * @param Z
 * @param stride_Z
 * @param n
 * 
 */                     
void add(mtx_t *X, size_t stride_X, mtx_t *Y, size_t stride_Y,
         mtx_t *Z, size_t stride_Z, size_t n);


/**
 * @brief
 * 
 * @param X
 * @param stride_X
 * @param Y
 * @param stride_Y
 * @param Z
 * @param stride_Z
 * @param n
 * 
 */                     
void sub(mtx_t *X, size_t stride_X, mtx_t *Y, size_t stride_Y,
         mtx_t *Z, size_t stride_Z, size_t n);

#endif // STRASSEN_H
