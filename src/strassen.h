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
 * @brief
 * 
 * @param n
 * 
 */
void strassen(unsigned n);

/**
 * @brief
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
void compute_strassen(mtx_t *A, unsigned stride_A, mtx_t *B, unsigned stride_B,
                      mtx_t *C, unsigned stride_C, unsigned n);

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
void add(mtx_t *X, unsigned stride_X, mtx_t *Y, unsigned stride_Y,
         mtx_t *Z, unsigned stride_Z, unsigned n);

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
void sub(mtx_t *X, unsigned stride_X, mtx_t *Y, unsigned stride_Y,
         mtx_t *Z, unsigned stride_Z, unsigned n);

#endif // STRASSEN_H
