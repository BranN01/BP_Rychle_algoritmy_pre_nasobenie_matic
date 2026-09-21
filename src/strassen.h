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
#define BASIC_START 16

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
void compute_strassen(float *A, unsigned stride_A, float *B, unsigned stride_B,
                      float *C, unsigned stride_C, unsigned n);

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
void add(float *X, unsigned stride_X, float *Y, unsigned stride_Y,
         float *Z, unsigned stride_Z, unsigned n);

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
void sub(float *X, unsigned stride_X, float *Y, unsigned stride_Y,
         float *Z, unsigned stride_Z, unsigned n);

#endif // STRASSEN_H
