/** 
 * BT - Fast algorithms for matrix multiplication
 * 
 * @file divide_conquer.h
 * @details divide and conquer algorithm for matrix multiplication
 * 
 * @author xpetkob00 Branislav Peťko
*/

#ifndef DIVIDE_CONQUER_H
#define DIVIDE_CONQUER_H

#define _POSIX_C_SOURCE 199309L

#include "main.h"


/**
 * @brief function to run divide and conquer algorithm for matrix multiplication
 * 
 * @param n matrix size (n x n)
 * 
 */
void divide_conquer(unsigned n);


/**
 * @brief main divide and conquer algorithm implementation
 * 
 * @param A pointer to the first matrix
 * @param B pointer to the second matrix
 * @param C pointer to the result matrix
 * @param stride row lenght of the original matrices, used to index submatrices 
 * @param n current submatrix size (n x n)
 * 
 */
void compute_divide_conquer(mtx_t *A, mtx_t *B, mtx_t *C, unsigned stride, unsigned n);

#endif // DIVIDE_CONQUER_H
