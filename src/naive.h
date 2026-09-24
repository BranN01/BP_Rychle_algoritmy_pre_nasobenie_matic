/** 
 * BT - Fast algorithms for matrix multiplication
 * 
 * @file naive.h
 * @details naive algorithm for matrix multiplication
 * 
 * @author xpetkob00 Branislav Peťko
*/

#ifndef NAIVE_H
#define NAIVE_H

#define _POSIX_C_SOURCE 199309L

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "main.h"


/**
 * @brief function to run naive algorithm for matrix multiplication
 * 
 * @param n matrix size (n x n)
 * 
 */
void naive(unsigned n);

/**
 * @brief function to run naive optimized algorithm for matrix multiplication
 * 
 * @param n matrix size (n x n)
 * 
 */
void naive_optimized(unsigned n);

#endif // NAIVE_H
