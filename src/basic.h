/** 
 * BT - Fast algorithms for matrix multiplication
 * 
 * @file basic.h
 * @details basic algorithm for matrix multiplication
 * 
 * @author xpetkob00 Branislav Peťko
*/

#ifndef BASIC_H
#define BASIC_H

#define _POSIX_C_SOURCE 199309L

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "main.h"


/**
 * @brief function to run basic algorithm for matrix multiplication
 * 
 * @param n matrix size (n x n)
 * 
 */
void basic(unsigned n);

#endif // BASIC_H
