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

#endif // STRASSEN_H
