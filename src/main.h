/**
 * BT - Fast algorithms for matrix multiplication
 * 
 * @file main.h
 * @details header file for main - defined functions and libraries
 * 
 * @author xpetkob00 Branislav Peťko
*/

#ifndef MAIN_H
#define MAIN_H

#define ALGORITHM_COUNT 2

#include <stdio.h>
#include <stdlib.h>
#include <signal.h>
#include <time.h>
#include <string.h>
#include <ctype.h>
#include <stdbool.h>

#include "error.h"
#include "basic.h"
#include "strassen.h"

typedef struct global_struct {
    float *mtx_A;
    float *mtx_B;
    float *mtx_C;
} global;

extern global g;


/**
 * @brief parses program arguments and controls the whole program
 * 
 * @param argc number of arguments (argument count)
 * @param argv list of arguments
 * 
 */
void arg_handler(int argc, char *argv[]);

/**
 * @brief prints help message
 */
void print_help();

/**
 * @brief handles SIGINT and SIGTERM signals - cleanup of all resources
 * 
 * @param sig signal number receiver (required by signal function)
 * 
 */
void handle_signal(int sig);

/**
 * @brief main starting point for program
 * 
 * @param argc number of arguments (argument count)
 * @param argv list of arguments
 * 
 * @return 0 if successful
 */
int main(int argc, char *argv[]);


#endif // MAIN_H
