/** 
 * BT - Fast algorithms for matrix multiplication
 * 
 * @file error.h
 * @details functions for handling errors
 * 
 * @author xpetkob00 Branislav Peťko
*/

#ifndef ERROR_H
#define ERROR_H

#define SUCCESS 0
#define ARG_ERROR 1
#define SIG_ERROR 90
#define UNEXPECTED_ERROR 99

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>


/**
 * @brief function to print a warning message to stderr
 * 
 * @param fmt format, in which to print warning message
 * 
 * @note Usage: warning("Something bad happened: %s", "something");
 */
void warning(const char* fmt, ...);


/**
 * @brief function to print an error message to stderr and exit with specific exit code
 * 
 * @param error_code error number to exit program with
 * @param fmt format, in which to print error message
 * 
 * @note Usage: error_exit(2, "Something bad happened: %s", "something");
 */
void error_exit(int error_code, const char* fmt, ...);

#endif // ERROR_H
