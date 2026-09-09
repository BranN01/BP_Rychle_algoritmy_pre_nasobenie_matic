/** 
 * BT - Fast algorithms for matrix multiplication
 * 
 * @file error.c
 * @details implementation of error handler
 * 
 * @author xpetkob00 Branislav Peťko
*/

#include "error.h"


void warning(const char *fmt, ...) {

    va_list args;
    va_start(args, fmt);

    fprintf(stderr, "Warning: ");
    vfprintf(stderr, fmt, args);
    fprintf(stderr, "\n");

    va_end(args);
}

void error_exit(int error_code, const char *fmt, ...) {
    va_list args;
    va_start(args, fmt);

    fprintf(stderr, "Error: ");
    vfprintf(stderr, fmt, args);
    fprintf(stderr, "\n");

    va_end(args);
    exit(error_code);
}
