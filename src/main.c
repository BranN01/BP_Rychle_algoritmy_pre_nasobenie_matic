/**
 * BT - Fast algorithms for matrix multiplication
 * 
 * @file main.c
 * @details main starting point for program
 * 
 * @author xpetkob00 Branislav Peťko
*/

#include "main.h"

// global variables - needed for cleanup
global g = {
    .mtx_A = NULL,
    .mtx_B = NULL,
    .mtx_C = NULL
};


void arg_handler(int argc, char *argv[]) {

    
}

void print_help() {

    printf("==============================    HELP    ==============================\n");
    printf("Desc\n");
    printf("\n\n");
    printf("\nUSAGE:\n");
    printf("\n\n\n");
    printf("OPTIONS:\n");
    printf("-h / --help \t\t Writes usage instructions\n");
    printf("...\n");
    printf("EXAMPLES:\n");
    printf("\n");
    printf("\nNOTES:\n");
    printf("- ...\n");
}


void handle_signal(int sig) {
    (void)sig; // unused parameter which is required by signal function

    if(g.mtx_A != NULL) {
        free(g.mtx_A);
    }
    if(g.mtx_B != NULL) {
        free(g.mtx_B);
    }
    if(g.mtx_C != NULL) {
        free(g.mtx_C);
    }
    
    warning("Program has been terminated by SIGTERM or SIGINT! Everything has been cleaned up!");
    error_exit(SIG_ERROR, "Exiting...");
}


int main(int argc, char *argv[]) {

    signal(SIGTERM, handle_signal);
    signal(SIGINT, handle_signal);
    srand(time(NULL));

    arg_handler(argc, argv);

    return SUCCESS;
}

