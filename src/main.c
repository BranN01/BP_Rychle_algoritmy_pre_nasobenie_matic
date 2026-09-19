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
    
    int temp = 0;
    bool algorithms[ALGORITHM_COUNT] = {false};

    // checking arguments validity
    for(int idx = 1; idx < argc; idx++) {

        // convert all characters to lower
        char *str = argv[idx];
        for(int jdx = 0; str[jdx] != '\0'; jdx++) {
            str[jdx] = tolower(str[jdx]);
        }

        if(strcmp(str, "-h") == 0 || strcmp(str, "--help") == 0) {
            print_help();
            return;
        }
        else if(strcmp(str, "-a") == 0) {
            for(int jdx = 0; jdx < ALGORITHM_COUNT; jdx++) {
                algorithms[jdx] = true;
            }
        }
        else if(strcmp(str, "-b") == 0) {
            algorithms[0] = true;
        }
        else if(strcmp(str, "-s") == 0) {
            algorithms[1] = true;
        }
        else {
            temp = atoi(str);
            if(temp == 0) {
                error_exit(ARG_ERROR, "Invalid argument!");
            }
        }
    }

    // check for size validity
    if(temp < 2 || ((temp & (temp - 1)) != 0)) {
        error_exit(ARG_ERROR, "Invalid matrix size! (Must be power of 2 and minimum 2)");
    }
    unsigned size = (unsigned)temp;

    // ensure at least one algorithm is selected
    for(int idx = 0; idx <= ALGORITHM_COUNT; idx++) {

        if(idx == ALGORITHM_COUNT) {
            error_exit(ARG_ERROR, "At least one algorithm must be selected!");
        }
        if(algorithms[idx]) {
            break;
        }
    }


    // allocating space for matrices
    g.mtx_A = calloc(size * size, sizeof(float));
    g.mtx_B = calloc(size * size, sizeof(float));
    g.mtx_C = calloc(size * size, sizeof(float));
    if(g.mtx_C == NULL || g.mtx_C == NULL || g.mtx_C == NULL) {
        error_exit(UNEXPECTED_ERROR, "Memory allocation failed!");
    }

    // filling matrices A and B with random numbers
    fprintf(stderr, "[INFO] Filling matrices with random data...\n");
    for(unsigned idx = 0; idx < size * size; idx++) {
        g.mtx_A[idx] = (rand() % 100001) / 1000.0f;
        g.mtx_B[idx] = (rand() % 100001) / 1000.0f;
    }
    debug_print(g.mtx_A, size);
    debug_print(g.mtx_B, size);

    
    // running selected algorithms
    fprintf(stderr, "[INFO] Running selected algorithms...\n\n\n");
    fprintf(stderr, "========================================================================\n");
    for(int idx = 0; idx < ALGORITHM_COUNT; idx++) {
        if(algorithms[idx]) {
            switch (idx) {
                case 0: basic(size); basic_optimized(size); break;
                case 1: strassen(size); break;
                // TODO
            }
        }
    }
    

    free(g.mtx_A);
    g.mtx_A = NULL;
    free(g.mtx_B);
    g.mtx_B = NULL;
    free(g.mtx_C);
    g.mtx_C = NULL;
}

void print_help() {

    printf("==============================    HELP    ==============================\n");
    printf("Executes selected matrix multiplication algorithms on randomly generated\n");
    printf("matrices of the specified size and prints the execution time for each.\n\n");
    printf("\nUSAGE:\n");
    printf("./mtx-multiplication [ALGORITHMS] SIZE\n\n\n");
    printf("OPTIONS:\n");
    printf("-h / --help \t\t Writes usage instructions\n");
    printf("-a \t\t\t Executes all algorithms\n");
    printf("-b \t\t\t Executes the basic algorithm\n"); // [TODO] add later
    printf("-s \t\t\t Executes the Strassen algorithm\n");
    printf("SIZE \t\t\t Size of the square matrix\n\n\n");
    printf("EXAMPLES:\n");
    printf("./mtx-multiplication --help\n");
    printf("./mtx-multiplication -a 64\n");
    printf("./mtx-multiplication -b ... 128\n\n"); // [TODO] later
    printf("\nNOTES:\n");
    printf("- At least one algorithm must be specified!\n");
    printf("- SIZE must be power of 2 and minimum 2\n");
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
