#include <stdio.h>
#include <omp.h>
#include "pi/pi.h"
#include "benchmarking.h"
#include "fixit/fixit.h"
#include "matrix_vector/mv.h"
#include "MPI/integration/integration.h"
#include "utils/utils.h"



void strassen_example(){
    float* matrixA, *matrixB;
    int n, m;

}

int main(int argc, char** argv) {
    int n = 150000000;
    int max_threads = omp_get_max_threads();
    printf("runtime generic sequential:\n");
    benchmark(generic_pi,  n);
    benchmark(sequential_pi,  n);
    printf("runtime parallel:\n");
    for(int i = 1; i < max_threads + 1; i ++){
        benchmarki(parallel_pi, i, n);
    }
//    dot();
    //integration_dispatch(argc, argv);
    //call();
    return 0;
}
