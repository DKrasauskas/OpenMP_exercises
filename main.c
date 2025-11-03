
#include <stdio.h>
#include <omp.h>
#include "pi/pi.h"
#include "benchmarking.h"
#include "fixit/fixit.h"
#include "matrix_vector/mv.h"
#include "MPI/integration/integration.h"
#include "utils/utils.h"
#include "MPI/ping_pong/pingPong-a.h"


void strassen_example(){
    float* matrixA, *matrixB;
    int n, m;

}

int main(int argc, char** argv) {
    int n = 1500000000;
    int max_threads = omp_get_max_threads();
    printf("runtime generic sequential:\n");
    benchmark(generic_pi,  n);
    #ifdef MVX_ARCH
    benchmark(sequential_pi,  n);
    printf("runtime parallel:\n");
    #endif
    for(int i = 1; i < max_threads + 1; i ++){
        benchmarki(reduction_parallel_pi, i, n);
    }

    return 0;
}
