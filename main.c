#include <stdio.h>
#include <omp.h>
#include "pi/pi.h"
#include "benchmarking.h"



int main() {
    benchmarki(parallel_pi, 20, 1000000000);
    benchmark(generic_pi,  1000000000);
    benchmark(sequential_pi,  1000000000);
    return 0;
}
