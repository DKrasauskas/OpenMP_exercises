//
// Created by dominykas on 11/3/25.
//

#ifndef OPENMP_EXERCISES_PINGPONG_A_H
#define OPENMP_EXERCISES_PINGPONG_A_H
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#define MAX_ARRAY_SIZE (1<<10)
void ping_pong_a(int argc, char **argv);
void ping_pong_b(int argc, char **argv);
void ping_pong_c(int argc, char **argv);
void ping_pong_d(int argc, char **argv);
void generate_ab(int argc, char **argv);
#define MAX_ARRAY_SIZE (1<<10)



#endif //OPENMP_EXERCISES_PINGPONG_A_H
