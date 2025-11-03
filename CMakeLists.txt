cmake_minimum_required(VERSION 3.26)
project(OpenMP_exercises C)
set(CMAKE_C_STANDARD 11)

find_package(OpenMP REQUIRED)
find_package(MPI REQUIRED)

add_executable(OpenMP_exercises
        pi/pi.c
        main.c
        benchmarking.h
        benchmarking.c
        typedefs.h
        fixit/fixit.c
        fixit/fixit.h
        matrix_vector/mv.h
        matrix_vector/mv.c
        MPI/integration/integration.c
        MPI/integration/integration.h
        utils/utils.h
        MPI/MatrixVector/matrix_vector_mpi.c
        MPI/MatrixVector/matrix_vector_mpi.h
        MPI/ping_pong/pingPong-a.h
        MPI/ping_pong/pingPong-a.c
        MPI/ping_pong/pingPong-b.c

)



target_link_libraries(OpenMP_exercises PRIVATE
        OpenMP::OpenMP_C
        MPI::MPI_C
        m
)

#target_compile_options(OpenMP_exercises PRIVATE
#        -mavx2  -mfma -march=native
#)

target_compile_options(OpenMP_exercises PRIVATE
         -O2 -mavx2 -mfma
)