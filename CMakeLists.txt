cmake_minimum_required(VERSION 3.26)
project(OpenMP_exercises C)
set(CMAKE_C_STANDARD 11)
set(CMAKE_EXE_LINKER_FLAGS "-static-libgcc -static-libstdc++")

find_package(OpenMP REQUIRED)
find_package(MPI REQUIRED)

add_executable(OpenMP_exercises
        openMP/pi/pi.c
        main.c
        benchmarking.h
        benchmarking.c
        typedefs.h
        openMP/fixit/fixit.c
        openMP/fixit/fixit.h
        openMP/matrix_vector/mv.h
        openMP/matrix_vector/mv.c
        MPI/integration/integration.c
        MPI/integration/integration.h
        utils/utils.h
        MPI/MatrixVector/matrix_vector_mpi.c
        MPI/MatrixVector/matrix_vector_mpi.h
        MPI/ping_pong/pingPong-a.h
        MPI/ping_pong/pingPong-a.c
        MPI/ping_pong/pingPong-b.c
        MPI/Jacobi/jacobi.c
        MPI/Jacobi/jacobi.h
        openMP/matrix_matrix/matrix-matrix.c
        openMP/matrix_matrix/matrix-matrix.h
        MPI/ping_pong/pingPong-e.c

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
         -mavx2 -mfma
)