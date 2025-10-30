cmake_minimum_required(VERSION 3.26)
project(OpenMP_exercises C)
set(CMAKE_C_STANDARD 11)

find_package(OpenMP REQUIRED)

add_executable(OpenMP_exercises
        pi/pi.c
        main.c
        benchmarking.h
        benchmarking.c
        typedefs.h

)



target_link_libraries(OpenMP_exercises PRIVATE
        OpenMP::OpenMP_C
)

target_compile_options(OpenMP_exercises PRIVATE
         -mavx2  -mfma -march=native
)
