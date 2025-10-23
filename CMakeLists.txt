cmake_minimum_required(VERSION 4.0)
project(algdat C)

set(CMAKE_C_STANDARD 17)

add_executable(algdat main.c
        )
add_executable(algdat2 main2.c)
add_executable(algdat3 main3.c)
add_executable(algdat4 main4.c)
add_executable(algdat5 main5.c)
add_executable(algdat6 main6.c)
add_executable(merge merge.c)
add_executable(moore  boyer-moore.c)
add_executable(test test.c)
