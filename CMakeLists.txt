cmake_minimum_required(VERSION 3.25)
project(solver)
set(CMAKE_CXX_STANDARD 17)



set(deallii_dir "/home/gordan/Programs/dealii/dealii-install")
include_directories(${deallii_dir}/include/)

add_executable(solver linear.cpp linear.h)

target_link_libraries(solver debug
        /home/gordan/Programs/dealii/dealii-install/lib/libdeal_II.so
        tbb
        )
