cmake_minimum_required(VERSION 3.25)
project(solver)
set(CMAKE_CXX_STANDARD 17)


add_library(solver SHARED
        src/LinearSolver.cpp
        src/PicardSolver.cpp)
target_include_directories(solver PUBLIC include)

set(deallii_dir "/home/gordan/Programs/dealii/dealii-install")
include_directories(${deallii_dir}/include/)

target_link_libraries(solver debug
        /home/gordan/Programs/dealii/dealii-install/lib/libdeal_II.so
        tbb)

# Locate GTest
find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})

add_executable(test_solver
        test/test_LinearSolver.cpp
        test/test_PicardSolver.cpp
        test/test_gmsh_mesh.cpp test/test_PicardSolver.cpp)

target_link_libraries(test_solver ${GTEST_LIBRARIES} gtest_main pthread solver)

add_executable(aml AppMagnetostaticLinear.cpp src/PicardSolver.cpp)

set(json_dir "/home/gordan/Programs/json")
include_directories(${json_dir}/include)

target_link_libraries(aml optimized
        /home/gordan/Programs/dealii/dealii-install/lib/libdeal_II.so
        tbb
        solver)