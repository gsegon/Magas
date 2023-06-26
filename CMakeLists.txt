cmake_minimum_required(VERSION 3.25)
project(solver)
set(CMAKE_CXX_STANDARD 17)

enable_testing()

add_subdirectory(src/preprocessing)
add_subdirectory(src/postprocessing)
add_subdirectory(src/solver)
add_subdirectory(src)

add_subdirectory(examples)
