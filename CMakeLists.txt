cmake_minimum_required(VERSION 3.16)
project(Magas)
set(CMAKE_CXX_STANDARD 17)

enable_testing()

# Modules
add_subdirectory(src/material)
add_subdirectory(src/sources)
add_subdirectory(src/preprocessing)
add_subdirectory(src/postprocessing)
add_subdirectory(src/solver)

# Apps
add_subdirectory(src)

# Examples
add_subdirectory(examples)
