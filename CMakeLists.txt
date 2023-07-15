cmake_minimum_required(VERSION 3.16)
project(Magas)
set(CMAKE_CXX_STANDARD 17)

enable_testing()

set(DEAL_II_DIR ~/Programs/dealii/dealii-9.5.0/build)
#set(DEAL_II_DIR ~/Programs/dealii/dealii-9.5.0/build-static)
find_package(deal.II 9.4 REQUIRED HINTS ${DEAL_II_DIR} ../ ../../ $ENV{DEAL_II_DIR})

# Modules
add_subdirectory(src/material)
add_subdirectory(src/sources)
add_subdirectory(src/input)
add_subdirectory(src/preprocessing)
add_subdirectory(src/postprocessing)
add_subdirectory(src/solver)

# Apps
add_subdirectory(src)

# Examples
add_subdirectory(examples)
