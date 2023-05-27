cmake_minimum_required(VERSION 3.25)
project(solver)
set(CMAKE_CXX_STANDARD 17)

# External libs
include_directories(demo_vtu11 /home/gordan/Programs/vtu11)
set(json_dir "/home/gordan/Programs/json")
include_directories(${json_dir}/include)
set(deallii_dir "/home/gordan/Programs/dealii/dealii-install")
include_directories(${deallii_dir}/include/)
include_directories(/home/gordan/Programs/cxxopts/include)

# solver
add_library(solver SHARED
        src/LinearSolver.cpp
        src/ExportVtu.cpp
        src/MagneticFluxPostprocessor.cpp
        src/MatIDPostprocessor.cpp
        src/MagneticEnergyPostprocessor.cpp src/MagneticEnergyDensityPostprocessor.cpp src/NewtonSolver.cpp)
target_include_directories(solver PUBLIC include)
target_link_libraries(solver debug
        /home/gordan/Programs/dealii/dealii-install/lib/libdeal_II.so
        tbb)

# Locate GTest
find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})

# test solver
add_executable(test_solver
        test/test_LinearSolver.cpp
        test/test_MagnetComputations.cpp test/test_NewtonSolver.cpp)

target_link_libraries(test_solver ${GTEST_LIBRARIES} gtest_main pthread solver)

# test export
add_executable(test_export
        test/test_ExportVtu.cpp test/test_MagneticFluxPostprocessor.cpp test/test_MagneticEnergyPostprocessor.cpp test/test_MagneticEnergyDensityPostprocessor.cpp)
target_link_libraries(test_export ${GTEST_LIBRARIES} gtest_main pthread solver)

# test mesh
add_executable(test_mesh
                test/test_gmsh_mesh.cpp)
target_link_libraries(test_mesh ${GTEST_LIBRARIES} gtest_main pthread solver)

# aml
add_executable(aml AppMagnetostaticLinear.cpp)
target_link_libraries(aml optimized
        /home/gordan/Programs/dealii/dealii-install/lib/libdeal_II.so
        tbb
        solver)

# amp
add_executable(amp AppMagnetostaticPicard.cpp)
target_link_libraries(amp optimized
        /home/gordan/Programs/dealii/dealii-install/lib/libdeal_II.so
        tbb
        solver)