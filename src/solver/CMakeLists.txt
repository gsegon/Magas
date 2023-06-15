# solver
add_library(solver OBJECT
        LinearSolver.cpp
        PicardSolver.cpp)

target_compile_options(solver PRIVATE -Wall -Wextra)
target_include_directories(solver PUBLIC include)
set(deallii_dir "/home/gordan/Programs/dealii/dealii-install")
target_include_directories(solver PUBLIC ${deallii_dir}/include/)
target_link_libraries(solver debug
        preprocessing
        ${deallii_dir}/lib/libdeal_II.so
        tbb)

# test solver
# Locate GTest
find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})

add_executable(test_solver
        test/test_LinearSolver.cpp
        test/test_PicardSolver.cpp)

target_link_libraries(test_solver ${GTEST_LIBRARIES} gtest_main pthread solver)