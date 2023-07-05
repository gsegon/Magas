# solver
add_library(solver OBJECT
        LinearSolver.cpp
        PicardSolver.cpp
        NewtonSolver.cpp)

find_package(deal.II 9.4.0 REQUIRED)

DEAL_II_INITIALIZE_CACHED_VARIABLES()
DEAL_II_SETUP_TARGET(solver)

target_compile_options(solver PRIVATE -Wall -Wextra)
target_include_directories(solver PUBLIC include)
target_link_libraries(solver
        optimized
        preprocessing
        tbb)

# test solver
# Locate GTest
find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})

add_executable(test_solver
        test/test_LinearSolver.cpp
        test/test_PicardSolver.cpp
        test/test_NewtonSolver.cpp)

target_link_libraries(test_solver ${GTEST_LIBRARIES} gtest_main pthread solver preprocessing)
DEAL_II_SETUP_TARGET(test_solver)

add_test(test_solver_gtest test_solver)
