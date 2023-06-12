# test preprocessing
find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})
add_executable(test_preprocessing
        test/test_gmsh_mesh.cpp)
target_link_libraries(test_preprocessing ${GTEST_LIBRARIES} gtest_main pthread solver)
