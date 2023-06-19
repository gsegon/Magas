add_library(preprocessing OBJECT
        Dummy.cpp)

target_include_directories(preprocessing PUBLIC include)

set(deallii_dir "/home/gordan/Programs/dealii/dealii-install")
target_include_directories(preprocessing PUBLIC ${deallii_dir}/include)


set(eigen_dir "/home/gordan/Programs/eigen-install")
target_include_directories(preprocessing PUBLIC ${eigen_dir}/include/eigen3)

target_link_libraries(preprocessing debug
        ${deallii_dir}/lib/libdeal_II.so
        tbb)

# test preprocessing
find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})
add_executable(test_preprocessing
        test/test_gmsh_mesh.cpp
        test/test_PeriodicityMapper.cpp
        test/test_OverlapPointsTransformation.cpp)
target_include_directories(test_preprocessing PUBLIC include)
target_include_directories(test_preprocessing PUBLIC .)
target_link_libraries(test_preprocessing ${GTEST_LIBRARIES} gtest_main pthread preprocessing)

add_test(test_preprocessing_gtest test_preprocessing)
