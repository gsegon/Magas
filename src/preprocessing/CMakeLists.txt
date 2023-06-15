add_library(preprocessing OBJECT
        PeriodicityMapper.cpp)

target_include_directories(preprocessing PUBLIC include)
set(deallii_dir "/home/gordan/Programs/dealii/dealii-install")
target_include_directories(preprocessing PUBLIC ${deallii_dir}/include/)
target_link_libraries(preprocessing debug
        ${deallii_dir}/lib/libdeal_II.so
        tbb)

# test preprocessing
find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})
add_executable(test_preprocessing
        test/test_gmsh_mesh.cpp
        test/test_PeriodicityMapper.cpp)
target_include_directories(test_preprocessing PUBLIC include)
target_include_directories(test_preprocessing PUBLIC .)
target_link_libraries(test_preprocessing ${GTEST_LIBRARIES} gtest_main pthread preprocessing)
