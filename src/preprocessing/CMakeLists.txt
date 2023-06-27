add_library(preprocessing OBJECT
        Dummy.cpp
        include/IPeriodicityMapper.h
        include/LinePeriodicityMapper.h
        include/CirclePeriodicityMapper.h include/ConfigParser.h)

target_include_directories(preprocessing PUBLIC include)
target_include_directories(preprocessing PUBLIC ~/Programs/vtu11)
target_include_directories(preprocessing PUBLIC ~/Programs/json/include)
target_link_libraries(preprocessing optimized
        tbb)
target_compile_features(preprocessing PRIVATE cxx_std_17)

# test preprocessing
find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})
add_executable(test_preprocessing
#        test/test_gmsh_mesh.cpp
        test/test_LinePeriodicityMapper.cpp
        test/test_CirclePeriodicityMapper.cpp test/test_ConfigParser.cpp)
target_include_directories(test_preprocessing PUBLIC include)
target_include_directories(test_preprocessing PUBLIC .)
target_link_libraries(test_preprocessing ${GTEST_LIBRARIES} gtest_main pthread preprocessing)
target_compile_features(test_preprocessing PRIVATE cxx_std_17)

add_test(test_preprocessing_gtest test_preprocessing)
