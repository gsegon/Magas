add_library(preprocessing OBJECT
        Dummy.cpp
        src/PeriodicityMapperFactory.cpp
        src/LinePeriodicityMapper.cpp
        src/CirclePeriodicityMapper.cpp
        src/SlidingRotation.cpp
        )

target_include_directories(preprocessing PUBLIC include)
target_include_directories(preprocessing PRIVATE ../external/vtu11)
target_include_directories(preprocessing PRIVATE ../external/json/include)
target_link_libraries(preprocessing optimized
        tbb
        gsl)

# test preprocessing
find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})
add_executable(test_preprocessing
        test/test_LinePeriodicityMapper.cpp
        test/test_CirclePeriodicityMapper.cpp
        test/test_PeriodicityMapperFactory.cpp
        test/test_SlidingRotation.cpp
        )
target_include_directories(test_preprocessing PUBLIC include)
target_include_directories(test_preprocessing PUBLIC .)
target_link_libraries(test_preprocessing ${GTEST_LIBRARIES} gtest_main pthread preprocessing)

add_test(test_preprocessing_gtest test_preprocessing)
