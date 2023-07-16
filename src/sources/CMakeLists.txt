add_library(sources OBJECT
        src/ConstFSource.cpp src/ExprFSource.cpp src/FSourceFactory.cpp)

target_include_directories(sources PUBLIC include)
target_include_directories(sources PRIVATE ../../external/exprtk)

# test material
find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})

add_executable(test_sources
        test/test_ConstFSource.cpp
        test/test_ExprFSource.cpp
        test/test_FSourceFactory.cpp)

target_link_libraries(test_sources ${GTEST_LIBRARIES} gtest_main pthread sources)
target_include_directories(test_sources PRIVATE include)
target_include_directories(test_sources PRIVATE ../../external/exprtk)


# CMake test
add_test(test_sources_gtest test_sources)
