add_library(input OBJECT
        Dummy.cpp)

target_include_directories(input PUBLIC include)
target_include_directories(input PUBLIC ~/Programs/json/include)

# test material
find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})

add_executable(test_input
        test/test_ConfigParser.cpp)
target_link_libraries(test_input ${GTEST_LIBRARIES} gtest_main pthread input)
target_include_directories(test_input PUBLIC include)

# CMake test
add_test(test_input_gtest test_input)
