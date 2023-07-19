add_library(input OBJECT
        Dummy.cpp
        src/JsonInputTranslator.cpp)

target_include_directories(input PUBLIC include)
target_include_directories(input PUBLIC ../../external/json/include)
target_link_libraries(input PUBLIC material sources)

# test material
find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})

add_executable(test_input
        test/test_ConfigParser.cpp test/test_JsonInputTranslator.cpp)
target_link_libraries(test_input ${GTEST_LIBRARIES} gtest_main pthread input sources material)
target_include_directories(test_input PRIVATE include)

# CMake test
add_test(test_input_gtest test_input)
