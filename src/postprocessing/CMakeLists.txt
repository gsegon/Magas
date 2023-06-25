

add_library(postprocessor OBJECT
        ExportVtu.cpp
        MatIDPostprocessor.cpp
        ExpressionPostprocessor.cpp
        ExpressionPostprocessorSum.cpp)

target_include_directories(postprocessor PUBLIC include)
target_compile_options(postprocessor PRIVATE -Wall -Wextra)
set(deallii_dir "/home/gordan/Programs/dealii/dealii-install")
target_include_directories(postprocessor PUBLIC ${deallii_dir}/include/)
target_include_directories(postprocessor PUBLIC /home/gordan/Programs/vtu11)
target_include_directories(postprocessor PUBLIC /home/gordan/Programs/exprtk)

# test export
find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})

add_executable(test_postprocessor
        test/test_export/test_ExportVtu.cpp
        test/test_processors/test_ExpressionPostprocessor.cpp
        test/test_processors/test_ExpressionPostprocessorSum.cpp)
target_link_libraries(test_postprocessor ${GTEST_LIBRARIES} gtest_main pthread solver postprocessor)

add_test(test_postprocessor_gtest test_postprocessor)