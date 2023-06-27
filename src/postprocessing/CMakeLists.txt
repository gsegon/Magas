add_library(postprocessor OBJECT
        src/export/ExportVtu.cpp
        src/processors/MatIDPostprocessor.cpp
        src/processors/ExpressionCellPostprocessor.cpp
        src/processors/ExpressionScalarPostprocessor.cpp
        src/processors/ArkkioScalarPostprocessor.cpp
        src/processors/PointBabsScalarPostprocessor.cpp)

find_package(deal.II 9.4.0 REQUIRED)
DEAL_II_INITIALIZE_CACHED_VARIABLES()
DEAL_II_SETUP_TARGET(postprocessor)

target_include_directories(postprocessor PUBLIC include)
target_compile_options(postprocessor PRIVATE -Wall -Wextra)
target_include_directories(postprocessor PUBLIC ~/Programs/vtu11)
target_include_directories(postprocessor PUBLIC ~/Programs/exprtk)

# test export
find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})

add_executable(test_postprocessor
        test/test_export/test_ExportVtu.cpp
        test/test_processors/test_ExpressionCellPostprocessor.cpp
        test/test_processors/test_ExpressionScalarPostprocessor.cpp test/test_processors/test_ArkkioScalarPostprocessor.cpp test/test_processors/test_PointBabsScalarPostprocessor.cpp)
target_link_libraries(test_postprocessor ${GTEST_LIBRARIES} gtest_main pthread solver postprocessor)
DEAL_II_SETUP_TARGET(test_postprocessor)

add_test(test_postprocessor_gtest test_postprocessor)
