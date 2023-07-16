add_library(postprocessor OBJECT
        src/export/ExportVtu.cpp
        src/processors/MatIDPostprocessor.cpp
        src/processors/ExpressionCellPostprocessor.cpp
        src/processors/ExpressionScalarPostprocessor.cpp
        src/processors/ArkkioScalarPostprocessor.cpp
        src/processors/PointBabsScalarPostprocessor.cpp
        src/processors/ScalarPostprocessorFactory.cpp
        src/processors/FluxLinkageScalarPostprocessor.cpp
        src/processors/EggShellPostprocessor.cpp
        src/processors/TorqueEggShellCellPostprocessor.cpp
        src/processors/TorqueEggShellScalarPostprocessor.cpp
        src/processors/ForceEggShellScalarPostprocessor.cpp)

set(DEAL_II_DIR ~/Programs/dealii/dealii/build)
find_package(deal.II 9.4.0 REQUIRED HINTS ${DEAL_II_DIR} ../ ../../ $ENV{DEAL_II_DIR})
DEAL_II_INITIALIZE_CACHED_VARIABLES()
DEAL_II_SETUP_TARGET(postprocessor)

target_include_directories(postprocessor PUBLIC include)
target_compile_options(postprocessor PRIVATE -Wall -Wextra)
target_include_directories(postprocessor PRIVATE ../../external/vtu11)
target_include_directories(postprocessor PRIVATE ../../external/exprtk)
target_link_libraries(postprocessor PRIVATE material)
target_link_libraries(postprocessor PRIVATE sources)

# test export
find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})

add_executable(test_postprocessor
        test/test_export/test_ExportVtu.cpp
        test/test_processors/test_ExpressionCellPostprocessor.cpp
        test/test_processors/test_ExpressionScalarPostprocessor.cpp
        test/test_processors/test_ArkkioScalarPostprocessor.cpp
        test/test_processors/test_PointBabsScalarPostprocessor.cpp
        test/test_processors/test_ScalarPostprocessorFactory.cpp
        test/test_processors/test_FluxLinkageScalarPostprocessor.cpp
        test/test_processors/test_EggShellCellPostprocessor.cpp
        test/test_processors/test_TorqueEggShellPostprocessor.cpp
        test/test_processors/test_TorqueEggShellScalarPostprocessor.cpp
        test/test_processors/test_ForceEggShellScalarPostprocessor.cpp)
target_link_libraries(test_postprocessor ${GTEST_LIBRARIES} gtest_main pthread solver postprocessor material preprocessing sources)
DEAL_II_SETUP_TARGET(test_postprocessor)

add_test(test_postprocessor_gtest test_postprocessor)
