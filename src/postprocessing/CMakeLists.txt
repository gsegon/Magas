add_library(postprocessor OBJECT
        ExportVtu.cpp
        MagneticFluxPostprocessor.cpp
        MatIDPostprocessor.cpp
        MagneticEnergyPostprocessor.cpp
        MagneticEnergyDensityPostprocessor.cpp
        ExpressionPostprocessor.cpp
        ExpressionPostprocessorSum.cpp)

target_compile_options(postprocessor PRIVATE -Wall -Wextra)
target_include_directories(postprocessor PUBLIC include)
set(deallii_dir "/home/gordan/Programs/dealii/dealii-install")
target_include_directories(postprocessor PUBLIC ${deallii_dir}/include/)
target_include_directories(postprocessor PUBLIC /home/gordan/Programs/vtu11)
target_include_directories(postprocessor PUBLIC /home/gordan/Programs/exprtk)

# test export
find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})

add_executable(test_postprocessor
        test/test_ExportVtu.cpp
        test/test_MagneticFluxPostprocessor.cpp
        test/test_MagneticEnergyPostprocessor.cpp
        test/test_MagneticEnergyDensityPostprocessor.cpp
        test/test_ExpressionPostprocessor.cpp
        test/test_ExpressionPostprocessorSum.cpp)
target_link_libraries(test_postprocessor ${GTEST_LIBRARIES} gtest_main pthread solver postprocessor)

