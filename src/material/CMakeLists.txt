add_library(material OBJECT
        src/LinearBHCurve.cpp
        src/InterpolatedBHCurve.cpp
        src/AnalyticBHCurve.cpp
        src/BHCurveFactory.cpp)

target_include_directories(material PUBLIC include)
target_include_directories(material PRIVATE ../../external/strtk)
target_link_libraries(material PRIVATE gsl)

# test material
find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})

add_executable(test_material
        test/test_InterpolatedBHCurve.cpp
        test/test_AnalyticBHCurve.cpp
        test/test_LinearBHCurve.cpp
        test/test_BHCurveFactory.cpp)
target_link_libraries(test_material ${GTEST_LIBRARIES} gtest_main pthread material)
target_include_directories(test_material PRIVATE include)
target_include_directories(test_material PRIVATE ../../external/strtk)

# CMake test
add_test(test_material_gtest test_material)
