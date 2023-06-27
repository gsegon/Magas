# aml
add_executable(aml AppMagnetostaticLinear.cpp)
target_link_libraries(aml optimized
        tbb
        solver
        postprocessor)

set(DEAL_II_DIR ~/Programs/dealii/dealii/build)
find_package(deal.II 9.4.0 REQUIRED HINTS ${DEAL_II_DIR} ../ ../../ $ENV{DEAL_II_DIR})
DEAL_II_INITIALIZE_CACHED_VARIABLES()
DEAL_II_SETUP_TARGET(aml)

target_compile_options(aml PRIVATE -Wall -Wextra)
target_include_directories(aml PUBLIC ~/Programs/cxxopts/include)
target_include_directories(aml PUBLIC ~/Programs/exprtk)
target_include_directories(aml PUBLIC ~/Programs/json/include)

# amp
add_executable(amp AppMagnetostaticPicard.cpp)
target_link_libraries(amp optimized
        tbb
        solver
        postprocessor)

DEAL_II_SETUP_TARGET(amp)

target_compile_options(amp PRIVATE -Wall -Wextra)
target_include_directories(amp PUBLIC ~/Programs/cxxopts/include)
target_include_directories(amp PUBLIC ~/Programs/exprtk)
target_include_directories(amp PUBLIC ~/Programs/json/include)

