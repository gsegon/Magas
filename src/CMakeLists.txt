# magas
add_executable(magas AppMagnetostaticSolver.cpp)
target_link_libraries(magas optimized
        tbb
        solver
        postprocessor
        preprocessing
        material
        sources
        input)

set(DEAL_II_DIR ~/Programs/dealii/dealii/build)
find_package(deal.II 9.4.0 REQUIRED HINTS ${DEAL_II_DIR} ../ ../../ $ENV{DEAL_II_DIR})
DEAL_II_INITIALIZE_CACHED_VARIABLES()
DEAL_II_SETUP_TARGET(magas)

target_compile_options(magas PRIVATE -Wall -Wextra)
target_include_directories(magas PUBLIC ~/Programs/cxxopts/include)
target_include_directories(magas PUBLIC ~/Programs/exprtk)
target_include_directories(magas PUBLIC ~/Programs/json/include)

add_custom_command(
        TARGET magas POST_BUILD
        DEPENDS magas
        COMMAND $<$<CONFIG:release>:${CMAKE_STRIP}>
        ARGS magas)