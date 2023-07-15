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

DEAL_II_SETUP_TARGET(magas)

target_compile_options(magas PUBLIC -Wall -Wextra)
target_include_directories(magas PUBLIC ~/Programs/cxxopts/include)
target_include_directories(magas PUBLIC ~/Programs/exprtk)
target_include_directories(magas PUBLIC ~/Programs/json/include)
