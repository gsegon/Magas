# magas
add_executable(magas AppMagnetostaticLinear.cpp)
target_link_libraries(magas optimized
        tbb
        solver
        postprocessor
        preprocessing
        material
        sources)

set(DEAL_II_DIR ~/Programs/dealii/dealii/build)
find_package(deal.II 9.4.0 REQUIRED HINTS ${DEAL_II_DIR} ../ ../../ $ENV{DEAL_II_DIR})
DEAL_II_INITIALIZE_CACHED_VARIABLES()
DEAL_II_SETUP_TARGET(magas)

target_compile_options(magas PRIVATE -Wall -Wextra)
target_include_directories(magas PUBLIC ~/Programs/cxxopts/include)
target_include_directories(magas PUBLIC ~/Programs/exprtk)
target_include_directories(magas PUBLIC ~/Programs/json/include)


# amnl
add_executable(amnl AppMagnetostaticNonLinear.cpp)
target_link_libraries(amnl optimized
        tbb
        solver
        postprocessor
        preprocessing
        material
        sources)

set(DEAL_II_DIR ~/Programs/dealii/dealii/build)
find_package(deal.II 9.4.0 REQUIRED HINTS ${DEAL_II_DIR} ../ ../../ $ENV{DEAL_II_DIR})
DEAL_II_INITIALIZE_CACHED_VARIABLES()
DEAL_II_SETUP_TARGET(amnl)

target_compile_options(amnl PRIVATE -Wall -Wextra)
target_include_directories(amnl PUBLIC ~/Programs/cxxopts/include)
target_include_directories(amnl PUBLIC ~/Programs/exprtk)
target_include_directories(amnl PUBLIC ~/Programs/json/include)
#
## amp
#add_executable(amp AppMagnetostaticPicard.cpp)
#target_link_libraries(amp optimized
#        tbb
#        solver
#        postprocessor)
#
#DEAL_II_SETUP_TARGET(amp)
#
#target_compile_options(amp PRIVATE -Wall -Wextra)
#target_include_directories(amp PUBLIC ~/Programs/cxxopts/include)
#target_include_directories(amp PUBLIC ~/Programs/exprtk)
#target_include_directories(amp PUBLIC ~/Programs/json/include)

