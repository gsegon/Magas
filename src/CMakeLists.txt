# aml
add_executable(aml AppMagnetostaticLinear.cpp)
target_link_libraries(aml optimized
        /home/gordan/Programs/dealii/dealii-install/lib/libdeal_II.so
        tbb
        solver
        postprocessor)
target_include_directories(aml PUBLIC /home/gordan/Programs/cxxopts/include)
target_include_directories(aml PUBLIC /home/gordan/Programs/exprtk)
target_include_directories(aml PUBLIC "/home/gordan/Programs/json/include")

# amp
add_executable(amp AppMagnetostaticPicard.cpp)
target_link_libraries(amp optimized
        /home/gordan/Programs/dealii/dealii-install/lib/libdeal_II.so
        tbb
        solver
        postprocessor)
target_include_directories(amp PUBLIC /home/gordan/Programs/cxxopts/include)
target_include_directories(amp PUBLIC /home/gordan/Programs/exprtk)
target_include_directories(amp PUBLIC "/home/gordan/Programs/json/include")