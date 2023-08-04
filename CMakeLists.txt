cmake_minimum_required(VERSION 3.16)
project(Magas VERSION 0.3.0 LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)

enable_testing()

# Modules
add_subdirectory(src/material)
add_subdirectory(src/sources)
add_subdirectory(src/input)
add_subdirectory(src/preprocessing)
add_subdirectory(src/postprocessing)
add_subdirectory(src/solver)

# Apps
add_subdirectory(src)

# Examples
add_subdirectory(examples)

# install
install(DIRECTORY examples DESTINATION $ENV{HOME}/magas/ PATTERN "CMakeLists.txt" EXCLUDE PERMISSIONS GROUP_EXECUTE GROUP_READ GROUP_WRITE GROUP_READ GROUP_WRITE WORLD_EXECUTE WORLD_READ WORLD_WRITE)
install(TARGETS magas)

# package
set(CPACK_PACKAGE_VENDOR "Gordan Segon")
set(CPACK_PACKAGE_CONTACT "gordan.segon@gmail.com")
set(CPACK_PACKAGE_DESCRIPTION "Magnetostatic analysis suite")
set(CPACK_DEBIAN_PACKAGE_DEPENDS "libdeal.ii-dev (>=9.4.0)")
include(CPack)
