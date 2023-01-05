include_directories(include)

set(CMAKE_BUILD_TYPE "Release")
find_package(PkgConfig)
pkg_search_module(Eigen3 REQUIRED eigen3)
include_directories(${Eigen3_INCLUDE_DIRS})
pkg_search_module(casadi REQUIRED casadi)
include_directories(${casadi_INCLUDE_DIRS})

# LINK_DIRECTORIES(lib)
add_library(dmoc_quadrotor_lib src/dmoc_quadrotor.cpp)

add_executable(quadrotor_mpc src/quadrotor_controller.cpp )
target_link_libraries(quadrotor_mpc dmoc_quadrotor_lib casadi)