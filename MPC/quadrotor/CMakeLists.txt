include_directories(include)

set(CMAKE_BUILD_TYPE "Release")
find_package(PkgConfig)
pkg_search_module(Eigen3 REQUIRED eigen3)
include_directories(${Eigen3_INCLUDE_DIRS})
pkg_search_module(casadi REQUIRED casadi)
include_directories(${casadi_INCLUDE_DIRS})

# LINK_DIRECTORIES(lib)
add_library(dmoc_quadrotor_lib src/model/dmoc_quadrotor.cpp)
add_library(mpc_quadrotor_lib src/model/mpc_quadrotor.cpp)

add_executable(quadrotor_dmoc src/dmoc_quadrotor_controller.cpp)
target_link_libraries(quadrotor_dmoc dmoc_quadrotor_lib casadi)

add_executable(quadrotor_mpc src/mpc_quadrotor_controller.cpp)
target_link_libraries(quadrotor_mpc mpc_quadrotor_lib casadi)