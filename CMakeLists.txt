cmake_minimum_required(VERSION 3.0)

project(casadi_mpe_mpc)

set(CMAKE_CXX_STANDARD 17)

find_package(casadi)
find_package(PkgConfig)
pkg_search_module(Eigen3 REQUIRED eigen3)
include_directories(${Eigen3_INCLUDE_DIRS})

add_subdirectory(MPC)

