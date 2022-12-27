include_directories(include)

set(CMAKE_BUILD_TYPE "Release")
find_package(PkgConfig)
pkg_search_module(Eigen3 REQUIRED eigen3)
include_directories(${Eigen3_INCLUDE_DIRS})
pkg_search_module(casadi REQUIRED casadi)
include_directories(${casadi_INCLUDE_DIRS})

# LINK_DIRECTORIES(lib)
add_library(dmcc_uav_manipulator src/models/dmcc_uav_manipulator.cpp)

add_executable(cpc_planner src/uav_manipulator_cpc_planner.cpp )
target_link_libraries(cpc_planner dmcc_uav_manipulator casadi)
