include_directories(include )
LINK_DIRECTORIES(lib)
add_executable(multi_shooting src/ca_multi_shooting.cpp)
target_link_libraries(multi_shooting casadi )

# add_executable(multi_shooting_lib src/ca_multi_shooting_lib.cpp)
# target_link_libraries(multi_shooting_lib casadi ${CMAKE_CURRENT_SOURCE_DIR}/lib/robot_mpc.so)