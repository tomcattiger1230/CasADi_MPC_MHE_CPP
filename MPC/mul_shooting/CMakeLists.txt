include_directories(include )

add_executable(multi_shooting src/ca_multi_shooting.cpp)
target_link_libraries(multi_shooting casadi)

