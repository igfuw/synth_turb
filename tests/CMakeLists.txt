cmake_minimum_required(VERSION 3.9.0)
project(st_test CXX)
add_executable(api_test api_test.cpp)
add_executable(dispersion_test dispersion_test.cpp)

find_package(synth_turb)
target_link_libraries(api_test synth_turb)
target_link_libraries(dispersion_test synth_turb)
