cmake_minimum_required(VERSION 3.9.0)
project(st_test CXX)
add_executable(periodicity_test periodicity_test.cpp)
add_executable(pair_separation_test pair_separation_test.cpp)

find_package(synth_turb)
target_link_libraries(periodicity_test synth_turb)
target_link_libraries(pair_separation_test synth_turb)
