cmake_minimum_required(VERSION 3.9.0)
project(synth_turb CXX)

add_library(synth_turb INTERFACE) # INTERFACE because it's header-only

# generate a header file with git revision id
if (EXISTS "${CMAKE_SOURCE_DIR}/../.git")
  execute_process(COMMAND
   bash -c "git log -1 --format=\"format:#define SYNTH_TURB_GIT_REVISION \\\"%H\\\"%n\" HEAD > git_revision.hpp"
   WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  )
endif()

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(synth_turb INTERFACE OpenMP::OpenMP_CXX)
endif()

#target_compile_options(synth_turb PUBLIC $<$<COMPILE_LANGUAGE:CXX>: ${OpenMP_CXX_FLAGS}>)
#target_link_libraries(synth_turb PUBLIC ${OpenMP_CXX_LIBRARIES})

target_include_directories(synth_turb INTERFACE
        $<INSTALL_INTERFACE:include>)
#        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>


# install files and export targets
#install(TARGETS cloudphxx_lgrngn EXPORT cloudphxx_lgrngn-targets
#  LIBRARY DESTINATION lib
#  INCLUDES DESTINATION include)

install(TARGETS synth_turb EXPORT synth_turb-targets
        INCLUDES DESTINATION include
        )

install(EXPORT synth_turb-targets DESTINATION share/synth_turb)

# generate and install a config file
include(CMakePackageConfigHelpers)

#export(TARGETS cloudphxx_lgrngn FILE MyLibraryConfig.cmake)

configure_package_config_file(
  synth_turb-config.cmake.in
  "${CMAKE_CURRENT_BINARY_DIR}/synth_turb-config.cmake"
  INSTALL_DESTINATION
    share/synth_turb
)
install(
  FILES
    "${CMAKE_CURRENT_BINARY_DIR}/synth_turb-config.cmake"
  DESTINATION
    share/synth_turb
)

install(
  DIRECTORY
    include/synth_turb
  DESTINATION
    include
)



#install(
#  FILES
#    synthetic_isotropic_turbulence_3d.hpp git_revision.hpp
#  DESTINATION
#    include/libmpdata++
#)
#
#install(
#  FILES
#    ../libmpdata++-config.cmake
#  DESTINATION
#    share/libmpdata++
#)
#
