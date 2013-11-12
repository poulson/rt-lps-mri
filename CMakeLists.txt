#
#  Copyright (c) 2013, Jack Poulson, Ricardo Otazo, and Emmanuel Candes
#  All rights reserved.
#
#  This file is part of RealTime-MRI and is under the BSD 2-Clause License, 
#  which can be found in the LICENSE file in the root directory, or at 
#  http://opensource.org/licenses/BSD-2-Clause
#
cmake_minimum_required(VERSION 2.8.8) 
project(RT-MRI)

set(RTMRI_VERSION_MAJOR 0)
set(RTMRI_VERSION_MINOR 1)

option(RTMRI_TESTS "Build a collection of test executables" ON)
option(RTMRI_EXAMPLES "Build a few example drivers" ON)

add_subdirectory(${PROJECT_SOURCE_DIR}/external/elemental)
include_directories(${PROJECT_BINARY_DIR}/external/elemental/include)
include_directories(${MPI_CXX_INCLUDE_PATH})
set(ELEM_CMAKE_DIR ${PROJECT_SOURCE_DIR}/external/elemental/cmake)
if(HAVE_QT5)
  include_directories(${Qt5Widgets_INCLUDE_DIRS})
  add_definitions(${Qt5Widgets_DEFINITIONS})
endif()

# TODO: Test for libnfft3.a

# Create the RT-MRI configuration header
configure_file( 
  ${PROJECT_SOURCE_DIR}/cmake/config.h.cmake
  ${PROJECT_BINARY_DIR}/include/RT-MRI/config.h)
install(FILES ${PROJECT_BINARY_DIR}/include/RT-MRI/config.h 
  DESTINATION include/RT-MRI)

# TODO: rtmrivariables Makefile include

# Grab all of the .c, .cpp, .h, and .hpp RT-MRI files
file(GLOB_RECURSE RTMRI_CPP RELATIVE ${PROJECT_SOURCE_DIR} "src/*.cpp")
file(GLOB_RECURSE RTMRI_HEADERS RELATIVE ${PROJECT_SOURCE_DIR} 
  "include/*.h" "include/*.hpp")
set(RTMRI_SRC "${RTMRI_CPP};${RTMRI_HEADERS}")

# The main library
add_library(rtmri ${LIBRARY_TYPE} ${RTMRI_SRC})
target_link_libraries(rtmri elemental)
endif()
install(TARGETS rtmri DESTINATION lib)

# Define the header-file preparation rules
set(PREPARED_HEADERS)
foreach(HEADER ${CLIQUE_HEADERS})
  add_custom_command(OUTPUT ${PROJECT_BINARY_DIR}/${HEADER}
    COMMAND ${CMAKE_COMMAND} -E copy ${PROJECT_SOURCE_DIR}/${HEADER}
            ${PROJECT_BINARY_DIR}/${HEADER}
    DEPENDS "${PROJECT_SOURCE_DIR}/${HEADER}")
  list(APPEND PREPARED_HEADERS ${PROJECT_BINARY_DIR}/${HEADER})

  get_filename_component(HEADER_PATH ${HEADER} PATH)
  install(FILES ${PROJECT_BINARY_DIR}/${HEADER} DESTINATION ${HEADER_PATH})
endforeach()
add_custom_target(prepare_rtmri_headers DEPENDS ${PREPARED_HEADERS})
add_dependencies(rtmri prepare_rtmri_headers)

# Make sure the RT-MRI headers can be found
include_directories("${PROJECT_BINARY_DIR}/include")

# Build the test drivers if necessary
if(CLIQ_TESTS)
  set(TEST_DIR ${PROJECT_SOURCE_DIR}/tests)
  set(TESTS)

  # Build the tests
  set(OUTPUT_DIR "${PROJECT_BINARY_DIR}/bin/tests")
  foreach(TEST ${TESTS})
    add_executable(tests-${TEST} ${TEST_DIR}/${TEST}.cpp)
    target_link_libraries(tests-${TEST} rtmri)
    set_target_properties(tests-${TEST} PROPERTIES
      OUTPUT_NAME ${TEST} RUNTIME_OUTPUT_DIRECTORY ${OUTPUT_DIR})
    if(MPI_LINK_FLAGS)
      set_target_properties(tests-${TEST} PROPERTIES
        LINK_FLAGS ${MPI_LINK_FLAGS})
    endif()
    install(TARGETS tests-${TEST} DESTINATION bin/tests)
  endforeach()
endif()

# Build the example drivers if necessary
if(RTMRI_EXAMPLES)
  set(EXAMPLE_DIR ${PROJECT_SOURCE_DIR}/examples)
  set(EXAMPLES Version)

  # Build the examples
  set(OUTPUT_DIR "${PROJECT_BINARY_DIR}/bin/examples")
  foreach(EXAMPLE ${EXAMPLES})
    add_executable(examples-${EXAMPLE} ${EXAMPLE_DIR}/${EXAMPLE}.cpp)
    target_link_libraries(examples-${EXAMPLE} rtmri)
    set_target_properties(examples-${EXAMPLE} PROPERTIES
      OUTPUT_NAME ${EXAMPLE} RUNTIME_OUTPUT_DIRECTORY ${OUTPUT_DIR})
    if(MPI_LINK_FLAGS)
      set_target_properties(examples-${EXAMPLE} PROPERTIES
        LINK_FLAGS ${MPI_LINK_FLAGS})
    endif()
    install(TARGETS examples-${EXAMPLE} DESTINATION bin/examples)
  endforeach()
endif()

# If RT-MRI is a subproject, then pass some variables to the parent
if(NOT CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
  include(./cmake/rtmri_sub.cmake)
endif()
