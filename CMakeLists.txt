#
#  Copyright (c) 2013, Jack Poulson, Ricardo Otazo, and Emmanuel Candes
#  All rights reserved.
#
#  This file is part of Real-Time Low-rank Plus Sparse MRI (RT-LPS-MRI) and is 
#  under the BSD 2-Clause License, which can be found in the LICENSE file in the
#  root directory, or at http://opensource.org/licenses/BSD-2-Clause
#
cmake_minimum_required(VERSION 2.8.8) 
project(RT-LPS-MRI)

set(RTLPSMRI_VERSION_MAJOR 0)
set(RTLPSMRI_VERSION_MINOR 1)

option(RTLPSMRI_TESTS "Build a collection of test executables" ON)
option(RTLPSMRI_EXAMPLES "Build a few example drivers" ON)

add_subdirectory(${PROJECT_SOURCE_DIR}/external/elemental)
include_directories(${PROJECT_BINARY_DIR}/external/elemental/include)
include_directories(${MPI_CXX_INCLUDE_PATH})
set(ELEM_CMAKE_DIR ${PROJECT_SOURCE_DIR}/external/elemental/cmake)
if(HAVE_QT5)
  include_directories(${Qt5Widgets_INCLUDE_DIRS})
  add_definitions(${Qt5Widgets_DEFINITIONS})
endif()

# Test whether or not MATH_LIBS contains NFFT3 (and FFTW3?)
include(CheckFunctionExists)
set(CMAKE_REQUIRED_LIBRARIES ${MATH_LIBS})
check_function_exists(fftw_execute HAVE_FFTW_EXECUTE)
check_function_exists(nfft_trafo HAVE_NFFT_TRAFO)
if(NOT HAVE_FFTW_EXECUTE)
  message(FATAL_ERROR "MATH_LIBS did not contain FFTW3")
endif()
if(NOT HAVE_NFFT_TRAFO)
  message(FATAL_ERROR "MATH_LIBS did not contain NFFT3")
endif()

# Check if the include directories work for NFFT
set(NFFT_CODE
    "#include \"nfft3.h\"
     int main()
     {
       nfft_plan p;
       nfft_init_guru(&p, 2, NULL, 5, NULL, 2, 0, 0);
       nfft_precompute_one_psi(&p);
       nfft_trafo_2d(&p);
       nfft_adjoint_2d(&p);
       nfft_finalize(&p);
       return 0;
     }")
set(CMAKE_REQUIRED_INCLUDES ${NFFT_INC_DIR})
set(CMAKE_REQUIRED_LIBRARIES ${MATH_LIBS})
check_cxx_source_compiles("${NFFT_CODE}" HAVE_NFFT_COMPILE)
if(NOT HAVE_NFFT_COMPILE)
  message(FATAL_ERROR "Could not compile NFFT program using NFFT_INC_DIR=${NFFT_INC_DIR} and MATH_LIBS=${MATH_LIBS}")
endif()
include_directories(${NFFT_INC_DIR})

# Create the RT-LPS-MRI configuration header
configure_file( 
  ${PROJECT_SOURCE_DIR}/cmake/config.h.cmake
  ${PROJECT_BINARY_DIR}/include/rt-lps-mri/config.h)
install(FILES ${PROJECT_BINARY_DIR}/include/rt-lps-mri/config.h 
  DESTINATION include/rt-lps-mri)

# TODO: rtmrivariables Makefile include

# Grab all of the .c, .cpp, .h, and .hpp RT-LPS-MRI files
file(GLOB_RECURSE RTLPSMRI_CPP RELATIVE ${PROJECT_SOURCE_DIR} "src/*.cpp")
file(GLOB_RECURSE RTLPSMRI_HEADERS RELATIVE ${PROJECT_SOURCE_DIR} 
  "include/*.h" "include/*.hpp")
set(RTLPSMRI_SRC "${RTLPSMRI_CPP};${RTLPSMRI_HEADERS}")

# The main library
add_library(rtlpsmri ${LIBRARY_TYPE} ${RTLPSMRI_SRC})
target_link_libraries(rtlpsmri elemental)
install(TARGETS rtlpsmri DESTINATION lib)

# Define the header-file preparation rules
set(PREPARED_HEADERS)
foreach(HEADER ${RTLPSMRI_HEADERS})
  add_custom_command(OUTPUT ${PROJECT_BINARY_DIR}/${HEADER}
    COMMAND ${CMAKE_COMMAND} -E copy ${PROJECT_SOURCE_DIR}/${HEADER}
            ${PROJECT_BINARY_DIR}/${HEADER}
    DEPENDS "${PROJECT_SOURCE_DIR}/${HEADER}")
  list(APPEND PREPARED_HEADERS ${PROJECT_BINARY_DIR}/${HEADER})

  get_filename_component(HEADER_PATH ${HEADER} PATH)
  install(FILES ${PROJECT_BINARY_DIR}/${HEADER} DESTINATION ${HEADER_PATH})
endforeach()
add_custom_target(prepare_rtlpsmri_headers DEPENDS ${PREPARED_HEADERS})
add_dependencies(rtlpsmri prepare_rtlpsmri_headers)

# Make sure the RT-LPS-MRI headers can be found
include_directories("${PROJECT_BINARY_DIR}/include")

# Build the test drivers if necessary
if(RTLPSMRI_TESTS)
  set(TEST_DIR ${PROJECT_SOURCE_DIR}/tests)
  set(TESTS Acquisition CoilAwareNFFT NFFT)

  # Build the tests
  set(OUTPUT_DIR "${PROJECT_BINARY_DIR}/bin/tests")
  foreach(TEST ${TESTS})
    add_executable(tests-${TEST} ${TEST_DIR}/${TEST}.cpp)
    target_link_libraries(tests-${TEST} rtlpsmri)
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
if(RTLPSMRI_EXAMPLES)
  set(EXAMPLE_DIR ${PROJECT_SOURCE_DIR}/examples)
  set(EXAMPLES Version)

  # Build the examples
  set(OUTPUT_DIR "${PROJECT_BINARY_DIR}/bin/examples")
  foreach(EXAMPLE ${EXAMPLES})
    add_executable(examples-${EXAMPLE} ${EXAMPLE_DIR}/${EXAMPLE}.cpp)
    target_link_libraries(examples-${EXAMPLE} rtlpsmri)
    set_target_properties(examples-${EXAMPLE} PROPERTIES
      OUTPUT_NAME ${EXAMPLE} RUNTIME_OUTPUT_DIRECTORY ${OUTPUT_DIR})
    if(MPI_LINK_FLAGS)
      set_target_properties(examples-${EXAMPLE} PROPERTIES
        LINK_FLAGS ${MPI_LINK_FLAGS})
    endif()
    install(TARGETS examples-${EXAMPLE} DESTINATION bin/examples)
  endforeach()
endif()

# If RT-LPS-MRI is a subproject, then pass some variables to the parent
if(NOT CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
  include(./cmake/rtlpsmri_sub.cmake)
endif()
