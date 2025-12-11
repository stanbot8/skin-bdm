# -----------------------------------------------------------------------------
#
# Copyright (C) 2021 CERN & University of Surrey for the benefit of the
# BioDynaMo collaboration. All Rights Reserved.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
#
# See the LICENSE file distributed with this work for details.
# See the NOTICE file distributed with this work for additional information
# regarding copyright ownership.
#
# -----------------------------------------------------------------------------
cmake_minimum_required(VERSION 3.19.3)
project(skibidy)

# BioDynaMo curretly uses the C++17 standard.
set(CMAKE_CXX_STANDARD 17)

# Use BioDynaMo in this project.
find_package(BioDynaMo REQUIRED)

# See UseBioDynaMo.cmake in your BioDynaMo build folder for details.
# Note that BioDynaMo provides gtest header/libraries in its include/lib dir.
include(${BDM_USE_FILE})

# Consider all files in src/ and modules/ for BioDynaMo simulation.
include_directories("src")
include_directories("modules")
include_directories("third_party")
file(GLOB_RECURSE PROJECT_HEADERS src/*.h modules/*.h)
file(GLOB_RECURSE PROJECT_SOURCES src/*.cc)

# main.cc lives at project root (outside src/) so libskibidy.so has no
# main() symbol -- the test binary provides its own gtest main instead.
bdm_add_executable(${CMAKE_PROJECT_NAME}
                   HEADERS ${PROJECT_HEADERS}
                   SOURCES ${PROJECT_SOURCES} ${CMAKE_SOURCE_DIR}/main.cc
                   LIBRARIES ${BDM_REQUIRED_LIBRARIES})

# Consider all files in tests/ for GoogleTests.
# Note: we avoid bdm_add_test / gtest_discover_tests because BioDynaMo's
# CommandLineOptions parser rejects --gtest_list_tests at discovery time.
# Using gtest_add_tests (source-based discovery) instead.
include_directories("tests")
file(GLOB_RECURSE TEST_SOURCES tests/*.cc)
file(GLOB_RECURSE TEST_HEADERS tests/*.h)

add_library(libgtest IMPORTED STATIC GLOBAL)
set_target_properties(libgtest PROPERTIES
    IMPORTED_LOCATION "$ENV{BDMSYS}/lib/libgtest.a"
    IMPORTED_LINK_INTERFACE_LIBRARIES "${CMAKE_THREAD_LIBS_INIT}")

add_executable(${CMAKE_PROJECT_NAME}-test ${TEST_SOURCES})
target_link_libraries(${CMAKE_PROJECT_NAME}-test
                      ${BDM_REQUIRED_LIBRARIES} ${CMAKE_PROJECT_NAME} libgtest)
# Register as a single test -- BioDynaMo's CommandLineOptions parser rejects
# --gtest_filter and --gtest_list_tests, so per-test CTest entries don't work.
add_test(NAME ${CMAKE_PROJECT_NAME}-tests
         COMMAND ${CMAKE_PROJECT_NAME}-test)
