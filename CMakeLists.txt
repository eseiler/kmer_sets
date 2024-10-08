# SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: CC0-1.0

cmake_minimum_required (VERSION 3.25)

# Define the application name and version.
project (k-mer_set_difference
         LANGUAGES CXX
         VERSION 1.0.0
         DESCRIPTION "My application description"
)

if (NOT DEFINED CMAKE_CXX_STANDARD)
    set (CMAKE_CXX_STANDARD 23)
endif ()

if (NOT DEFINED CMAKE_CXX_STANDARD_REQUIRED)
    set (CMAKE_CXX_STANDARD_REQUIRED OFF)
endif ()

if (NOT DEFINED CMAKE_CXX_EXTENSIONS)
    set (CMAKE_CXX_EXTENSIONS OFF)
endif ()

# This allows including `*.cmake` files from the `cmake` directory without specifying the full path.
list (APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")

# Specify the directories where to store the built archives, libraries and executables.
include (output_directories)

# Add packages.
# We use CPM for package management: https://github.com/cpm-cmake/CPM.cmake
# The added packages (e.g., hibf, sharg, seqan3) are defined in the `cmake/package-lock.cmake` file.
include (CPM)
CPMUsePackageLock (${CMAKE_CURRENT_LIST_DIR}/cmake/package-lock.cmake)

# Use ccache if available. This speeds up the build process by caching files that have been compiled before.
CPMGetPackage (use_ccache)
# CPMGetPackage (hibf)
CPMGetPackage (sharg)
CPMGetPackage (seqan3)
# CPMGetPackage (fmindex)

# Add the application. This will include `src/CMakeLists.txt`.
add_subdirectory (src)

# An option to disable configuring and building the tests. Tests are enabled by default.
# If your project-name (line 8 of this file) is `app-template`, the option will be `app-template_TEST`.
# It can be used when calling CMake: `cmake .. -Dapp-template_TEST=OFF`.
# It is good practice to allow disabling tests. If another project includes your application,
# it might not want to build your tests.
option (${PROJECT_NAME}_TEST "Enable testing for ${PROJECT_NAME}." OFF)

if (${PROJECT_NAME}_TEST)
    # Add the tests. This will include `test/CMakeLists.txt`.
    add_subdirectory (test EXCLUDE_FROM_ALL)
endif ()
