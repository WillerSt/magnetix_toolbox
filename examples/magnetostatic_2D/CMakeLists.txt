# SPDX-FileCopyrightText: 2025 Stephan Willerich
# SPDX-License-Identifier: MIT License

cmake_minimum_required(VERSION 3.16)

set(PROJECT_NAME magnetostatic_2D_exec)
project(${PROJECT_NAME} LANGUAGES C CXX)

# Set C++20 standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
add_compile_options(-O3)

if(NOT TARGET dolfinx)
  find_package(DOLFINX REQUIRED)
endif()

find_package(fenicsx_magnetics_toolbox REQUIRED)

include(CheckSymbolExists)
set(CMAKE_REQUIRED_INCLUDES ${PETSC_INCLUDE_DIRS})
check_symbol_exists(PETSC_USE_COMPLEX petscsystypes.h PETSC_SCALAR_COMPLEX)

# Add target to compile UFL files
if(PETSC_SCALAR_COMPLEX EQUAL 1)
  set(SCALAR_TYPE "--scalar_type=double _Complex")
endif()

set(CMAKE_INCLUDE_CURRENT_DIR ON)

add_executable(${PROJECT_NAME}	main.cpp)
								
target_link_libraries(${PROJECT_NAME} PUBLIC fenicsx_magnetics_toolbox)


include(CheckCXXCompilerFlag)
check_cxx_compiler_flag("-Wno-comment" HAVE_NO_MULTLINE)

target_compile_definitions(${PROJECT_NAME} PUBLIC GLIBCXX_USE_CXX20_ABI)


