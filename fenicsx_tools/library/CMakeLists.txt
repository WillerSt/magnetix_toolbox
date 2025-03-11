# SPDX-FileCopyrightText: 2025 Stephan Willerich
# SPDX-License-Identifier: MIT License

cmake_minimum_required(VERSION 3.16)

set(PROJECT_NAME fenicsx_magnetics_toolbox)
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
project(${PROJECT_NAME} LANGUAGES C CXX)

#-----------------------------------------------------------------------------------
# set install directory

if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    message(STATUS "Setting default install path")
    SET(CMAKE_INSTALL_PREFIX "${CMAKE_CURRENT_LIST_DIR}/../../install" CACHE PATH "install prefix" FORCE)
endif()

message(STATUS "installation directory=${CMAKE_INSTALL_PREFIX}")

#-----------------------------------------------------------------------------------
# define library version 

set(MAG_TOOLS_VERSION 0 CACHE STRING "major version" FORCE)
set(MAG_TOOLS_VERSION_MINOR 9 CACHE STRING "minor version" FORCE)
set(MAG_TOOLS_VERSION ${MAG_TOOLS_VERSION}.${MAG_TOOLS_VERSION_MINOR} CACHE STRING "version" FORCE)


# Set specific directories
set(MAG_TOOLS_DIR "${CMAKE_INSTALL_PREFIX}/fenicsx_magnetics_toolbox-${MAG_TOOLS_VERSION}" CACHE STRING "top installation directory")
set(MAG_TOOLS_INCLUDE_DIR "${MAG_TOOLS_DIR}/include")
set(MAG_TOOLS_FFCX_DIR "${MAG_TOOLS_DIR}/ffcx_code")
set(MAG_TOOLS_FFCX_INCLUDE_DIR "${MAG_TOOLS_FFCX_DIR}/include")
set(MAG_TOOLS_LIB_DIR "${MAG_TOOLS_DIR}/lib")

# Set C++20 standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# required packages
find_package(DOLFINX REQUIRED)
find_package(DPC_Hysteresis REQUIRED)

# optional packages
find_package(DOLFINX_MPC QUIET)


include(CheckSymbolExists)
set(CMAKE_REQUIRED_INCLUDES ${PETSC_INCLUDE_DIRS})
check_symbol_exists(PETSC_USE_COMPLEX petscsystypes.h PETSC_SCALAR_COMPLEX)

if(PETSC_SCALAR_COMPLEX EQUAL 1)
  set(SCALAR_TYPE "--scalar_type=double _Complex")
endif()

set(SOURCES ffcx_code/src/A_form_NR.c 
			ffcx_code/src/curl_evaluation.c 
			ffcx_code/src/curl_projection.c 
			ffcx_code/src/quadVec_dgVec.c 			
			ffcx_code/src/curlcg2_dgVec.c 
			ffcx_code/src/quadSca_dgSca.c 
			ffcx_code/src/magForce_virtualWork.c
			ffcx_code/src/quadVec_dgVec_rotated.c
			src/file_management.cpp 
			src/scenario_description.cpp 
			src/maxwell_solvers.cpp
			src/quadrature_dof_coupler.cpp
			src/linear_solver.cpp
			src/function_processing.cpp
			src/post_processing.cpp
			src/field_sources.cpp
			src/field_quantity.cpp
			src/coupled_material_model_atan.cpp
			src/coupled_material_model_hgm.cpp
			src/coupled_material_model_spline.cpp
			src/input_interpreter.cpp
		)

set (PYTHON_FILES 	python-scripts/fmt_scen.py
					python-scripts/msh_to_xdmf.py
		)

if (DOLFINX_MPC_FOUND)
	message("-- Found optional package DOLFINX_MPC at " ${DOLFINX_MPC_DIR})
	set(SOURCES ${SOURCES}
			src/test_rotation_solver.cpp
		)
endif()

if (DPC_Hysteresis_FOUND)
	message("-- Found required package DPC_Hysteresis at " ${DPC_Hysteresis_DIR})
endif()


#add_executable(${PROJECT_NAME} main.cpp ${SOURCES})
add_library(${PROJECT_NAME} SHARED ${SOURCES})

target_include_directories(${PROJECT_NAME} PUBLIC 	$<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>
													$<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/ffcx_code/include>
													$<INSTALL_INTERFACE:${MAG_TOOLS_INCLUDE_DIR}>
													$<INSTALL_INTERFACE:${MAG_TOOLS_FFCX_INCLUDE_DIR}>)
								
target_link_libraries(${PROJECT_NAME} PUBLIC 	dolfinx
												DPC_Hysteresis
												boost_filesystem
												boost_program_options)

if (DOLFINX_MPC_FOUND)
	target_link_libraries(${PROJECT_NAME} PUBLIC 	dolfinx_mpc)
endif()

target_compile_definitions(${PROJECT_NAME} PUBLIC GLIBCXX_USE_CXX20_ABI)
target_compile_options(${PROJECT_NAME} PUBLIC -O3 -Wall -Wextra -pedantic)


# installation and packaging 
install(TARGETS ${PROJECT_NAME} EXPORT ${PROJECT_NAME} DESTINATION ${MAG_TOOLS_LIB_DIR})
install(DIRECTORY  ${CMAKE_SOURCE_DIR}/include DESTINATION ${MAG_TOOLS_DIR})
install(DIRECTORY  ${CMAKE_SOURCE_DIR}/ffcx_code/include DESTINATION ${MAG_TOOLS_FFCX_DIR})
install(FILES fenicsx_magnetics_toolbox-config.cmake DESTINATION ${MAG_TOOLS_DIR})
install(EXPORT ${PROJECT_NAME} DESTINATION "${MAG_TOOLS_DIR}")

install(FILES ${PYTHON_FILES} DESTINATION ${MAG_TOOLS_DIR}/python)



