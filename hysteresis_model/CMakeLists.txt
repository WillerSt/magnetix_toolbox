CMAKE_MINIMUM_REQUIRED (VERSION 3.16)
set(PROJECT_NAME DPC_Hysteresis)
PROJECT(${PROJECT_NAME})

#-----------------------------------------------------------------------------------
# set install directory

if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    message(STATUS "Setting default install path")
    SET(CMAKE_INSTALL_PREFIX "${CMAKE_CURRENT_LIST_DIR}/../install" CACHE PATH "install prefix" FORCE)
endif()

message(STATUS "installation directory=${CMAKE_INSTALL_PREFIX}")

#-----------------------------------------------------------------------------------
# define library version 

set(HYST_MODELS_VERSION 0 CACHE STRING "major version" FORCE)
set(HYST_MODELS_VERSION_MINOR 1 CACHE STRING "minor version" FORCE)
set(HYST_LIB_VERSION ${HYST_MODELS_VERSION}.${HYST_MODELS_VERSION_MINOR} CACHE STRING "version" FORCE)

#-----------------------------------------------------------------------------------
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

add_subdirectory(main_lib/src)