# Set specific directories
set(DPC_HYSTERESIS_DIR "${CMAKE_INSTALL_PREFIX}/dpc_hysteresis-${HYST_LIB_VERSION}" CACHE STRING "top installation directory")
set(DPC_HYSTERESIS_INCLUDE_DIR "${DPC_HYSTERESIS_DIR}/include")
set(DPC_HYSTERESIS_LIB_DIR "${DPC_HYSTERESIS_DIR}/lib")

# Search for required includes and libraries

find_package(Boost COMPONENTS thread REQUIRED)
message(STATUS "Found boost_thread")

find_package(Eigen3 REQUIRED)
message(STATUS "Found Eigen3 in: ${EIGEN3_INCLUDE_DIR}")

#find_package(CGAL REQUIRED)

#find_package(cxx_misc_tools REQUIRED)

# Language Standard and compiler options

# Use C++20
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Compiler options
add_compile_options(-O3)





#include_directories(${CMAKE_SOURCE_DIR}/main_lib/include)

set(SOURCES continuous_dpc_model.cpp
            dpc_hysteresis_model.cpp
            elementary_hysteron.cpp
            malloc_grid.cpp
            xml_grid_constructor.cpp
            critical_surface/shapes/CS_sphere_gen.cpp
            critical_surface/shapes/CS_sphere.cpp
            hystLib_functions.cpp
            third_party/tinyxml2.cpp
    )

set(HYST_HEADER_PATH ${CMAKE_SOURCE_DIR}/main_lib/include)

# target, includes and linked libraries
add_library(${PROJECT_NAME} SHARED ${SOURCES})
target_include_directories(${PROJECT_NAME} PUBLIC   $<BUILD_INTERFACE:${HYST_HEADER_PATH}>
                                                    $<INSTALL_INTERFACE:${DPC_HYSTERESIS_INCLUDE_DIR}>
                                                ${EIGEN3_INCLUDE_DIR}
                                                ${CGAL_INCLUDE_DIR})

target_link_libraries(${PROJECT_NAME} PRIVATE boost_thread 
                                            boost_system 
                                            boost_filesystem)          

# installation and packaging 
install(TARGETS ${PROJECT_NAME} EXPORT ${PROJECT_NAME} DESTINATION ${DPC_HYSTERESIS_LIB_DIR})
install(DIRECTORY  ${HYST_HEADER_PATH} DESTINATION ${DPC_HYSTERESIS_DIR})
install(FILES ../dpc_hysteresis-config.cmake DESTINATION ${DPC_HYSTERESIS_DIR})
install(EXPORT ${PROJECT_NAME} DESTINATION "${DPC_HYSTERESIS_DIR}")