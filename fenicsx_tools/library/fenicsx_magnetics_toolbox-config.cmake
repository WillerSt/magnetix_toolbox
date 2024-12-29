# my_library-config.cmake - package configuration file

get_filename_component(SELF_DIR "${CMAKE_CURRENT_LIST_FILE}" PATH)
include(${SELF_DIR}/fenicsx_magnetics_toolbox.cmake)

find_package(DOLFINX REQUIRED)
find_package(DPC_Hysteresis REQUIRED)
find_package(DOLFINX_MPC)
