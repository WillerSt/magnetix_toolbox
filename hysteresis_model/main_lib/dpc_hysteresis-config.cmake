# my_library-config.cmake - package configuration file

get_filename_component(SELF_DIR "${CMAKE_CURRENT_LIST_FILE}" PATH)
include(${SELF_DIR}/DPC_Hysteresis.cmake)

#find_package(cxx_misc_tools REQUIRED)