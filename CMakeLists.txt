cmake_minimum_required(VERSION 3.12)

set(MAIN_PROJECT_NAME LBM)

add_subdirectory(LBM_Core)


set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT LBM_Core)