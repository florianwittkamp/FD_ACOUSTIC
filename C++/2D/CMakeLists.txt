cmake_minimum_required(VERSION 2.8.9)
project(FD_2D_DX4_DT2)
 
file(GLOB SOURCES "*.cpp")
file(GLOB HEADER_FILES "*.h")
 
add_executable(FD_2D_DX4_DT2 ${SOURCES})