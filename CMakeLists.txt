cmake_minimum_required(VERSION 3.2)

project(Oersted)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++14")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O0 --coverage")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

find_package(Boost REQUIRED COMPONENTS system filesystem)
include_directories(${Boost_INCLUDE_DIR})

include_directories(./lib/)
include_directories(./lib/Eigen/)
include_directories(./lib/GoogleTest/googletest/include/gtest)

add_subdirectory(./lib/)

include_directories(./src/)
include_directories(./src/Sketch/include/)
include_directories(./src/Mesh/include/)

add_subdirectory(./src/)

add_subdirectory(./test/)
