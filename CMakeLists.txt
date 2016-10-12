cmake_minimum_required(VERSION 3.5)

project(Oersted)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O0 --coverage")

include_directories(./lib/ ./lib/Eigen/ ./lib/GoogleTest/googletest/include/gtest)
add_subdirectory(./lib/)

include_directories(./src/ ./src/Sketch/include)
add_subdirectory(./src/)

add_subdirectory(./test/)