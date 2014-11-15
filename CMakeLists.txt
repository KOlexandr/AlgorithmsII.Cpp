cmake_minimum_required(VERSION 2.8.4)
project(AlgorithmsII_Cpp)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -fopenmp")

set(SOURCE_FILES main.cpp)
add_executable(AlgorithmsII_Cpp ${SOURCE_FILES})