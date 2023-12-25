cmake_minimum_required(VERSION 3.26)
project(interpolation_project)

set(CMAKE_CXX_STANDARD 17)

add_executable(interpolation_project main.cpp
        gnuplot-iostream.h
)
