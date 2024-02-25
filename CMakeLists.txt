cmake_minimum_required(VERSION 3.26)
project(interpolation_project)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -fopenmp")

add_executable(interpolation_project main.cpp)

include_directories(B_spline_lib)
add_subdirectory(B_spline_lib)

target_link_libraries(interpolation_project B_spline_lib)

#enable_testing()
#add_subdirectory(ctest)