project(B_spline_lib)

set(HEADER_FILES
        Interpolant.h)

set(SOURCE_FILES
        Interpolant.cpp
)

add_library(B_spline_lib STATIC ${SOURCE_FILES} ${HEADER_FILES})