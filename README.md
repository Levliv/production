## Overview
This is a repository to support the experiments performed as part of the research final qualification work in 2023/2024 academic  year.

There are 2 blocks of conducted experiments:

1. experiments on the error of quasi-linear and linear interpolation;
2. experiments on parallelization of quasi-linear interpolation and time spent.
## To run the code:
1. clone this repository
2. Make sure you have MinGW (gcc and g++ compilers) and CMake installed and paths set
3. In project root create Build directory `mkdir Build`
4. Move to the newly created directory `cd Build`
5. Generate Makefile: `cmake ../ -G"MinGW Makefiles"`
6. Build the project: `cmake --build .`
7. Run the code: `interpolation_project.exe`
8. Results will be shown in CLI (Command line interface)