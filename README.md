This repository is for research work within the framework of the final qualifying work in the 2023/2024 academic year.

## To run the code:
1. clone this repository
2. Make sure you have MinGW (gcc and g++ compilers) and CMake installed and paths set
3. In project root create Build directory `mkdir Build`
4. Move to the newly created directory `cd Build`
5. Generate Makefile: `cmake ../ -G"MinGW Makefiles"`
6. Build the project: `cmake --build .`
7. Run the code: `interpolation_project.exe`
8. Results will be shown in CLI (Command line interface)