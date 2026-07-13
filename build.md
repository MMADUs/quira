# HOW TO BUILD THE PROJECT

NOTE: go the project root where CMakeLists.txt is

## build cmake command:
    cmake -S . -B build

-S .      source directory is current folder
-B build  put generated build files inside build/

## compile project inside build/ command:
    cmake --build build

## or build it directly using bash:
    chmod +x scripts/build.sh
then choose either:
- build: ./scripts/build.sh
- clean build: ./scripts/build.sh clean

## run tests command:
- regular: ctest --test-dir build
- more verbose output: ctest --test-dir build --output-on-failure

## run example:
- on linux/macos: ./build/examples/matrix_add
- on windows: .\build\examples\matrix_add.exe or .\build\examples\Debug\matrix_add.exe
