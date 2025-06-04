Library for simulating fluid flows using uniform and non uniform meshes with different algorithm approaches

# Building

## Prerequisites

* GCC 14/15 (MinGW or native)
or
* Clang 18/19
* CMake >= 3.20
* OpenMP >= 3.1
* [gmsh](https://gmsh.info) and python3 for building examples
* All other dependencies are fetched automatically

## CMake flags

* BUILD_EXAMPLES - add examples compilation
* BUILD_MESH_DRAWER - add tool for exemining mesh before running solver
By default all flags are disabled

example usage:

```
cmake .. -G 'Ninja' -DCMAKE_BUILD_TYPE=Release -DBUILD_EXAMPLES=1 -DBUILD_MESH_DRAWER=0
```

## Usage

```
add_executable(your_executable /* Your sources */)

add_subdirectory(FVM)

target_link_libraries(your_executable FVM)
```
