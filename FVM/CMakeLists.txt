cmake_minimum_required(VERSION 3.21)

set(PROJECT_NAME FVM)

project(PROJECT_NAME)

set(SOURCES
    main.cpp

    Mesh/Geometry.cpp
    Mesh/2D/Structured/CartesianMesh2D.cpp

    Discretization/LinearCombination/LinearCombination.hpp
    Discretization/Interpolation/Interpolation.hpp

    Solvers/SIMPLE/SimpleAlgorithm.cpp
)

set(HEADERS
    Utils/Types.h
    Utils/TypesOperations.h

    Config/PhysicalProperties.h
    Config/SolverControl.h

    Boundary/BoundaryCondition.h

    Mesh/MeshBase.h
    Mesh/Geometry.h
    Mesh/2D/Structured/CartesianMesh2D.h

    Discretization/LinearCombination/LinearCombination.h
    Discretization/Interpolation/Interpolation.h

    Solvers/SIMPLE/SimpleAlgorithm.h
)

set(CMAKE_CXX_STANDARD 20)

set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -march=native -O2 -Wall -fopenmp")

set(TARGET_NAME fvm_test)
add_executable(${TARGET_NAME} ${SOURCES} ${HEADERS})


add_subdirectory(Eigen)

target_link_libraries(${TARGET_NAME} eigen)
target_include_directories(${TARGET_NAME} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
