cmake_minimum_required(VERSION 3.21)

set(PROJECT_NAME FVM)

project(PROJECT_NAME)

set(SOURCES
    main.cpp

    Utils/MatrixSolver.cpp
    Utils/Timer.cpp

    Mesh/MeshBase.cpp
    Mesh/Geometry.cpp
    Mesh/2D/Structured/CartesianMesh2D.cpp
    Mesh/2D/Unstructured/PolyMesh2D.cpp
    Mesh/2D/Unstructured/Parse.cpp

    Discretization/LinearCombination/LinearCombination.hpp
    Discretization/Interpolation/Interpolation.hpp

    Solvers/SIMPLE/SimpleAlgorithm.cpp
)

set(HEADERS
    Utils/Types.h
    Utils/TypesOperations.h
    Utils/MatrixSolver.h
    Utils/Timer.h

    Config/PhysicalProperties.h
    Config/SolverControl.h

    Boundary/BoundaryCondition.h

    Mesh/MeshBase.h
    Mesh/Geometry.h
    Mesh/2D/Structured/CartesianMesh2D.h
    Mesh/2D/Unstructured/PolyMesh2D.h
    Mesh/2D/Unstructured/Parse.h

    Discretization/LinearCombination/LinearCombination.h
    Discretization/Interpolation/Interpolation.h

    Solvers/SIMPLE/SimpleAlgorithm.h
)

set(CMAKE_CXX_STANDARD 20)

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -Wextra -Wpedantic")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -march=native")

set(TARGET_NAME fvm_test)
add_executable(${TARGET_NAME} ${SOURCES} ${HEADERS})

target_include_directories(${TARGET_NAME} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

# OpenMP
find_package(OpenMP)
if (OpenMP_CXX_FOUND)
    target_link_libraries(${TARGET_NAME} OpenMP::OpenMP_CXX)
else()
    message(WARNING "CMake can't find OpenMP\n")
endif()

#Eigen
add_subdirectory(Eigen)
target_link_libraries(${TARGET_NAME} eigen)
