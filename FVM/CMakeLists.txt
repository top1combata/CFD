cmake_minimum_required(VERSION 3.21)

set(PROJECT_NAME FVM)

set(SOURCES
    main.cpp

    Mesh/2D/Structured/CartesianMesh2D.cpp

    Discretization/LinearCombination/LinearCombination.hpp
)

set(HEADERS
    Utils/Types.h

    Boundary/BoundaryCondition.h

    Mesh/2D/Structured/CartesianMesh2D.h

    Discretization/LinearCombination/LinearCombination.h
)

set(CMAKE_CXX_STANDARD 20)

set(TARGET_NAME fvm_test)
add_executable(${TARGET_NAME} ${SOURCES} ${HEADERS})


add_subdirectory(Eigen)

target_link_libraries(${TARGET_NAME} eigen)
target_include_directories(${TARGET_NAME} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
