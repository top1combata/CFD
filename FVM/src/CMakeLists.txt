
set(SOURCES
    Utils/MatrixSolver.cpp
    Utils/Timer.cpp

    Config/Config.cpp

    Mesh/MeshBase.cpp
    Mesh/Geometry.cpp
    Mesh/2D/Structured/CartesianMesh2D.cpp
    Mesh/2D/Unstructured/PolyMesh2D.cpp
    Mesh/2D/Unstructured/Parse.cpp

    Discretization/LinearCombination/LinearCombination.cpp
    Discretization/Interpolation/Interpolation.hpp

    Solvers/SIMPLE/SimpleAlgorithm.cpp
)

set(HEADERS
    Utils/Types.h
    Utils/TypesOperations.h
    Utils/MatrixSolver.h
    Utils/Timer.h

    Config/Config.h

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

add_library(${LIBRARY_NAME} ${SOURCES} ${HEADERS})

target_include_directories(${LIBRARY_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# Sanitizers
if (UNIX)
    if (${TSAN})
        set(SANITIZER_FLAGS -fsanitize=thread)
    else()
        set(SANITIZER_FLAGS -fsanitize=address)
    endif()
endif()

if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    message(STATUS "Sanitize flags: " ${SANITIZER_FLAGS})
endif()

target_compile_options(${LIBRARY_NAME} PUBLIC    
    $<$<CONFIG:Debug>: -Wall -Wextra -Wpedantic -Werror ${SANITIZER_FLAGS}>
    $<$<CONFIG:Release>: -march=native>
)

target_link_options(${LIBRARY_NAME} PUBLIC
    $<$<CONFIG:Debug>: ${SANITIZER_FLAGS}>
)

# OpenMP
find_package(OpenMP)
if (OpenMP_CXX_FOUND)
    target_link_libraries(${LIBRARY_NAME} PUBLIC OpenMP::OpenMP_CXX)
else()
    message(WARNING "CMake can't find OpenMP\n")
endif()

#Eigen
add_subdirectory(Eigen)
target_link_libraries(${LIBRARY_NAME} PUBLIC eigen)
