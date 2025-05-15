add_library(${LIBRARY_NAME})

add_subdirectory(Boundary)
add_subdirectory(Config)
add_subdirectory(Discretization)
add_subdirectory(Mesh)
add_subdirectory(Solvers)
add_subdirectory(Utils)
add_subdirectory(Visualization)
if (${BUILD_MESH_DRAWER})
    add_subdirectory(MeshDrawer)
endif()

target_include_directories(${LIBRARY_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# Sanitizers
if (UNIX)
    if (${TSAN})
        set(SANITIZER_FLAGS -fsanitize=thread)
    else()
        set(SANITIZER_FLAGS -fsanitize=undefined -fsanitize=address)
endif()
endif()

if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    message(STATUS "Sanitize flags: " ${SANITIZER_FLAGS})
endif()

target_compile_options(${LIBRARY_NAME} PUBLIC    
    $<$<CONFIG:Debug>: -Wall ${SANITIZER_FLAGS}>
    $<$<CONFIG:Release>: -march=native -flto>
)

target_compile_options(${LIBRARY_NAME} PRIVATE    
    $<$<CONFIG:Debug>: -Wall -Wextra -Wpedantic -Werror>
)

target_link_options(${LIBRARY_NAME} PUBLIC
    $<$<CONFIG:Debug>: ${SANITIZER_FLAGS}>
    $<$<CONFIG:Release>: -flto>
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


include(FetchContent)

FetchContent_Declare(
  mimalloc
  GIT_REPOSITORY https://github.com/microsoft/mimalloc.git
  GIT_TAG        v2.1.2
)
FetchContent_MakeAvailable(mimalloc)

target_link_libraries(${LIBRARY_NAME} PUBLIC 
    eigen 
    mimalloc
)
