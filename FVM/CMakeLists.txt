cmake_minimum_required(VERSION 3.21)

project(FVM)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(LIBRARY_NAME ${PROJECT_NAME})

add_subdirectory(src)
add_subdirectory(test)

if (${BUILD_EXAMPLES})
    add_subdirectory(examples)
endif()
