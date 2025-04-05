cmake_minimum_required(VERSION 3.21)

set(PROJECT_NAME FVM)
project(${PROJECT_NAME})

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(LIBRARY_NAME ${PROJECT_NAME})

add_subdirectory(src)
add_subdirectory(test)

add_executable(main main.cpp)

target_link_libraries(main ${LIBRARY_NAME})