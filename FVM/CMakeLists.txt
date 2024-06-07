cmake_minimum_required(VERSION 3.21)

set(PROJECT_NAME FVM)

set(SOURCES
    main.cpp
)

set(HEADERS
    Utils/Types.h
)

set(TARGET_NAME fvm_test)
add_executable(${TARGET_NAME} ${SOURCES} ${HEADERS})


add_subdirectory(Eigen)

target_link_libraries(${TARGET_NAME} eigen)
target_include_directories(${TARGET_NAME} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
