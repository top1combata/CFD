project(cylinder-flow)

add_executable(${PROJECT_NAME})

target_sources(${PROJECT_NAME} PRIVATE
    main.cpp
)

target_link_libraries(${PROJECT_NAME} ${LIBRARY_NAME})

include(${CMAKE_CURRENT_LIST_DIR}/../scripts/ExecutePython.cmake)

message("Generating Cylinder example mesh...")

execute_python_script(
    SCRIPT_PATH ${CMAKE_CURRENT_LIST_DIR}/mesh_gen.py
    IMPORT_DIRS ${CMAKE_CURRENT_LIST_DIR}/../scripts
)

message("Cylinder example mesh generated\n")
