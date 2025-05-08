project(mesh-drawer)

add_executable(${PROJECT_NAME})

target_sources(${PROJECT_NAME} PRIVATE
    main.cpp
    MeshDrawer.cpp
)

set_target_properties(${PROJECT_NAME} PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}"
)

target_link_libraries(${PROJECT_NAME} ${LIBRARY_NAME})
