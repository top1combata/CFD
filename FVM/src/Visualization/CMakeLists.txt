cmake_minimum_required(VERSION 3.20)

include(FetchContent)

FetchContent_Declare(
    sfml
    GIT_REPOSITORY https://github.com/SFML/SFML
    GIT_TAG 3.0.1
)
set(SFML_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
set(SFML_BUILD_DOCS OFF CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(sfml)

FetchContent_Declare(
    imgui
    GIT_REPOSITORY https://github.com/ocornut/imgui
    GIT_TAG af987eb
)
FetchContent_MakeAvailable(imgui)
set (IMGUI_DIR ${imgui_SOURCE_DIR})

FetchContent_Declare(
    imgui-sfml
    GIT_REPOSITORY https://github.com/SFML/imgui-sfml
    GIT_TAG 070d6ab
)
set(IMGUI_SFML_FIND_SFML OFF CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(imgui-sfml)


target_link_libraries(${LIBRARY_NAME} PRIVATE
    SFML::Graphics
    SFML::Window
    SFML::System
    ImGui-SFML::ImGui-SFML
)

target_sources(${LIBRARY_NAME} PRIVATE
    MeshDrawer.cpp
    PostProcessor.cpp
)

add_subdirectory(Window)
add_subdirectory(Drawing)
