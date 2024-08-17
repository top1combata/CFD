include(FetchContent)
FetchContent_Declare(
  googletest
  GIT_REPOSITORY "https://github.com/google/googletest"
  GIT_TAG "v1.15.2"
)

set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

enable_testing()

set(MESH_TEST_SOURCES
  Mesh/CartesianMesh2D.test.cpp
  Mesh/PolyMesh2D.test.cpp
)

set(DISCRETIZATION_TEST_SOURCES
  Discretization/LinearCombination.test.cpp
)

set(TEST_EXEC_NAME FVM-test)

add_executable(${TEST_EXEC_NAME} 
  ${MESH_TEST_SOURCES}
  ${DISCRETIZATION_TEST_SOURCES}
)

target_include_directories(${TEST_EXEC_NAME} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(${TEST_EXEC_NAME} GTest::gtest_main)
target_link_libraries(${TEST_EXEC_NAME} ${LIBRARY_NAME})

include(GoogleTest)
gtest_discover_tests(${TEST_EXEC_NAME})