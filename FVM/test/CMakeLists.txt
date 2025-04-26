include(FetchContent)
FetchContent_Declare(
  googletest
  GIT_REPOSITORY "https://github.com/google/googletest"
  GIT_TAG "v1.15.2"
)

set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

enable_testing()

set(TEST_EXECUTABLE FVM-test)

add_executable(${TEST_EXECUTABLE})

add_subdirectory(Discretization)
add_subdirectory(Mesh)
add_subdirectory(Solvers)

target_include_directories(${TEST_EXECUTABLE} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(${TEST_EXECUTABLE} GTest::gtest_main)
target_link_libraries(${TEST_EXECUTABLE} ${LIBRARY_NAME})

include(GoogleTest)
gtest_discover_tests(${TEST_EXECUTABLE})
