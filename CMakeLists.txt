cmake_minimum_required(VERSION 3.20)
project(HeatEquation2D_MPI LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(MPI REQUIRED)

add_subdirectory(extern/json)
add_subdirectory(extern/Catch2)

add_executable(heat_solver 
    src/main.cpp 
    src/solver/solve.cpp
    src/solver/error.cpp
    src/solver/temperature.cpp
    src/solver/apply_boundaries.cpp
    src/param.cpp
)

target_include_directories(heat_solver PRIVATE include)

target_link_libraries(heat_solver 
    PRIVATE 
    MPI::MPI_CXX
    nlohmann_json::nlohmann_json
    Catch2::Catch2WithMain
)


# Testing
enable_testing()

add_executable(tests
    tests/build_test.cpp
)

target_link_libraries(tests PRIVATE Catch2::Catch2WithMain)

include(CTest)
include(Catch)
catch_discover_tests(tests)