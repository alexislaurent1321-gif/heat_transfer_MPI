cmake_minimum_required(VERSION 3.20)
project(HeatEquation2D_MPI LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(MPI REQUIRED)

add_subdirectory(extern/fmt)
add_subdirectory(extern/json)

add_executable(heat_solver 
    src/main.cpp 
    src/solve.cpp
)

target_include_directories(heat_solver PRIVATE include)

target_link_libraries(heat_solver 
    PRIVATE 
    MPI::MPI_CXX
    fmt::fmt-header-only
    nlohmann_json::nlohmann_json
)