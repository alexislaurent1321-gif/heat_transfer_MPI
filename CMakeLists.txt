cmake_minimum_required(VERSION 3.20)
project(HeatEquation2D_MPI LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(MPI REQUIRED)

include(FetchContent)
FetchContent_Declare(
  fmt
  GIT_REPOSITORY https://github.com/fmtlib/fmt.git
  GIT_TAG 10.1.1
)
FetchContent_MakeAvailable(fmt)

include(FetchContent)
FetchContent_Declare(
  json
  URL https://github.com/nlohmann/json/releases/download/v3.11.2/json.tar.xz
)
FetchContent_MakeAvailable(json)

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