cmake_minimum_required(VERSION 3.16)

project(area_and_topology_preserving_polygon_simplification LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_executable(area_and_topology_preserving_polygon_simplification
    src/main.cpp
    src/geometry.cpp
    src/io_csv.cpp
    src/simplify.cpp
    src/validation.cpp
)

target_include_directories(area_and_topology_preserving_polygon_simplification
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/include
)
