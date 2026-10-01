# TIDE -- C++ implementation and benchmarks used for the accompanying manuscript.
#   Linux/macOS : cmake -S . -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build -j
#   Windows MSVC: cmake -S . -B build -A x64 && cmake --build build --config Release
#   Windows MinGW: cmake -S . -B build -G "MinGW Makefiles" -DCMAKE_BUILD_TYPE=Release && cmake --build build
# No third-party dependencies. OpenMP is optional (timings use OMP_NUM_THREADS=1).
# Fast-math is deliberately NOT enabled.
cmake_minimum_required(VERSION 3.15)
project(tide LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
endif()
if(MSVC)
  add_compile_options(/O2 /fp:precise /W4)
else()
  add_compile_options(-O3 -fno-fast-math -Wall -Wextra)
endif()
find_package(OpenMP)
add_library(tide_headers INTERFACE)
target_include_directories(tide_headers INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)
if(OpenMP_CXX_FOUND)
  target_link_libraries(tide_headers INTERFACE OpenMP::OpenMP_CXX)
  target_compile_definitions(tide_headers INTERFACE TIDE_OPENMP)
endif()
add_executable(tide_tools src/tide/tide_tools.cpp)
add_executable(bench_scaling src/benchmarks/bench_scaling.cpp)
target_link_libraries(tide_tools PRIVATE tide_headers)
target_link_libraries(bench_scaling PRIVATE tide_headers)
# Print the compiler identity so it can be quoted in the paper
file(WRITE ${CMAKE_BINARY_DIR}/compiler.txt "${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}\n")
