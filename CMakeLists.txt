cmake_minimum_required(VERSION 3.12)
project(regrider)

# Add the cmake folder so the FindSphinx module is found
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

# output paths
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# write compile_commands.json file
set(CMAKE_EXPORT_COMPILE_COMMANDS TRUE)

# compile flags
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native")
set(CMAKE_CXX_FLAGS_DEBUG "-O1 -g -fsanitize=address,undefined,leak -fno-omit-frame-pointer -DDEBUG")

set(SRC
    src/utils.cpp
    src/grid.cpp
    src/gbptrees.cpp
    )

add_library(regrider_lib STATIC ${SRC})
set_property(TARGET regrider_lib PROPERTY CXX_STANDARD 20)

target_include_directories(regrider_lib PUBLIC "${CMAKE_SOURCE_DIR}/include")

find_package(OpenMP REQUIRED)
target_link_libraries(regrider_lib PRIVATE OpenMP::OpenMP_CXX)

add_subdirectory(vendor/fmt EXCLUDE_FROM_ALL)
target_link_libraries(regrider_lib PUBLIC fmt::fmt)

target_link_libraries(regrider_lib PUBLIC libfftw3f_omp.so libfftw3f.so)
find_path(FFTW_INCLUDE_DIR NAMES fftw3f.h fftw3.h)
target_include_directories(regrider_lib PUBLIC "${FFTW_INCLUDE_DIR}")

add_executable(regrider src/main.cpp)
target_link_libraries(regrider PRIVATE regrider_lib)

add_subdirectory("docs")

enable_testing()
add_subdirectory(tests)
