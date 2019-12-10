cmake_minimum_required(VERSION 3.0)
project(regrider)

# output paths
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# write compile_commands.json file
set(CMAKE_EXPORT_COMPILE_COMMANDS TRUE)

# compile flags
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native")
set(CMAKE_CXX_FLAGS_DEBUG "-O1 -g -fsanitize=address,undefined,leak -fno-omit-frame-pointer")

set(SRC
    src/regrider.cpp
    )

add_executable(regrider ${SRC})
set_property(TARGET regrider PROPERTY CXX_STANDARD 20)

target_include_directories(regrider PRIVATE "${CMAKE_SOURCE_DIR}/include")

find_package(fmt)
target_link_libraries(regrider PRIVATE fmt::fmt)

find_path(FFTW_LIBRARY_DIR libfftw3f.a)
target_link_libraries(regrider PRIVATE "${FFTW_LIBRARY_DIR}/libfftw3f.a")
find_path(FFTW_INCLUDE_DIR fftw3f.h)
target_include_directories(regrider PRIVATE "${FFTW_INCLUDE_DIR}")
