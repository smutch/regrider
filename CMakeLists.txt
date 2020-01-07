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
set(CMAKE_CXX_FLAGS_DEBUG "-O1 -g -fsanitize=address,undefined,leak -fno-omit-frame-pointer")

set(SRC
    src/main.cpp
    src/grid.cpp
    src/gbptrees.cpp
    )

add_executable(regrider ${SRC})
set_property(TARGET regrider PROPERTY CXX_STANDARD 20)

target_include_directories(regrider PRIVATE "${CMAKE_SOURCE_DIR}/include")

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(regrider PRIVATE OpenMP::OpenMP_CXX)
endif()

add_subdirectory(vendor/fmt EXCLUDE_FROM_ALL)
target_link_libraries(regrider PRIVATE fmt::fmt)

find_path(FFTW_LIBRARY_DIR libfftw3f_omp.a)
target_link_libraries(regrider PRIVATE "${FFTW_LIBRARY_DIR}/libfftw3f_omp.a" "${FFTW_LIBRARY_DIR}/libfftw3f.a")
find_path(FFTW_INCLUDE_DIR fftw3f.h)
target_include_directories(regrider PRIVATE "${FFTW_INCLUDE_DIR}")

add_subdirectory("docs")
