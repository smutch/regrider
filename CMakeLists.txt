cmake_minimum_required(VERSION 3.12)
project(regrider CXX)

# Add the cmake folder so the FindSphinx module is found
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

# output paths
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# write compile_commands.json file
set(CMAKE_EXPORT_COMPILE_COMMANDS TRUE)

set(default_build_type "Debug")
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
  set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE
      STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

# compile flags
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native")
set(CMAKE_CXX_FLAGS_DEBUG "-O1 -g -fsanitize=address,undefined,leak -fno-omit-frame-pointer -DDEBUG")

set(SRC
    src/utils.cpp
    src/grid.cpp
    src/gbptrees.cpp
    src/velociraptor.cpp
    )

add_library(regrider_lib STATIC ${SRC})
set_property(TARGET regrider_lib PROPERTY CXX_STANDARD 11)  # WARNING: This must be 11 for hdf5 to work!!!

target_include_directories(regrider_lib PUBLIC "${CMAKE_SOURCE_DIR}/include")

find_package(OpenMP REQUIRED)
target_link_libraries(regrider_lib PRIVATE OpenMP::OpenMP_CXX)

add_subdirectory(vendor/fmt EXCLUDE_FROM_ALL)
target_link_libraries(regrider_lib PUBLIC fmt::fmt)

find_package(FFTW REQUIRED)
target_include_directories(regrider_lib PUBLIC ${FFTW_INCLUDE_DIRS})
target_link_libraries(regrider_lib PUBLIC ${FFTW_LIBRARIES})

find_package(HDF5 REQUIRED COMPONENTS CXX HL)
target_link_libraries(regrider_lib PUBLIC HDF5::HDF5)

add_executable(regrider src/main.cpp)
target_link_libraries(regrider PRIVATE regrider_lib)

add_subdirectory("docs")

enable_testing()
add_subdirectory(tests)
