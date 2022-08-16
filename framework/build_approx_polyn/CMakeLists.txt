# author: Ivan Novikau
cmake_minimum_required(VERSION 3.17)

# name of the executable file:
set(OUTPUT_EXE "approx_polyn")

# DEBUG or RELEASE:
set(CMAKE_BUILD_TYPE Debug)

# Set the project:
project(APP LANGUAGES CUDA CXX VERSION 0.0)

# if(NOT DEFINED CMAKE_CUDA_STANDARD)
#     set(CMAKE_CUDA_STANDARD 11)
#     set(CMAKE_CUDA_STANDARD_REQUIRED ON)
# endif()

SET(CMAKE_CUDA_COMPILER /usr/local/cuda-11.6/bin/nvcc)
SET(CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES /usr/local/cuda-11.6/targets/x86_64-linux/include)


# output some information:
message(STATUS "CMake version: ${CMAKE_VERSION}")
message(STATUS "CMake build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "CUDA Toolkit include: ${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES}")
message(STATUS "CUDA COMPILER ID: ${CMAKE_CUDA_COMPILER_ID}")

# set directories, where the executable file will be stored: 
SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE "./")
SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG   "./")

# find CUDA Toolkit:
find_package(CUDAToolkit)

# # find MPI
# find_package(MPI REQUIRED  COMPONENTS MPICXX)

# find HDF5
find_package(HDF5 REQUIRED  COMPONENTS CXX HL)

# create the target:
add_executable(${OUTPUT_EXE} ../src/ApprPolynFourier.cu)

# set C++ compiler version:
set_property(TARGET ${OUTPUT_EXE} PROPERTY CXX_STANDARD 17)

# set CUDA architecture:
set_property(TARGET ${OUTPUT_EXE} PROPERTY CUDA_ARCHITECTURES 75)

# link CUDA libraries:
target_link_libraries(${OUTPUT_EXE} PRIVATE CUDA::cudart_static)

# #link MPI library
# target_link_libraries(${OUTPUT_EXE} PRIVATE MPI::MPI_CXX)

#link HDF5 library
target_link_libraries(${OUTPUT_EXE} PRIVATE HDF5::HDF5)

# for correct debuging on device:
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    # target_compile_options(${OUTPUT_EXE} PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:-G>)
    target_compile_options(${OUTPUT_EXE} PRIVATE "$<$<AND:$<CONFIG:Debug>,$<COMPILE_LANGUAGE:CUDA>>:-G;-src-in-ptx>")
endif()