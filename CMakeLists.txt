

#--------------------------------------------
#
#   CMAKE CONFIG FILE FOR HEAT
#
#   Author: Roman Keller
#    Start: 22.03.2017
#     Last: 22.03.2017
#--------------------------------------------
cmake_minimum_required(VERSION 3.7)
project (HEAT)
enable_language (Fortran)

#---------------------
# Set up the compiler
#---------------------
set (CMAKE_Fortran_COMPILER gfortran)

# set the preporcessor directives for compilation
set (PREPROC_DIRECTIVES "-DVERBOSE")
# set Compiler Directives to use with every BUILD_TYPE
set (GLOBAL_FLAGS "-ffree-form -Wextra -Wconversion -Wimplicit-interface -Wsurprising -ffree-line-length-130")
# -std=f2008 -Wall  -Werror
# set the compiler flags
set (CMAKE_Fortran_FLAGS_RELEASE "${PREPROC_DIRECTIVES} ${GLOBAL_FLAGS}  -O3")
set (CMAKE_Fortran_FLAGS_DEBUG   "${PREPROC_DIRECTIVES} ${GLOBAL_FLAGS} -fbacktrace -ffpe-trap=invalid,zero,overflow -fcheck=all -fbounds-check -g")
#set (CMAKE_BUILD_TYPE Debug)
set (CMAKE_BUILD_TYPE Release)

#------------------------
# Define the source path
#------------------------
file (GLOB_RECURSE sources src/*.F90)

#--------------------------
# Include HDF5 Library
#--------------------------
# On DLR MAchine hdf5 is loaded by module 
#set (CMAKE_PREFIX_PATH /media/LDS/module/software/libs/hdf5/1.8.17_gnu)
set(HDF5_ROOT /media/LDS/module/software/libs/hdf5/1.8.17_gnu)
find_package(HDF5 REQUIRED)

include_directories(${HDF5_INCLUDE_DIRS})
# on DLR Machine pdf5.mod is in shared subdir
include_directories(${HDF5_INCLUDE_DIRS}/shared)

#-------------------------
# Define the build & link
#-------------------------
add_executable (heat.exe ${sources})
#target_link_libraries(heat.exe "${HDF5}" "${HDF5_Fortran_LIBRARIES}")
target_link_libraries(heat.exe "hdf5" "hdf5_fortran")

add_subdirectory(tools)
add_subdirectory(test)
