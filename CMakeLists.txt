cmake_minimum_required (VERSION 2.6)

#--------------------------------------------
#
#   CMAKE CONFIG FILE FOR    NS3DNEO
#
#--------------------------------------------
project (HEAT)
enable_language (Fortran)


#---------------------
# Set up the compiler
#---------------------
set (CMAKE_Fortran_COMPILER gfortran)

# set the preporcessor directives for compilation
set (PREPROC_DIRECTIVES "-DVERBOSE")
set (GLOBAL_FLAGS "-ffree-form -Wextra -Wconversion -Wimplicit-interface -Wsurprising -ffree-line-length-130")
# -std=f2008 -Wall 
# set the compiler flags
set (CMAKE_Fortran_FLAGS_RELEASE "${PREPROC_DIRECTIVES} ${GLOBAL_FLAGS}  -O3")
set (CMAKE_Fortran_FLAGS_DEBUG   "${PREPROC_DIRECTIVES} ${GLOBAL_FLAGS} -fbacktrace -ffpe-trap=invalid,zero,overflow -fcheck=all -fbounds-check -g")
set (CMAKE_BUILD_TYPE Debug)


#------------------------
# Define the source path
#------------------------
file (GLOB_RECURSE sources src/*.F90)


#--------------------------
# Set up the EAS3 paths :(
#--------------------------
#set (EAS3_LIBEXT "_x86_64_r8i4")
find_package(HDF5 REQUIRED PATHS /media/LDS/module/software/libs/hdf5/1.8.17_gnu/share/cmake/ )
#find_library(EAS3_BASELIB   NAMES "baselib${EAS3_LIBEXT}"   HINTS "$ENV{BASELIB_DIR}/${EAS3_LIBEXT}")
#find_library(EAS3_N3DLIB    NAMES "n3dlib${EAS3_LIBEXT}"    HINTS "$ENV{N3DLIB_DIR}/${EAS3_LIBEXT}")
#find_library(EAS3_EAS3IOMOD NAMES "eas3iomod${EAS3_LIBEXT}" HINTS "$ENV{EAS3IOMOD_DIR}/${EAS3_LIBEXT}")

#include_directories("$ENV{BASELIB_DIR}/${EAS3_LIBEXT}")
#include_directories("$ENV{N3DLIB_DIR}/${EAS3_LIBEXT}")
#include_directories("$ENV{EAS3IOMOD_DIR}/${EAS3_LIBEXT}")
include_directories(${HDF5_INCLUDE_DIRS})

#-------------------------
# Define the build & link
#-------------------------
add_executable (heat.exe ${sources})
target_link_libraries(heat.exe  "${HDF5_LIBRARIES}")
