add_executable(gabel3D_gridgen.exe gridgen.F90 ../../src/const.F90)
#target_link_libraries(gridgen.exe "hdf5" "hdf5_fortran")
target_link_libraries(gabel3D_gridgen.exe module_gridgen)
#target_link_libraries(gridgen.exe PRIVATE "hdf5_fortran")
#add_custom_target(run ./gabel3D_gridgen.exe)
