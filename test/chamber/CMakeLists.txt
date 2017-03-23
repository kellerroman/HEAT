add_executable(chamber_gridgen.exe gridgen.F90 ../../src/const.F90)
target_link_libraries(chamber_gridgen.exe module_gridgen)

#add_custom_target(run ./gabel3D_gridgen.exe)

