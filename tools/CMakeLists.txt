add_executable(pre_config.exe pre_config.F90)
add_executable(create_xdmf.exe create_xdmf.F90)
target_link_libraries(create_xdmf.exe "hdf5_fortran")

add_library(module_gridgen module_gridgen.F90)
# this is required so the compiler finds the module when compiling a code using this module
target_include_directories (module_gridgen PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
target_link_libraries(module_gridgen INTERFACE "hdf5_fortran")

add_library(module_post module_post.F90)
