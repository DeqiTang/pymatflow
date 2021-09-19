cmake_minimum_required(VERSION 2.8.5)
project(atomscikitf)
enable_language(Fortran)
set(CMAKE_BUILD_TYPE Release)

# needed if the objects are later connected to cython to build for python import
add_compile_options(-fPIC) 

#set(BUILD_SHARED_LIBS ON)

set(CMAKE_Fortran_FLAGS "-g ${CMAKE_Fortran_FLAGS}")

file(GLOB_RECURSE LIBATOMSCIKITF_SRC RELATIVE ${PROJECT_SOURCE_DIR} src/*.f90)
file(GLOB_RECURSE LIBATOMSCIKITF_C_BINDING_SRC RELATIVE ${PROJECT_SOURCE_DIR} c_binding/*.f90)

# set output path 
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib) # static lib
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib) # 
set(LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib) # 
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# set install
install(DIRECTORY ${PROJECT_BINARY_DIR}/mod/ DESTINATION mod/atomscikitf)
install(DIRECTORY ${PROJECT_BINARY_DIR}/lib/ DESTINATION lib)
install(DIRECTORY ${PROJECT_BINARY_DIR}/bin/ DESTINATION bin)

add_library(atomscikitf ${LIBATOMSCIKITF_SRC})

add_library(atomscikitf-c-binding ${LIBATOMSCIKITF_C_BINDING_SRC})
target_link_libraries(atomscikitf-c-binding atomscikitf)

add_executable(askitf-cube-to-vtk.x
    "./cmd/cube_to_vtk_program.f90"
)
target_link_libraries(askitf-cube-to-vtk.x atomscikitf)

add_executable(askitf-cube-diff-1d.x
    "./cmd/cube_diff_1d_program.f90"
)
target_link_libraries(askitf-cube-diff-1d.x atomscikitf)

add_executable(askitf-cube-1d.x
    "./cmd/cube_1d_program.f90"
)
target_link_libraries(askitf-cube-1d.x atomscikitf)