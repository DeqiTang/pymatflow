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

add_library(atomscikitf ${LIBATOMSCIKITF_SRC})
add_library(atomscikitf-c-binding ${LIBATOMSCIKITF_SRC} ${LIBATOMSCIKITF_C_BINDING_SRC})

add_executable(askitf-cube-cube-to-vtk.x
    "./cmd/cube_to_vtk_program.f90"
    ${LIBATOMSCIKITF_SRC}
)
add_executable(askitf-cube-diff-1d.x
    "./cmd/cube_diff_1d_program.f90"
    ${LIBATOMSCIKITF_SRC}
)
add_executable(askitf-cube-1d.x
    "./cmd/cube_1d_program.f90"
    ${LIBATOMSCIKITF_SRC}
)



