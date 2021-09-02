cmake_minimum_required(VERSION 2.8.5)
project(askitf)
enable_language(Fortran)
set(CMAKE_BUILD_TYPE Release)

add_compile_options(-fPIC)

#set(BUILD_SHARED_LIBS ON)

set(CMAKE_Fortran_FLAGS "-g ${CMAKE_Fortran_FLAGS}")

file(GLOB_RECURSE LIBASKITF_SRC RELATIVE ${PROJECT_SOURCE_DIR} src/*.f90)
file(GLOB_RECURSE LIBASKITF_C_BINDING_SRC RELATIVE ${PROJECT_SOURCE_DIR} c_binding/*.f90)

add_library(askitf ${LIBASKITF_SRC})
add_library(askitf-c-binding ${LIBASKITF_SRC} ${LIBASKITF_C_BINDING_SRC})

add_executable(askit-cube-cube-to-vtk.x
    "./cmd/cube_to_vtk_program.f90"
    ${LIBASKITF_SRC}
)
add_executable(askit-cube-diff-1d.x
    "./cmd/cube_diff_1d_program.f90"
    ${LIBASKITF_SRC}
)
add_executable(askit-cube-1d.x
    "./cmd/cube_1d_program.f90"
    ${LIBASKITF_SRC}
)



