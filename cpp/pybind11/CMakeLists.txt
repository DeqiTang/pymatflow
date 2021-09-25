cmake_minimum_required(VERSION 3.5)

project(atomsciflow-pybind11 LANGUAGES CXX)

#set(CMAKE_INCLUDE_CURRENT_DIR ON)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC ")

find_package(Boost COMPONENTS program_options filesystem system REQUIRED)
include_directories( ${Boost_INCLUDE_DIR}  )

#include_directories("../atomsciflow/include")
include_directories("../atomsciflow/src")

# scikit-build
if(SKBUILD)
  message(STATUS "The project is built using scikit-build")
endif()


# pybind11
find_package(pybind11 REQUIRED)
include_directories(${PYTHON_INCLUDE_DIRS})

# OpenMP
find_package(OpenMP)
if(OPENMP_FOUND)
    message("OPENMP FOUND")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS"${CMAKE_EXE_LINKER_FLAGS}${OpenMP_EXE_LINKER_FLAGS}")
endif()

find_package(Armadillo REQUIRED)
include_directories(${ARMADILLO_INCLUDE_DIRS})
link_libraries(${ARMADILLO_LIBRARIES})

find_package(LAPACK REQUIRED)
find_package(BLAS REQUIRED)

# -----------------------------------------------------
# for pybind11 and scikit-build
# -----------------------------------------------------
file(GLOB CPPTEST_SOURCES
    "./cpptest/*.cpp"
    "../atomsciflow/src/cpptest/*.cpp"
)
pybind11_add_module(cpptest MODULE ${CPPTEST_SOURCES})

file(GLOB CUBE_HANDLE_SOURCES
    "./cube_handle/*.cpp"
    "../atomsciflow/src/atomsciflow/cube_handle/*.cpp"
    "../atomsciflow/src/atomsciflow/parser/*.cpp"
    "../atomsciflow/src/atomsciflow/base/*.cpp"
)
pybind11_add_module(cube_handle MODULE ${CUBE_HANDLE_SOURCES})

# scikit-build
if(SKBUILD)
    install(TARGETS cpptest DESTINATION .)
    install(TARGETS cube_handle DESTINATION .)
endif()
