cmake_minimum_required(VERSION 3.5)

project(atomscikit-pybind11 LANGUAGES CXX)

#set(CMAKE_INCLUDE_CURRENT_DIR ON)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


include_directories("../atomscikit/include")

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


# -----------------------------------------------------
# for pybind11 and scikit-build
# -----------------------------------------------------
file(GLOB CPPTEST_SOURCES
    "./cpptest/*.cpp"
    "../atomscikit/src/cpptest/*.cpp"
)

pybind11_add_module(cpptest MODULE ${CPPTEST_SOURCES})

# scikit-build
if(SKBUILD)
    install(TARGETS cpptest DESTINATION .)
endif()