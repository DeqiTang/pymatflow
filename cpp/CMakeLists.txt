cmake_minimum_required(VERSION 3.5)

project(cpptest LANGUAGES CXX)

#set(CMAKE_INCLUDE_CURRENT_DIR ON)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


include_directories("./include")

# scikit-build
if(SKBUILD)
  message(STATUS "The project is built using scikit-build")
endif()


# pybind11
find_package(pybind11 REQUIRED)
include_directories(${PYTHON_INCLUDE_DIRS})

file(GLOB PROJECT_LIB_SOURCES
    "./src/cpptest/*.cpp"
    )

#set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
#SET (LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)

#add_library(cpptest STATIC
#    ${PROJECT_LIB_SOURCES}
#    )
#add_library(cpptest_shared SHARED
#    ${PROJECT_LIB_SOURCES}
#    )

#SET_TARGET_PROPERTIES (cpptest_shared PROPERTIES OUTPUT_NAME "cpptest")
#GET_TARGET_PROPERTY (OUTPUT_VALUE cpptest_shared OUTPUT_NAME)
#SET_TARGET_PROPERTIES (cpptest_shared PROPERTIES CLEAN_DIRECT_OUTPUT 1)
#SET_TARGET_PROPERTIES (cpptest PROPERTIES CLEAN_DIRECT_OUTPUT 1)

#pybind11_add_module(pyaskit SHARED ${PROJECT_LIB_SOURCES})
pybind11_add_module(pyaskit MODULE ${PROJECT_LIB_SOURCES})


# scikit-build
if(SKBUILD)
    install(TARGETS pyaskit DESTINATION .)
endif()