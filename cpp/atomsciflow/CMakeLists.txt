# Usage:
#   find_package(atomsciflow)
#   include_directories(${atomsciflow_INCLUDE_DIR})
#   link_libraries(${atomsciflow_LIBRARIES})
# if atomsciflow install dir is not in CMAKE MODULE PATH
# cmake -DCMAKE_PREFIX_PATH=/PATH/TO/atomsciflow/CMAKE_INSTALL_PREFIX
# or add it to the CMAKE_MODULE_PATH environment variable


cmake_minimum_required(VERSION 3.5)

project(atomsciflow LANGUAGES CXX VERSION 0.0.1)

#set(CMAKE_INCLUDE_CURRENT_DIR ON)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC ")

find_package(Boost COMPONENTS program_options filesystem system REQUIRED)
include_directories( ${Boost_INCLUDE_DIR}  )

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

#include_directories("./include")
include_directories("./src")

#AUX_SOURCE_DIRECTORY(./src PROJECT_LIB_SOURCES)
#file(GLOB_RECURSE PROJECT_LIB_SOURCES 
#    "./src/*" 
#    "./include/*"
#    )
file(GLOB PROJECT_LIB_SOURCES
    "./src/atomsciflow/*.cpp"
    "./src/atomsciflow/base/*.cpp"
    "./src/atomsciflow/parser/*.cpp"
    "./src/atomsciflow/abinit/*.cpp"
    "./src/atomsciflow/cp2k/*.cpp"
    "./src/atomsciflow/qe/*.cpp"
    "./src/atomsciflow/cube_handle/*.cpp"
)

file(GLOB PROJECT_HEADER_FILES
    "./src/atomsciflow/*.h"
    "./src/atomsciflow/base/*.h"
    "./src/atomsciflow/parser/*.h"
    "./src/atomsciflow/abinit/*.h"
    "./src/atomsciflow/cp2k/*.h"
    "./src/atomsciflow/qe/*.h"
    "./src/atomsciflow/cube_handle/*.h"
)

# set output path 
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib) # static lib
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib) # 
set(LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib) # 
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)


add_library(atomsciflow STATIC
    ${PROJECT_LIB_SOURCES}
)


add_executable(asflow.x 
    "./cmd/asflow.cpp"
)
target_link_libraries(asflow.x
    atomsciflow
    ${Boost_LIBRARIES}
)


add_executable(asflow-traj.x
    "./cmd/asflow_traj.cpp"
)
target_link_libraries(asflow-traj.x
    atomsciflow
    ${Boost_LIBRARIES}
)

add_executable(asflow-cube.x
    "./cmd/asflow_cube.cpp"
)
target_link_libraries(asflow-cube.x
    atomsciflow
    ${Boost_LIBRARIES}
)

add_executable(asflow-calc.x
    "./cmd/asflow_calc.cpp"
)
target_link_libraries(asflow-calc.x
    atomsciflow
    ${Boost_LIBRARIES}
)

# set install
install(
    DIRECTORY ${CMAKE_SOURCE_DIR}/src/
    DESTINATION include #include/atomsciflow
    FILES_MATCHING
    PATTERN "*.h"
)
#install(DIRECTORY ${PROJECT_BINARY_DIR}/lib/ DESTINATION lib)
#install(TARGETS atomsciflow EXPORT atomsciflow DESTINATION lib)
#install(DIRECTORY ${PROJECT_BINARY_DIR}/bin/ DESTINATION bin)
#install(TARGETS asflow.x asflow-traj.x asflow-cube.x asflow-calc.x  DESTINATION bin)

install(
    TARGETS atomsciflow asflow.x asflow-traj.x asflow-cube.x asflow-calc.x
    EXPORT atomsciflow
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    RUNTIME DESTINATION bin
)

install(
    EXPORT atomsciflow
    FILE atomsciflow.cmake
    DESTINATION lib/cmake/atomsciflow
)

configure_file(./cmake/atomsciflow-config.cmake.in ${PROJECT_BINARY_DIR}/atomsciflow-config.cmake @ONLY)
install(
    FILES ${CMAKE_BINARY_DIR}/atomsciflow-config.cmake
    DESTINATION lib/cmake/atomsciflow
)
