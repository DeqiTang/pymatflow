cmake_minimum_required(VERSION 3.5)

project(atomscikit LANGUAGES CXX)

#set(CMAKE_INCLUDE_CURRENT_DIR ON)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

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

include_directories("./include")

#AUX_SOURCE_DIRECTORY(./src PROJECT_LIB_SOURCES)
#file(GLOB_RECURSE PROJECT_LIB_SOURCES 
#    "./src/*" 
#    "./include/*"
#    )
file(GLOB PROJECT_LIB_SOURCES
    "./src/*.cpp"
    "./src/base/*.cpp"
    "./src/parser/*.cpp"
    "./src/abinit/*.cpp"
    "./src/cp2k/*.cpp"
    "./src/qe/*.cpp"
    "./include/askit/*.h"
    "./include/askit/base/*.h"
    "./include/askit/parser/*.h"
    "./include/askit/abinit/*.h"
    "./include/askit/cp2k/*.h"
    "./include/askit/qe/*.h"
    )

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)

add_library(atomscikit STATIC
    ${PROJECT_LIB_SOURCES}
    )


add_executable(skit.x 
    ${PROJECT_LIB_SOURCES}
    "./cmd/skit.cpp"
    )
target_link_libraries(skit.x
    ${Boost_LIBRARIES}
    )


add_executable(skit-traj.x
    "./cmd/skit_traj.cpp"
    ${PROJECT_LIB_SOURCES}
    )
target_link_libraries(skit-traj.x
    ${Boost_LIBRARIES}
    )

add_executable(skit-cube.x
    "./cmd/skit_cube.cpp"
    ${PROJECT_LIB_SOURCES}
    )
target_link_libraries(skit-cube.x
    ${Boost_LIBRARIES}
    )

add_executable(askit-calc.x
    "./cmd/askit_calc.cpp"
    ${PROJECT_LIB_SOURCES}
    )
target_link_libraries(askit-calc.x
    ${Boost_LIBRARIES}
    )
