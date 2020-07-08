cmake_minimum_required(VERSION 3.10)
project(pulse_cpp)

file(GLOB SOURCES "src/*.cpp")

set(GNU_HOST x86_64-w64-mingw32)
set(CMAKE_SYSTEM_PROCESSOR "x86_64")

set(COMPILER_PREFIX "${GNU_HOST}-")

set(CMAKE_SYSTEM_NAME "Windows")
set(CMAKE_CROSSCOMPILING TRUE)
set(WIN32 TRUE)
set(MINGW TRUE)

set(CMAKE_C_COMPILER ${COMPILER_PREFIX}gcc)
set(CMAKE_CXX_COMPILER ${COMPILER_PREFIX}g++)
set(CMAKE_RC_COMPILER ${COMPILER_PREFIX}windres)

set(CMAKE_FIND_ROOT_PATH_MODE_PROGRAM NEVER)
set(CMAKE_FIND_ROOT_PATH_MODE_LIBRARY ONLY)
set(CMAKE_FIND_ROOT_PATH_MODE_INCLUDE ONLY)


#SET(CMAKE_FIND_ROOT_PATH  /usr/x86_64-w64-mingw32 /home/dcosta/mingw-install )

#set(CXX g++)
#set(ARMA_USE_CXX11)
set(ARMA_USE_LAPACK true)

set(CMAKE_CXX_STANDARD 11)
#set(CMAKE_CXX_FLAGS "-O3")
set(CMAKE_CXX_FLAGS "-g3 -Wall -pedantic")
#set(CMAKE_CXX_FLAGS "-g ")
#set(CMAKE_BUILD_TYPE RelWithDebInfo)

find_package(Armadillo REQUIRED)
find_package(OpenMP)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

include_directories(${AMRADILLO_INCLUDE_DIR})

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

add_executable(pulse_cpp
        ${SOURCES})

target_link_libraries(
        pulse_cpp
        ${ARMADILLO_LIBRARIES} )



