cmake_minimum_required(VERSION 3.10)
project(pulse_cpp)

file(GLOB SOURCES "src/*.cpp")

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



