cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
project(BSG)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_MATH} "${CMAKE_SOURCE_DIR}")

FIND_PACKAGE(GSL REQUIRED)
include_directories(${GSL_INCLUDE_DIRS})

FIND_PACKAGE(ROOT REQUIRED)
include_directories(${ROOT_INCLUDE_DIRS})

include_directories(${PROJECT_SOURCE_DIR}/include)

file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cpp)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.h)
file(COPY data DESTINATION .)
file(COPY init DESTINATION .)

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)

add_library(bsg SHARED ${sources})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -L/usr/lib -lgsl -lgslcblas -lm -std=c++11")

add_executable(BSG BSG.cc ${sources} ${headers})

target_link_libraries(BSG ${GSL_LIBRARIES} ${ROOT_LIBRARIES})

install(TARGETS BSG DESTINATION bin)
