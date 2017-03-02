cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
project(BDSM)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_MATH} "${CMAKE_SOURCE_DIR}")

FIND_PACKAGE(GSL REQUIRED)
include_directories(${GSL_INCLUDE_DIRS})

FIND_PACKAGE(ROOT REQUIRED)
include_directories(${ROOT_INCLUDE_DIRS})

FIND_PACKAGE(Boost COMPONENTS program_options REQUIRED)
INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIR})

include_directories(${PROJECT_SOURCE_DIR}/include)

file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cpp)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.h)
file(COPY data DESTINATION .)
file(COPY init DESTINATION .)

configure_file(${PROJECT_SOURCE_DIR}/config.txt ${PROJECT_BINARY_DIR}/config.txt COPYONLY)

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)

add_library(bdsm SHARED ${sources})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -L/usr/lib -lgsl -lgslcblas -lm -std=c++1y")

add_executable(BDSM BDSM.cc ${sources} ${headers})

target_link_libraries(BDSM ${GSL_LIBRARIES} ${ROOT_LIBRARIES} ${Boost_LIBRARIES})

install(TARGETS BDSM DESTINATION bin)
