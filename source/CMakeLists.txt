set(BSG_VERSION "1.1")
set(BSG_VERSION_MAJOR "1")
set(BSG_VERSION_MINOR "1")
set(BSG_LAST_UPDATE "April 30, 2019")

set(NME_VERSION "1.0")
set(NME_VERSION_MAJOR "1")
set(NME_VERSION_MINOR "0")
set(NME_LAST_UPDATE "April 30, 2018")

# Look for required packages
find_package(GSL REQUIRED)
include_directories(${GSL_INCLUDE_DIRS})
find_package(ROOT REQUIRED)
include_directories(${ROOT_INCLUDE_DIRS})
find_package(Boost COMPONENTS program_options REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})
find_package(Threads)

message(STATUS "BOOST: I: ${Boost_INCLUDE_DIRS} L: ${Boost_LIBRARIES}")

find_package(spdlog REQUIRED)
include_directories(${spdlog_INCLUDE_DIRS})

# Configure the version header files
configure_file("${PROJECT_SOURCE_DIR}/cmake/Templates/BSGConfig.h.in" "${PROJECT_BINARY_DIR}/BSGConfig.h")
configure_file("${PROJECT_SOURCE_DIR}/cmake/Templates/NMEConfig.h.in" "${PROJECT_BINARY_DIR}/NMEConfig.h")
include_directories("${PROJECT_BINARY_DIR}")

include_directories("${PROJECT_SOURCE_DIR}/source/nme/include")
include_directories("${PROJECT_SOURCE_DIR}/source/bsg/include")

set(EXTRA_COMPILE_FLAGS "-g -std=c++11")
set(EXTRA_LINKING_FLAGS "-lstdc++")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${EXTRA_COMPILE_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${EXTRA_LINKING_FLAGS}")

add_subdirectory(nme)
add_subdirectory(bsg)
add_subdirectory(executables)
