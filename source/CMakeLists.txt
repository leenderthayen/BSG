set(BSG_VERSION "1.0")
set(BSG_VERSION_MAJOR "1")
set(BSG_VERSION_MINOR "0")
set(BSG_LAST_UPDATE "November 23, 2018")

set(NME_VERSION "1.0")
set(NME_VERSION_MAJOR "1")
set(NME_VERSION_MINOR "0")
set(NME_LAST_UPDATE "November 23, 2018")

# Look for required packages
find_package(GSL REQUIRED)
include_directories(${GSL_INCLUDE_DIRS})
find_package(ROOT REQUIRED)
include_directories(${ROOT_INCLUDE_DIRS})
find_package(Boost COMPONENTS program_options REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})
find_package(Threads)

find_package(spdlog REQUIRED)
include_directories(${spdlog_INCLUDE_DIRS})

# Configure the version header files
configure_file("${PROJECT_SOURCE_DIR}/cmake/Templates/BSGConfig.h.in" "${PROJECT_BINARY_DIR}/BSGConfig.h")
configure_file("${PROJECT_SOURCE_DIR}/cmake/Templates/NMEConfig.h.in" "${PROJECT_BINARY_DIR}/NMEConfig.h")
include_directories("${PROJECT_BINARY_DIR}")

include_directories("${PROJECT_SOURCE_DIR}/source/nme/include")
include_directories("${PROJECT_SOURCE_DIR}/source/bsg/include")

# Copy the data directory
#configure_file("${PROJECT_SOURCE_DIR}/data/config.txt" "${PROJECT_BINARY_DIR}/bin/config.txt")
#configure_file("${PROJECT_SOURCE_DIR}/data/ExchangeData.dat" "${PROJECT_BINARY_DIR}/bin/ExchangeData.dat")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -std=c++11")

add_subdirectory(nme)
add_subdirectory(bsg)
add_subdirectory(executables)
