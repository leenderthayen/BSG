#-----------------------------------------------------------------------
# - Enforce an out-of-source builds before anything else
#
if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
    message(STATUS "BSG requires an out-of-source build.")
    message(STATUS "Please remove these files from ${CMAKE_BINARY_DIR} first:")
    message(STATUS "CMakeFiles")
    message(STATUS "Once these files are removed, create a separate directory")
    message(STATUS "and run CMake from there")
    message(FATAL_ERROR "in-source build detected")
endif()

cmake_minimum_required(VERSION 3.1)
project(BSG VERSION 1.0 LANGUAGES CXX)

# - Prepend our own CMake Modules to the search path
# NB: if our custom modules include others that we don't supply, those in
# the base path will be used, so watch for incompatibilities!!
#
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/Modules ${CMAKE_MODULE_PATH})
include(CMakeDependentOption)
include(GNUInstallDirs)
include(CMakeUninstallTarget)

#---------------------------------------------------------------------------------------
# compiler config
#---------------------------------------------------------------------------------------
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU" OR "${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
    add_compile_options("-Wall")
    add_compile_options("-Wextra")
    add_compile_options("-Wconversion")
    add_compile_options("-pedantic")
    add_compile_options("-Wfatal-errors")
endif()

#---------------------------------------------------------------------------------------
# Options
#---------------------------------------------------------------------------------------
option(BSG_BUILD_TESTS "Build tests" ON)
option(BSG_INSTALL "Generate the install target." OFF)
option(BSG_BUILD_EXAMPLES "Build examples" ON)

find_package(spdlog REQUIRED)
find_package(GSL REQUIRED)

#---------------------------------------------------------------------------------------
# BSG target
#---------------------------------------------------------------------------------------
add_library(BSG SHARED src/Generator.cc src/ConfigParser.cc src/SpectralFunctions.cc)

set(NHL_INCLUDE_DIRS "/Users/leenderthayen/git/NHL/include")
set(PDS_INCLUDE_DIRS "/Users/leenderthayen/git/PDS/include")
#set(NME_INCLUDE_DIRS "/home/leendert/git/NME/include")
#set(NME_LIBRARIES "-L/home/leendert/git/NME/build -lNME")

target_include_directories(
    BSG
    PUBLIC
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
    "$<INSTALL_INTERFACE:include>"
)
target_include_directories(BSG PUBLIC ${spdlog_INCLUDE_DIRS})
target_include_directories(BSG PUBLIC ${NHL_INCLUDE_DIRS})
#target_include_directories(BSG PUBLIC ${NME_INCLUDE_DIRS})
target_include_directories(BSG PUBLIC ${PDS_INCLUDE_DIRS})
#target_link_libraries(BSG PUBLIC ${NME_LIBRARIES})
target_include_directories(BSG PUBLIC ${GSL_INCLUDE_DIRS})

set(HEADER_BASE "${CMAKE_CURRENT_SOURCE_DIR}/include")

target_link_libraries(BSG ${GSL_LIBRARIES})

if(BSG_BUILD_TESTS)
    include(CTest)
    add_subdirectory(tests)
endif()

if(BSG_BUILD_EXAMPLES)
    add_subdirectory(examples)
endif()
