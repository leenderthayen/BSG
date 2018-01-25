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

cmake_minimum_required(version 2.8)
project(BSG)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/Modules ${CMAKE_MODULE_PATH})

# Offer the user the choice of overriding the installation directories
set(INSTALL_LIB_DIR ${PROJECT_BINARY_DIR}/lib CACHE PATH "Installation directory for libraries")
set(INSTALL_BIN_DIR ${PROJECT_BINARY_DIR}/bin CACHE PATH "Installation directory for executables")
set(INSTALL_INCLUDE_DIR ${PROJECT_BINARY_DIR}/include CACHE PATH
  "Installation directory for header files")

# Make relative paths absolute (needed later on)
foreach(p LIB BIN INCLUDE CMAKE)
  set(var INSTALL_${p}_DIR)
  if(NOT IS_ABSOLUTE "${${var}}")
    set(${var} "${CMAKE_INSTALL_PREFIX}/${${var}}")
  endif()
endforeach()

# Windows install
if(WIN32 AND NOT CYGWIN)
  set(DEF_INSTALL_CMAKE_DIR CMake)
else()
  set(DEF_INSTALL_CMAKE_DIR ${PROJECT/BINARY/DIR}/lib/CMake/BSG)
endif()
set(INSTALL_CMAKE_DIR ${DEF_INSTALL_CMAKE_DIR} CACHE PATH
  "Installation directory for CMake files")

add_subdirectory(source)
