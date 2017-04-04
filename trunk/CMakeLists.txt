cmake_minimum_required(VERSION 2.8)
project(BSG)
 
set(BSG_VERSION_MAJOR 0)
set(BSG_VERSION_MINOR 11)
set(BSG_VERSION
  ${BSG_VERION_MAJOR}.${BSG_VERSION_MINOR})
 
# Offer the user the choice of overriding the installation directories
set(INSTALL_LIB_DIR ${PROJECT_BINARY_DIR}/lib CACHE PATH "Installation directory for libraries")
set(INSTALL_BIN_DIR ${PROJECT_BINARY_DIR}/bin CACHE PATH "Installation directory for executables")
set(INSTALL_INCLUDE_DIR ${PROJECT_BINARY_DIR}/include CACHE PATH
  "Installation directory for header files")
if(WIN32 AND NOT CYGWIN)
  set(DEF_INSTALL_CMAKE_DIR CMake)
else()
  set(DEF_INSTALL_CMAKE_DIR ${PROJECT/BINARY/DIR}/lib/CMake/BSG)
endif()
set(INSTALL_CMAKE_DIR ${DEF_INSTALL_CMAKE_DIR} CACHE PATH
  "Installation directory for CMake files")
 
# Make relative paths absolute (needed later on)
foreach(p LIB BIN INCLUDE CMAKE)
  set(var INSTALL_${p}_DIR)
  if(NOT IS_ABSOLUTE "${${var}}")
    set(${var} "${CMAKE_INSTALL_PREFIX}/${${var}}")
  endif()
endforeach()

set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH};${PROJECT_SOURCE_DIR}/Utilities")

find_package(GSL REQUIRED)
include_directories(${GSL_INCLUDE_DIRS})
find_package(ROOT REQUIRED)
include_directories(${ROOT_INCLUDE_DIRS})
find_package(Boost COMPONENTS program_options REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

configure_file("${PROJECT_SOURCE_DIR}/BSGConfig.h.in" "${PROJECT_BINARY_DIR}/BSGConfig.h")
configure_file("${PROJECT_SOURCE_DIR}/Utilities/config.txt" "${PROJECT_BINARY_DIR}/NME_exec/config.txt")
configure_file("${PROJECT_SOURCE_DIR}/Utilities/config.txt" "${PROJECT_BINARY_DIR}/BSG_exec/config.txt")

include_directories("${PROJECT_BINARY_DIR}")
include_directories("${PROJECT_SOURCE_DIR}/NME/include")
include_directories("${PROJECT_SOURCE_DIR}/BSG/include")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -std=c++11")
 
# Add sub-directories
add_subdirectory(NME)
add_subdirectory(BSG)
add_subdirectory(BSG_exec)
add_subdirectory(NME_exec)
 
# The interesting stuff goes here
# ===============================
 
# Add all targets to the build-tree export set
export(TARGETS nme bsg BSG NME
  FILE "${PROJECT_BINARY_DIR}/BSGTargets.cmake")
 
# Export the package for use from the build-tree
# (this registers the build-tree with a global CMake-registry)
export(PACKAGE BSG)
 
# Create the BSGConfig.cmake and BSGConfigVersion files
file(RELATIVE_PATH REL_INCLUDE_DIR "${INSTALL_CMAKE_DIR}"
   "${INSTALL_INCLUDE_DIR}")
# ... for the build tree
set(CONF_INCLUDE_DIRS "${PROJECT_SOURCE_DIR}" "${PROJECT_BINARY_DIR}")
configure_file(BSGConfig.cmake.in
  "${PROJECT_BINARY_DIR}/BSGConfig.cmake" @ONLY)
# ... for the install tree
set(CONF_INCLUDE_DIRS "\${BSG_CMAKE_DIR}/${REL_INCLUDE_DIR}")
configure_file(BSGConfig.cmake.in
  "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/BSGConfig.cmake" @ONLY)
# ... for both
configure_file(BSGConfigVersion.cmake.in
  "${PROJECT_BINARY_DIR}/BSGConfigVersion.cmake" @ONLY)
 
# Install the BSGConfig.cmake and BSGConfigVersion.cmake
install(FILES
  "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/BSGConfig.cmake"
  "${PROJECT_BINARY_DIR}/BSGConfigVersion.cmake"
  DESTINATION "${INSTALL_CMAKE_DIR}" COMPONENT dev)
 
# Install the export set for use with the install-tree
install(EXPORT BSGTargets DESTINATION
  "${INSTALL_CMAKE_DIR}" COMPONENT dev)
