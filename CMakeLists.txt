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

cmake_minimum_required(VERSION 2.8)

project(BSG)

# - Prepend our own CMake Modules to the search path
# NB: if our custom modules include others that we don't supply, those in
# the base path will be used, so watch for incompatibilities!!
#
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/Modules ${CMAKE_MODULE_PATH})

include(CMakeUninstallTarget)

# Offer the user the choice of overriding the installation directories
set(INSTALL_LIB_DIR ${PROJECT_BINARY_DIR}/lib CACHE PATH "Installation directory for libraries")
set(INSTALL_BIN_DIR ${PROJECT_BINARY_DIR}/bin CACHE PATH "Installation directory for executables")
set(INSTALL_INCLUDE_DIR ${PROJECT_BINARY_DIR}/include CACHE PATH
  "Installation directory for header files")

add_subdirectory(source)

include(BSGInstallData)

# Install the targets file
install(EXPORT bsg-targets DESTINATION lib/bsg)
# TODO Make the project usable from the build-tree by creating a file explicitly
export(TARGETS bsg nme bsg_exec nme_exec FILE bsg-exports.cmake)
#configure_file(${PROJECT_SOURCE_DIR}/bsg-config.cmake.in ${PROJECT_BINARY_DIR}/bsg-config.cmake @ONLY)

# Make files for finding the package later using find_package in the install tree
configure_file(${PROJECT_SOURCE_DIR}/cmake/Templates/bsg-config.cmake.in ${PROJECT_BINARY_DIR}/config/bsg-config.cmake @ONLY)
configure_file(${PROJECT_SOURCE_DIR}/cmake/Templates/bsg-config-version.cmake.in ${PROJECT_BINARY_DIR}/config/bsg-config-version.cmake @ONLY)
install(FILES ${PROJECT_BINARY_DIR}/config/bsg-config.cmake ${PROJECT_BINARY_DIR}/config/bsg-config-version.cmake
DESTINATION lib/bsg)

install(DIRECTORY ${PROJECT_SOURCE_DIR}/data/ DESTINATION data)

# Copy over the documentation
install(DIRECTORY ${PROJECT_SOURCE_DIR}/docs/doxygen/html/ DESTINATION doc/bsg)
install(DIRECTORY ${PROJECT_SOURCE_DIR}/docs/build/html/ DESTINATION doc/bsg)

install(FILES ${PROJECT_BINARY_DIR}/BSGConfig.h ${PROJECT_BINARY_DIR}/NMEConfig.h DESTINATION include)
