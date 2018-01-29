set(BSG_VERSION_MAJOR 0)
set(BSG_VERSION_MINOR 13)
set(BSG_LAST_UPDATE "January 24, 2018")

set(NME_VERSION_MAJOR 0)
set(NME_VERSION_MINOR 10)
set(NME_LAST_UPDATE "January 24, 2018")
 
find_package(GSL REQUIRED)
include_directories(${GSL_INCLUDE_DIRS})
find_package(ROOT REQUIRED)
include_directories(${ROOT_INCLUDE_DIRS})
find_package(Boost COMPONENTS program_options REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})
find_package(Threads)

find_package(spdlog)
if(NOT ${spdlog_FOUND})
  message(STATUS "spdlog not found through find_package. Checking manually in PATH.")
  find_path(spdlog_INCLUDES spdlog.h PATHS
            /usr/local/include
            /usr/include)

  if(${SPDLOG_INCLUDES} STREQUAL "SPDLOG_INCLUDES-NOTFOUND")
  	message(FATAL_ERROR "The spdlog.h file was not found. Add it to your PATH.")
  endif()
endif()

include_directories(${spdlog_INCLUDES})

# Configure the version header files
configure_file("${PROJECT_SOURCE_DIR}/cmake/Templates/BSGConfig.h.in" "${PROJECT_BINARY_DIR}/BSGConfig.h")
configure_file("${PROJECT_SOURCE_DIR}/cmake/Templates/NMEConfig.h.in" "${PROJECT_BINARY_DIR}/NMEConfig.h")

include_directories("${PROJECT_SOURCE_DIR}/source/nme/include")
include_directories("${PROJECT_SOURCE_DIR}/source/bsg/include")

# Copy the data directory
#configure_file("${PROJECT_SOURCE_DIR}/Utilities/config.txt" "${PROJECT_BINARY_DIR}/NME_exec/config.txt")
#configure_file("${PROJECT_SOURCE_DIR}/Utilities/config.txt" "${PROJECT_BINARY_DIR}/BSG_exec/config.txt")
#configure_file("${PROJECT_SOURCE_DIR}/Utilities/data/ExchangeData.dat" "${PROJECT_BINARY_DIR}/BSG_exec/ExchangeData.dat")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -std=c++11")
 
# Add sub-directories
include_directories("${PROJECT_BINARY_DIR}")
add_subdirectory(nme)
add_subdirectory(bsg)
add_subdirectory(executables)

export(TARGETS nme bsg bsg_exec nme_exec FILE ${PROJECT_BINARY_DIR}/bsg-exports.cmake)

install(FILES ${PROJECT_BINARY_DIR}/BSGConfig.h ${PROJECT_BINARY_DIR}/NMEConfig.h DESTINATION include)
 
install(EXPORT BSGTargets DESTINATION lib)

# Export the package for use from the build-tree
# (this registers the build-tree with a global CMake-registry)
#export(PACKAGE BSG)
 
# Create the BSGConfig.cmake and BSGConfigVersion files
#file(RELATIVE_PATH REL_INCLUDE_DIR "${INSTALL_CMAKE_DIR}"
#   "${INSTALL_INCLUDE_DIR}")
# ... for the build tree
#set(CONF_INCLUDE_DIRS "${PROJECT_SOURCE_DIR}" "${PROJECT_BINARY_DIR}")
#configure_file(BSGConfig.cmake.in
#  "${PROJECT_BINARY_DIR}/BSGConfig.cmake" @ONLY)
# ... for the install tree
#set(CONF_INCLUDE_DIRS "\${BSG_CMAKE_DIR}/${REL_INCLUDE_DIR}")
#configure_file(BSGConfig.cmake.in
#  "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/BSGConfig.cmake" @ONLY)
# ... for both
#configure_file(BSGConfigVersion.cmake.in
#  "${PROJECT_BINARY_DIR}/BSGConfigVersion.cmake" @ONLY)
 
# Install the BSGConfig.cmake and BSGConfigVersion.cmake
#install(FILES
#  "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/BSGConfig.cmake"
#  "${PROJECT_BINARY_DIR}/BSGConfigVersion.cmake"
#  DESTINATION "${INSTALL_CMAKE_DIR}" COMPONENT dev)
 
# Install the export set for use with the install-tree
#install(EXPORT BSGTargets DESTINATION
#  "${INSTALL_CMAKE_DIR}" COMPONENT dev)
