# - Configure and install BSG Data Libraries
# The following functions and macros are defined to help define, install,
# reuse and export datasets and associated variables.
#
# FUNCTIONS AND MACROS
# ====================
# Public API
# ----------
# function bsg_get_datasetnames(<output variable>)
#          Store list of currently known dataset names in output variable
#
# function bsg_tupleize_datasets(<output variable>)
#          Set output variable to list of old-style dataset tuples.
#          A tuple has the format:
#            NAME/VERSION/FILENAME/EXTENSION/ENVVAR/MD5SUM
#          Provided for backward compatibility.
#
# function bsg_export_datasets([BUILD|INSTALL] <output variable>)
#          Set output variable to list of dataset tuples for export to
#          configuration scripts
#          A tuple has the format:
#            NAME/ENVVAR/PATH/FILENAME/MD5SUM
#          BUILD will set the PATH entry to the path to the dataset used
#          for the build of BSG.
#
#          INSTALL will set the PATH entry to the path to the dataset used
#          by an install of BSG.
#
# function bsg_add_dataset(NAME      <id>
#                             VERSION   <ver>
#                             FILENAME  <file>
#                             EXTENSION <ext>
#                             ENVVAR    <varname>
#                             MD5SUM    <md5>)
#          Define a new dataset with the following properties
#            NAME         common name of the dataset
#            VERSION      dataset version
#            FILENAME     name of packaged dataset file
#            EXTENSION    extension of packaged dataset file
#            ENVVAR       environment variable used by BSG code
#                         to locate this dataset
#            MD5SUM       MD5 hash of packaged dataset file
#
# function bsg_has_dataset(<name> <output variable>)
#          Set output variable to TRUE if the named dataset exists
#          in the list of known datasets.
#
# function bsg_get_dataset_property(<name> <property> <output variable>)
#          Set output variable to value of dataset "name"'s property.
#          If property does not exist, then value will be blank.
#
# function bsg_set_dataset_property(<name> <property> <value>)
#          Set value of dataset property to supplied value
#
# function bsg_configure_datasets(DESTINATION <dir>
#                                    DOWNLOAD    <installmissing>
#                                    TIMEOUT     <seconds>
#                                    )
#          Perform the actual heavy lifting to configure each declared
#          dataset for reuse or download as needed.
#
# function bsg_dataset_isinstalled(<name>
#                                     <root directory>
#                                     <output variable>)
#          Check if named dataset is installed under the root directory.
#          Set output variable to TRUE if it is, FALSE otherwise.
#
# function bsg_install_dataset(<name> <destination> <timeout>)
#          Download dataset with name to build directory, timing out the
#          download after timeout seconds, and install it
#          into its directory under the destination directory.
#
# function bsg_reuse_dataset(<name> <destination> <is installed>)
#          Reuse the dataset with name located at destination directory.
#          If it is not installed, warn user that it will need installing
#          manually in destination.
#

#-----------------------------------------------------------------------
# BSG PHYSICS DATA - GLOBAL CMAKE VARIABLES
#-----------------------------------------------------------------------
# URLs, directories and dataset entries
# We may want these as properties so we can have a small API for
# retrieving them globally
#-----------------------------------------------------------------------
# Where to install data in the build tree
set(BSG_BUILD_FULL_DATADIR ${PROJECT_BINARY_DIR}/data/databases)

# Where to install data in the install tree (a Default)
set(BSG_INSTALL_DATADIR_DEFAULT "${CMAKE_INSTALL_PREFIX}/data/databases")

# File containing dataset list
set(BSG_DATASETS_DEFINITIONS "BSGDatasetDefinitions")


#-----------------------------------------------------------------------
# BSG PHYSICS DATA - PUBLIC CMAKE API FOR DATASET HANDLING
#-----------------------------------------------------------------------
# Properties? Shouldn't clash with the tuplized variable...
define_property(GLOBAL PROPERTY "BSG_DATASETS"
  BRIEF_DOCS "List of all defined BSG dataset names"
  FULL_DOCS
  "Each element of the list gives the name defined for the dataset.
   This name can be used in other BSG Data API functions to
   extract other properties of the dataset"
  )

#-----------------------------------------------------------------------
# function bsg_get_datasetnames(<output variable>)
#          Store list of currently known dataset names in output variable
#
function(bsg_get_datasetnames _output)
  get_property(_tmp GLOBAL PROPERTY BSG_DATASETS)
  set(${_output} ${_tmp} PARENT_SCOPE)
endfunction()

#-----------------------------------------------------------------------
# function bsg_export_datasets([BUILD|INSTALL] <output variable>)
#          Set output variable to list of dataset tuples for export to
#          BSGConfig.cmake
#          A tuple has the format:
#            NAME/ENVVAR/PATH/FILENAME/MD5SUM
#          BUILD will set the PATH entry to the path to the dataset used
#          for the build of BSG.
#
#          INSTALL will set the PATH entry to the path to the dataset used
#          by an install of BSG.
#
#
# function(bsg_export_datasets _type _output)
#   bsg_get_datasetnames(_names)
#   set(_tmplist)
#
#   foreach(_ds ${_names})
#     set(_tuple ${_ds})
#     get_property(_tmpprop GLOBAL PROPERTY ${_ds}_ENVVAR)
#     list(APPEND _tuple ${_tmpprop})
#
#     if(${_type} STREQUAL "BUILD")
#       get_property(_tmpprop GLOBAL PROPERTY ${_ds}_BUILD_DIR)
#     elseif(${_type} STREQUAL "INSTALL")
#       get_property(_tmpprop GLOBAL PROPERTY ${_ds}_INSTALL_DIR)
#     else()
#       message(FATAL_ERROR "incorrect argument to bsg_export_datasets")
#     endif()
#     # Ensure CMake paths
#     file(TO_CMAKE_PATH "${_tmpprop}" _tmpprop)
#     list(APPEND _tuple ${_tmpprop})
#
#     get_property(_fname GLOBAL PROPERTY ${_ds}_FILENAME)
#     get_property(_fvers GLOBAL PROPERTY ${_ds}_VERSION)
#     get_property(_fextn GLOBAL PROPERTY ${_ds}_EXTENSION)
#     list(APPEND _tuple "${_fname}.${_fvers}.${_fextn}")
#
#     get_property(_tmpprop GLOBAL PROPERTY ${_ds}_MD5SUM)
#     list(APPEND _tuple ${_tmpprop})
#
#     # Because we have paths, use tuple separator that should not
#     # appear in a path.
#     string(REPLACE ";" "|" _tuple "${_tuple}")
#     list(APPEND _tmplist "${_tuple}")
#   endforeach()
#
#   set(${_output} ${_tmplist} PARENT_SCOPE)
# endfunction()

#-----------------------------------------------------------------------
# function bsg_add_dataset(NAME      <id>
#                             VERSION   <ver>
#                             URL  <url>
#                             ENVVAR    <varname>)
#          Define a new dataset with the following properties
#            NAME         common name of the dataset
#            VERSION      dataset version
#            URL          url of the data file
#            ENVVAR       environment variable used by BSG code
#                         to locate this dataset
#
function(bsg_add_dataset)
  # - Parse arguments and create variables
  set(oneValueArgs NAME VERSION URL ENVVAR NO_EXTRACT)
  cmake_parse_arguments(_BSGADDDATA "" "${oneValueArgs}" "" ${ARGN})

  # - Fail if already defined
  bsg_has_dataset(${_BSGADDDATA_NAME} _dsexists)
  if(_dsexists)
    message(FATAL_ERROR "Dataset ${_BSGADDDATA_NAME} is already defined")
  endif()

  # - Set properties as global props <NAME>_<PROP>
  # Simple properties...
  foreach(_prop VERSION URL ENVVAR NO_EXTRACT)
    set_property(GLOBAL PROPERTY ${_BSGADDDATA_NAME}_${_prop} ${_BSGADDDATA_${_prop}})
  endforeach()

  # DIRECTORY : the name of the directory that results from unpacking
  #             the packed dataset file.
  set_property(GLOBAL PROPERTY ${_BSGADDDATA_NAME}_DIRECTORY
    "${_BSGADDDATA_NAME}${_BSGADDDATA_VERSION}"
    )

  # - add it to the list of defined datasets
  set_property(GLOBAL APPEND PROPERTY BSG_DATASETS ${_BSGADDDATA_NAME})
endfunction()

#-----------------------------------------------------------------------
# function bsg_has_dataset(<name> <output variable>)
#          Set output variable to TRUE if the named dataset exists
#          in the list of known datasets.
#
function(bsg_has_dataset _name _output)
  get_property(_dslist GLOBAL PROPERTY BSG_DATASETS)
  list(FIND _dslist ${_name} _index)
  if(_index GREATER -1)
    set(${_output} TRUE PARENT_SCOPE)
  else()
    set(${_output} FALSE PARENT_SCOPE)
  endif()
endfunction()

#-----------------------------------------------------------------------
# function bsg_get_dataset_property(<name> <property> <output variable>)
#          Set output variable to value of dataset "name"'s property.
#          If property does not exist, then value will be blank.
#
function(bsg_get_dataset_property _name _prop _output)
  bsg_has_dataset(${_name} _dsexists)
  if(NOT _dsexists)
    message(FATAL_ERROR "non-existent dataset ${_name}")
  endif()

  get_property(_tmp GLOBAL PROPERTY ${_name}_${_prop})
  set(${_output} ${_tmp} PARENT_SCOPE)
endfunction()

#-----------------------------------------------------------------------
# function bsg_set_dataset_property(<name> <property> <value>)
#          Set value of dataset property to supplied value
#
function(bsg_set_dataset_property _name _prop _value)
  bsg_has_dataset(${_name} _dsexists)
  if(NOT _dsexists)
    message(FATAL_ERROR "non-existent dataset ${_name}")
  endif()
  set_property(GLOBAL PROPERTY ${_name}_${_prop} "${_value}")
endfunction()

#-----------------------------------------------------------------------
# function bsg_configure_datasets(DESTINATION <dir>
#                                    DOWNLOAD    <installmissing>
#                                    TIMEOUT     <seconds>
#                                    )
#          Perform the actual heavy lifting to configure each declared
#          dataset for reuse or download as needed.
#
function(bsg_configure_datasets)
  # - Parse arguments and create variables
  set(oneValueArgs DESTINATION DOWNLOAD TIMEOUT)
  cmake_parse_arguments(_BSGCFGDSS "" "${oneValueArgs}" "" ${ARGN})

  # - Load configuration
  include(${BSG_DATASETS_DEFINITIONS})
  bsg_get_datasetnames(_dsnames)
  set(_notinstalled )

  foreach(_ds ${_dsnames})
    bsg_dataset_isinstalled(${_ds} "${_BSGCFGDSS_DESTINATION}" _installed)
    if(NOT _installed AND _BSGCFGDSS_DOWNLOAD)
      bsg_install_dataset(${_ds} "${_BSGCFGDSS_DESTINATION}" ${_BSGCFGDSS_TIMEOUT})
    else()
      bsg_reuse_dataset(${_ds} "${_BSGCFGDSS_DESTINATION}" ${_installed})
      if(NOT _installed)
        list(APPEND _notinstalled ${_ds})
      endif()
    endif()
  endforeach()

  # - Produce report on datasets needing manual install, advising
  # user on how to handle these.
  # Yes, it's long, but at least it's clear :-)
  if(_notinstalled)
    message("  *WARNING*")
    message("    BSG has been pre-configured to look for datasets")
    message("    in the directory:")
    message(" ")
    message("    ${_BSGCFGDSS_DESTINATION}")
    message(" ")
    message("    but the following datasets are NOT present on disk at")
    message("    that location:")
    message(" ")
    foreach(_missing ${_notinstalled})
      bsg_get_dataset_property(${_missing} VERSION _vers)
      message("    ${_missing} (${_vers})")
    endforeach()
    message(" ")
    message("    If you want to have these datasets installed automatically")
    message("    simply re-run cmake and set the BSG_INSTALL_DATA")
    message("    variable to ON. This will configure the build to download")
    message("    and install these datasets for you. For example, on the")
    message("    command line, do:")
    message(" ")
    message("    cmake -DBSG_INSTALL_DATA=ON <otherargs>")
    message(" ")
    message("    The variable can also be toggled in ccmake or cmake-gui.")
    message("")
    message("    If you're running on a Windows system, this is the best")
    message("    solution as CMake will unpack the datasets for you")
    message("    without any further software being required")
    message(" ")
    message("    Alternatively, you can install these datasets manually")
    message("    now or after you have installed BSG. To do this,")
    message("    download the following files:")
    message(" ")
    foreach(_missing ${_notinstalled})
      bsg_get_dataset_property(${_missing} URL _url)
      message("    ${_url}")
    endforeach()
    message(" ")
    message("    and unpack them under the directory:")
    message(" ")
    message("    ${_BSGCFGDSS_DESTINATION}")
    message(" ")
    message("    As we supply the datasets packed in gzipped tar files,")
    message("    you will need the 'tar' utility to unpack them.")
    message(" ")
    message("    Nota bene: Missing datasets will not affect or break")
    message("               compilation and installation of the BSG")
    message("               libraries.")
    message(" ")
  endif()
endfunction()

#-----------------------------------------------------------------------
# function bsg_dataset_isinstalled(<name>
#                                     <root directory>
#                                     <output variable>)
#          Check if named dataset is installed under the root directory.
#          Set output variable to TRUE if it is, FALSE otherwise.
#
function(bsg_dataset_isinstalled _name _rdirectory _output)
  bsg_get_dataset_property(${_name} DIRECTORY _dsdir)
  set(_expectedpath ${_rdirectory}/${_dsdir})

  if(IS_DIRECTORY ${_expectedpath})
    set(${_output} TRUE PARENT_SCOPE)
  else()
    set(${_output} FALSE PARENT_SCOPE)
  endif()
endfunction()

#-----------------------------------------------------------------------
# function bsg_install_dataset(<name> <destination> <timeout>)
#          Download dataset with name to build directory, timing out the
#          download after timeout seconds, and install it
#          into its directory under the destination directory.
#
function(bsg_install_dataset _name _destination _timeout)
  # - Extract needed dataset properties
  bsg_get_dataset_property(${_name} DIRECTORY _ds_dir)
  bsg_get_dataset_property(${_name} VERSION _ds_version)
  bsg_get_dataset_property(${_name} URL _ds_url)
  bsg_get_dataset_property(${_name} NO_EXTRACT _ds_extract)

  message(STATUS "Configuring download of missing dataset ${_name} (${_ds_version})")

  if(${_ds_extract})
    file(DOWNLOAD ${_ds_url} ${BSG_BUILD_FULL_DATADIR}/${_ds_dir}/${_name}.txt SHOW_PROGRESS)
    message(STATUS "Downloaded ${_name}.txt to ${BSG_BUILD_FULL_DATADIR}/${_ds_dir}.")
  else()
    include(ExternalProject)
    ExternalProject_Add(${_name}
      PREFIX Externals/${_name}-${_ds_version}
      SOURCE_DIR ${BSG_BUILD_FULL_DATADIR}/${_ds_dir}
      URL ${_ds_url}
      TIMEOUT ${_timeout}
      CONFIGURE_COMMAND ""
      BUILD_COMMAND ""
      INSTALL_COMMAND ""
      )
  endif()

  # - Configure the dataset's build and install locations
  bsg_set_dataset_property(${_name} BUILD_DIR "${PROJECT_BINARY_DIR}/data/databases/${_ds_dir}")
  bsg_set_dataset_property(${_name} INSTALL_DIR "${_destination}/${_ds_dir}")

  # - Add install target, and report back paths...
  install(DIRECTORY ${PROJECT_BINARY_DIR}/data/databases/${_ds_dir}
    DESTINATION ${_destination}
    COMPONENT Data
    )
endfunction()

#-----------------------------------------------------------------------
# function bsg_reuse_dataset(<name> <destination> <is installed>)
#          Reuse the dataset with name located at destination directory.
#          If it is not installed, warn user that it will need installing
#          manually in destination.
#
function(bsg_reuse_dataset _name _destination _ispresent)
  bsg_get_dataset_property(${_name} VERSION _ds_ver)
  bsg_get_dataset_property(${_name} DIRECTORY _ds_dir)
  if(_ispresent)
    message(STATUS "Reusing dataset ${_name} (${_ds_ver})")
  else()
    message(STATUS "Pre-configuring dataset ${_name} (${_ds_ver})")
  endif()

  # - In both cases, the build and install dirs are identical
  bsg_set_dataset_property(${_name} BUILD_DIR "${_destination}/${_ds_dir}")
  bsg_set_dataset_property(${_name} INSTALL_DIR "${_destination}/${_ds_dir}")
endfunction()


#-----------------------------------------------------------------------
# BSG PHYSICS DATA - USER INTERFACE AND PROCESSING
#-----------------------------------------------------------------------
# User options for installing data
# - Choose a directory under which to install the data.
# - Choose whether to download and install missing datasets.
# - Change download timeout for problematic connections
#-----------------------------------------------------------------------
# Choose Physics Data Install Dir
# This follows the pattern for interface and setting as in GNUInstallDirs
if(NOT BSG_INSTALL_DATADIR)
  set(BSG_INSTALL_DATADIR "" CACHE PATH "read-only architecture independent BSG physics data (DATAROOTDIR/${BSG_INSTALL_DATADIR_DEFAULT}")
  set(BSG_INSTALL_DATADIR "${BSG_INSTALL_DATADIR_DEFAULT}")
endif()

if(NOT IS_ABSOLUTE ${BSG_INSTALL_DATADIR})
  set(BSG_INSTALL_FULL_DATADIR "${CMAKE_INSTALL_PREFIX}/${BSG_INSTALL_DATADIR}")
else()
  set(BSG_INSTALL_FULL_DATADIR "${BSG_INSTALL_DATADIR}")
endif()

#-----------------------------------------------------------------------
# Select whether to download and install missing datasets
option(BSG_INSTALL_DATA "Download/Install datasets missing from BSG_INSTALL_DATADIR" OFF)

#-----------------------------------------------------------------------
# Provide an option for increasing the download timeout
# Helps with large datasets over slow connections.
set(BSG_INSTALL_DATA_TIMEOUT 1500 CACHE STRING "Timeout for Data Library downloads")
mark_as_advanced(BSG_INSTALL_DATA_TIMEOUT)

#-----------------------------------------------------------------------
# Set up check, download and install of needed data
#
bsg_configure_datasets(
  DESTINATION ${BSG_INSTALL_FULL_DATADIR}
  DOWNLOAD    ${BSG_INSTALL_DATA}
  TIMEOUT     ${BSG_INSTALL_DATA_TIMEOUT}
  )
