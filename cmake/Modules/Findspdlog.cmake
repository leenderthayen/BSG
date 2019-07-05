# Find path to the spdlog library (https://github.com/gabime/spdlog)

# This module defines
# spdlod_INCLUDE_DIRS, where to find the headers
# spd_FOUND, whether or not the library was found

find_path(spdlog_INCLUDE_DIRS spdlog/spdlog.h
    ~/spdlog/include
    /usr/local/include
    /usr/include
    ~/git/spdlog/include
  ${PROJECT_SOURCE_DIR})
if(${spdlog_INCLUDE_DIRS})
	set(spdlog_FOUND "YES")
endif()
