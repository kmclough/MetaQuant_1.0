cmake_minimum_required (VERSION 2.6)

project(metaquant)

set(${PROJECT_NAME}_VERSION_MAJOR 0)
set(${PROJECT_NAME}_VERSION_MINOR 1)
set(${PROJECT_NAME}_VERSION_PATCH 0)

if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
    set(CMAKE_CXX_FLAGS "-Wall -DLinux --std=gnu++0x -lrt -ltcmalloc")
else(CMAKE_SYSTEM_NAME STREQUAL "Linux")
    set(CMAKE_CXX_FLAGS "-Wall -DDarwin")
endif(CMAKE_SYSTEM_NAME STREQUAL "Linux")

set(CMAKE_CXX_FLAGS_DEBUG "-g ${CMAKE_CXX_FLAGS}")
set(CMAKE_CXX_FLAGS_RHDEBINFO "-O3 -g ${CMAKE_CXX_FLAGS}")
set(CMAKE_CXX_FLAGS_MINSIZEREL "-Os ${CMAKE_CXX_FLAGS}")

set(CMAKE_BUILD_TYPE Release)
set(Boost_USE_STATIC_LIBS ON)

find_package(Boost 1.56
    COMPONENTS
        thread
   	system
    	filesystem
    	program_options
    	date_time
        chrono
        atomic
    REQUIRED)

find_library(GPERFTOOLS_TCMALLOC_LIB tcmalloc)
if (GPERFTOOLS_TCMALLOC_LIB)
   message(STATUS "Found GPERFTOOLS_TCMALLOC: ${GPERFTOOLS_TCMALLOC_LIB}")
else (GPERFTOOLS_TCMALLOC_LIB)
   message(STATUS "Could NOT find GPERFTOOLS_TCMALLOC: Install to improve speed.")
endif(GPERFTOOLS_TCMALLOC_LIB)

include_directories(${Boost_INCLUDE_DIRS} "${CMAKE_CURRENT_SOURCE_DIR}/bamtools/include") 
set(PROTO_INT 0)

if(WIN32)
	set(CMAKE_CXX_FLAGS "/EHsc")
	set(WIN32_INT 1)
else(WIN32)
	find_package(ZLIB REQUIRED)
	set(WIN32_INT 0)
endif(WIN32)

configure_file (
  "${PROJECT_SOURCE_DIR}/src/config.h.in"
  "${PROJECT_SOURCE_DIR}/src/config.h"
  )


add_subdirectory(src)
