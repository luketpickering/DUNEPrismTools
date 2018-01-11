set(CXX_WARNINGS -Wall )

if(CMAKE_CXX_COMPILER_VERSION VERSION_GREATER 5.0)
  LIST(APPEND EXTRA_CXX_FLAGS -std=c++1y)
elseif(CMAKE_CXX_COMPILER_VERSION VERSION_GREATER 4.8)
  LIST(APPEND EXTRA_CXX_FLAGS -std=c++11)
endif()

cmessage(DEBUG "EXTRA_CXX_FLAGS: ${EXTRA_CXX_FLAGS}")
string(REPLACE ";" " " STR_EXTRA_CXX_FLAGS "${EXTRA_CXX_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${STR_EXTRA_CXX_FLAGS} ${CXX_WARNINGS}")

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -fPIC -O3")

if(CMAKE_BUILD_TYPE MATCHES DEBUG)
  set(CURRENT_CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS_DEBUG})
elseif(CMAKE_BUILD_TYPE MATCHES RELEASE)
  set(CURRENT_CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS_RELEASE})
else()
  cmessage(FATAL_ERROR "[ERROR]: Unknown CMAKE_BUILD_TYPE (\"${CMAKE_BUILD_TYPE}\"): Should be \"DEBUG\" or \"RELEASE\".")
endif()


SET(STR_EXTRA_LINK_DIRS)
if(NOT EXTRA_LINK_DIRS STREQUAL "")
  string(REPLACE ";" " -L" STR_EXTRA_LINK_DIRS "-L${EXTRA_LINK_DIRS}")
endif()
SET(STR_EXTRA_LIBS)
if(NOT EXTRA_LIBS STREQUAL "")
  string(REPLACE ";" " -l" STR_EXTRA_LIBS "-l${EXTRA_LIBS}")
endif()
SET(STR_EXTRA_SHAREDOBJS)
if(NOT EXTRA_SHAREDOBJS STREQUAL "")
  string(REPLACE ";" " " STR_EXTRA_SHAREDOBJS "${EXTRA_SHAREDOBJS}")
endif()
SET(STR_EXTRA_LINK_FLAGS)
if(NOT EXTRA_LINK_FLAGS STREQUAL "")
  string(REPLACE ";" " " STR_EXTRA_LINK_FLAGS "${EXTRA_LINK_FLAGS}")
endif()

cmessage(DEBUG "EXTRA_LINK_DIRS: ${STR_EXTRA_LINK_DIRS}")
cmessage(DEBUG "EXTRA_LIBS: ${STR_EXTRA_LIBS}")
cmessage(DEBUG "EXTRA_SHAREDOBJS: ${STR_EXTRA_SHAREDOBJS}")
cmessage(DEBUG "EXTRA_LINK_FLAGS: ${STR_EXTRA_LINK_FLAGS}")

if(NOT STR_EXTRA_LINK_DIRS STREQUAL "" AND NOT STR_EXTRA_LIBS STREQUAL "")
  SET(CMAKE_DEPENDLIB_FLAGS "${STR_EXTRA_LINK_DIRS} ${STR_EXTRA_LIBS}")
endif()

if(NOT EXTRA_SHAREDOBJS STREQUAL "")
  if(NOT STR_EXTRA_LINK_FLAGS STREQUAL "")
    SET(STR_EXTRA_LINK_FLAGS "${STR_EXTRA_SHAREDOBJS} ${STR_EXTRA_LINK_FLAGS}")
  else()
    SET(STR_EXTRA_LINK_FLAGS "${STR_EXTRA_SHAREDOBJS}")
  endif()
endif()

if(NOT EXTRA_LINK_FLAGS STREQUAL "")
  if(NOT CMAKE_LINK_FLAGS STREQUAL "")
    SET(CMAKE_LINK_FLAGS "${CMAKE_LINK_FLAGS} ${STR_EXTRA_LINK_FLAGS}")
  else()
    SET(CMAKE_LINK_FLAGS "${STR_EXTRA_LINK_FLAGS}")
  endif()
endif()

if (VERBOSE)
  cmessage (STATUS "C++ Compiler      : ${CXX_COMPILER_NAME}")
  cmessage (STATUS "    flags         : ${CMAKE_CXX_FLAGS}")
  cmessage (STATUS "    Release flags : ${CMAKE_CXX_FLAGS_RELEASE}")
  cmessage (STATUS "    Debug flags   : ${CMAKE_CXX_FLAGS_DEBUG}")
  cmessage (STATUS "    Link Flags    : ${CMAKE_LINK_FLAGS}")
  cmessage (STATUS "    Lib Flags     : ${CMAKE_DEPENDLIB_FLAGS}")
endif()
