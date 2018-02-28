cmake_minimum_required (VERSION 2.6 FATAL_ERROR)

project(DUNEPrismTools)

include(ExternalProject)

set (DUNEPrismTools_VERSION_MAJOR 0)
set (DUNEPrismTools_VERSION_MINOR 2)
set (DUNEPrismTools_VERSION_REVISION 0)

set (DUNEPrismTools_VERSION_STRING "v${DUNEPrismTools_VERSION_MAJOR}r${DUNEPrismTools_VERSION_MINOR}p${DUNEPrismTools_VERSION_REVISION}")

#Set this to TRUE to enable build debugging messages
set(BUILD_DEBUG_MSGS TRUE)

include(${CMAKE_SOURCE_DIR}/cmake/cmessage.cmake)

#Changes default install path to be a subdirectory of the build dir.
#Can set build dir at configure time with -DCMAKE_INSTALL_PREFIX=/install/path
if(CMAKE_INSTALL_PREFIX STREQUAL "" OR CMAKE_INSTALL_PREFIX STREQUAL
  "/usr/local")
  set(CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}/${CMAKE_SYSTEM_NAME}")
elseif(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}/${CMAKE_SYSTEM_NAME}")
endif()

if(CMAKE_BUILD_TYPE STREQUAL "")
  set(CMAKE_BUILD_TYPE DEBUG)
elseif(NOT DEFINED CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE DEBUG)
endif()

################################################################################
#                            Check Dependencies
################################################################################

SET(EXTRA_CXX_FLAGS)
SET(INCLUDE_DIRECTORIES)
SET(EXTRA_LINK_DIRS)
SET(EXTRA_LIBS)

##################################  ROOT  ######################################
include(${CMAKE_SOURCE_DIR}/cmake/ROOTSetup.cmake)

################################  Prob3++  #####################################
include(${CMAKE_SOURCE_DIR}/cmake/Prob3++Setup.cmake)

################################## COMPILER ####################################
include(${CMAKE_SOURCE_DIR}/cmake/c++CompilerSetup.cmake)

################################## DOXYGEN ####################################
include(${CMAKE_SOURCE_DIR}/cmake/doxSetup.cmake)


LIST(APPEND INCLUDE_DIRECTORIES ${CMAKE_SOURCE_DIR}/src)

add_subdirectory(src/app)

add_subdirectory(G4ArSimulationTools)

configure_file(cmake/setup.sh.in
  "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/setup.sh" @ONLY)
install(FILES
  "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/setup.sh" DESTINATION
  ${CMAKE_INSTALL_PREFIX})

file(GLOB JOBSUB_SCRIPTS scripts/jobsub/*.sh)
install(PROGRAMS ${JOBSUB_SCRIPTS} DESTINATION ${CMAKE_INSTALL_PREFIX}/scripts)

file(GLOB CONFIG_FILES configs/*.*ml)
install(FILES ${CONFIG_FILES} DESTINATION ${CMAKE_INSTALL_PREFIX}/configs)
