cmake_minimum_required (VERSION 2.6 FATAL_ERROR)

project(DUNEPrismTools)

include(ExternalProject)

set (DUNEPrismTools_VERSION_MAJOR 0)
set (DUNEPrismTools_VERSION_MINOR 3)
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

if(USE_FHICL STREQUAL "")
  set(USE_FHICL TRUE)
elseif(NOT DEFINED USE_FHICL)
  set(USE_FHICL TRUE)
endif()

if(USE_EIGEN STREQUAL "")
  set(USE_EIGEN TRUE)
elseif(NOT DEFINED USE_EIGEN)
  set(USE_EIGEN TRUE)
endif()

################################################################################
#                            Check Dependencies
################################################################################

SET(EXTRA_CXX_FLAGS)
SET(INCDIRS)
SET(EXTRA_LINK_DIRS)
SET(EXTRA_LIBS)

##################################  ROOT  ######################################
include(${CMAKE_SOURCE_DIR}/cmake/ROOTSetup.cmake)

################################  Prob3++  #####################################
include(${CMAKE_SOURCE_DIR}/cmake/Prob3++Setup.cmake)

###############################  fhiclcpp  #####################################
include(${CMAKE_SOURCE_DIR}/cmake/fhiclcppSetup.cmake)

##########################  eigen and spectra  #################################
include(${CMAKE_SOURCE_DIR}/cmake/eigenSetup.cmake)

#################################  TH2Jagged  #################################
include(${CMAKE_SOURCE_DIR}/cmake/TH2JaggedSetup.cmake)

#################################  TH2Jagged  #################################
include(${CMAKE_SOURCE_DIR}/cmake/dk2nuSetup.cmake)

################################## COMPILER ####################################
LIST(APPEND EXTRA_CXX_FLAGS -fopenmp)
include(${CMAKE_SOURCE_DIR}/cmake/c++CompilerSetup.cmake)

################################## DOXYGEN ####################################
include(${CMAKE_SOURCE_DIR}/cmake/doxSetup.cmake)

add_subdirectory(src)

cmessage(DEBUG "SRC INCDIRS: ${INCDIRS}")
cmessage(DEBUG "Module Targets: ${MODULETargets}")

add_subdirectory(app)

add_subdirectory(configs)
add_subdirectory(scripts)

add_subdirectory(dox)

configure_file(cmake/setup.sh.in
  "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/setup.sh" @ONLY)
install(FILES
  "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/setup.sh" DESTINATION
  ${CMAKE_INSTALL_PREFIX})
