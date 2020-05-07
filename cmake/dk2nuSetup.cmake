if(DEFINED USE_DK2NU AND USE_DK2NU)

  #Check if it already set up (probably by ups)
  if(DEFINED ENV{DK2NU_LIB} AND DEFINED ENV{DK2NU_INC})
    LIST(APPEND INCDIRS $ENV{DK2NU_INC})
    LIST(APPEND EXTRA_LINK_DIRS $ENV{DK2NU_LIB})
    LIST(APPEND EXTRA_LIBS dk2nuTree)
  else()

  SET(DK2NU_EXTRA_CXX_FLAGS)
  SET(WITH_GENIE,false)
    if(DEFINED ENV{GENIE})
      SET(WITH_GENIE,true)

      include(${CMAKE_SOURCE_DIR}/cmake/parseConfigApp.cmake)

      ################################  LIBXML  ######################################
      if("${LIBXML2_LIB} " STREQUAL " ")
        GETLIBDIR(xml2-config --libs LIBXML2_LIB IGNORE_EMPTY_RESPONSE)
        if("${LIBXML2_LIB}  " STREQUAL " ")
          message(WARNING "Variable LIBXML2_LIB is not defined, as xml2-config was found and didn't report a library include path, it is likely that libxml2.so can be found in the standard system location, lets hope so. Alternativly, a location can be forced by configering with -DLIBXML2_LIB=/path/to/LIBXML2_libraries or as an environment variable LIBXML2_LIB.")
        endif()
      endif()
      if(NOT "${LIBXML2_LIB} " STREQUAL " ")
        SET(DK2NU_EXTRA_CXX_FLAGS "${DK2NU_EXTRA_CXX_FLAGS} -L${LIBXML2_LIB}")
      endif()

      if("${LIBXML2_INC} " STREQUAL " ")
        GETINCDIR(xml2-config --cflags LIBXML2_INC IGNORE_EMPTY_RESPONSE)
        if("${LIBXML2_INC} " STREQUAL " ")
          message(WARNING "Variable LIBXML2_INC is not defined, as xml2-config was found and didn't report an include path, it is likely that libxml2.so can be found in the standard system location, lets hope so. Alternativly, a location can be forced by configering with -DLIBXML2_INC=/path/to/LIBXML2_includes or as an environment variable LIBXML2_INC.")
        endif()
      endif()
      if(NOT "${LIBXML2_INC} " STREQUAL " ")
        SET(DK2NU_EXTRA_CXX_FLAGS "${DK2NU_EXTRA_CXX_FLAGS} -I${LIBXML2_INC}")
      endif()

      ###############################  log4cpp  ######################################
      if("${LOG4CPP_LIB} " STREQUAL " ")
        GETLIBDIR(log4cpp-config --libs LOG4CPP_LIB IGNORE_EMPTY_RESPONSE)
        if("${LOG4CPP_LIB} " STREQUAL " ")
          message(WARNING "Variable LOG4CPP_LIB is not defined, as xml2-config was found and didn't report a library include path, it is likely that liblog4cpp.so can be found in the standard system location, lets hope so. Alternativly, a location can be forced by configering with -DLOG4CPP_LIB=/path/to/LOG4CPP_libraries or as an environment variable LOG4CPP_LIB.")
        endif()
      endif()
      if(NOT "${LOG4CPP_LIB} " STREQUAL " ")
        SET(DK2NU_EXTRA_CXX_FLAGS "${DK2NU_EXTRA_CXX_FLAGS} -L${LOG4CPP_LIB}")
      endif()

      if("${LOG4CPP_INC} " STREQUAL " ")
        GETINCDIR(log4cpp-config --cflags LOG4CPP_INC IGNORE_EMPTY_RESPONSE)
        if("${LOG4CPP_INC} " STREQUAL " ")
          message(WARNING "Variable LOG4CPP_INC is not defined, as xml2-config was found and didn't report an include path, it is likely that log4cpp headers can be found in the standard system location, lets hope so. Alternativly, a location can be forced by configering with -DLOG4CPP_INC=/path/to/LOG4CPP_includes or as an environment variable LOG4CPP_INC.")
        endif()
      endif()
      if(NOT "${LOG4CPP_INC} " STREQUAL " ")
        SET(DK2NU_EXTRA_CXX_FLAGS "${DK2NU_EXTRA_CXX_FLAGS} -I${LOG4CPP_INC}")
      endif()

      ################################################################################


    endif()

    ###### dk2nu set up
    ExternalProject_Add(dk2nu
      PREFIX "${PROJECT_BINARY_DIR}/Ext"
      SVN_REPOSITORY https://cdcvs.fnal.gov/subversion/dk2nu/tags/v01_05_01
      UPDATE_DISCONNECTED 1
      CONFIGURE_COMMAND ${CMAKE_COMMAND}
        ${PROJECT_BINARY_DIR}/Ext/src/dk2nu/dk2nu
        -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}
        -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
        -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
        -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
        -DCMAKE_CXX_FLAGS=${DK2NU_EXTRA_CXX_FLAGS}
        -DWITH_GENIE=${WITH_GENIE} -DWITH_TBB=OFF
      )

    LIST(APPEND INCDIRS ${PROJECT_BINARY_DIR}/Ext/src/dk2nu/)
    LIST(APPEND EXTRA_LIBS dk2nuTree)
  endif()

  LIST(APPEND EXTRA_CXX_FLAGS -DUSE_DK2NU)

  SET(USE_DK2NU TRUE)
else()
  SET(USE_DK2NU FALSE)
endif()
