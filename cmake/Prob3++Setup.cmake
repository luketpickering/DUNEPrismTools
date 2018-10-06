ExternalProject_Add(prob3pp
  PREFIX "${CMAKE_BINARY_DIR}/Ext"
  URL "http://webhome.phy.duke.edu/~raw22/public/Prob3++/Prob3++.20121225.tar.gz"
  CONFIGURE_COMMAND ""
  BUILD_IN_SOURCE 1
  UPDATE_COMMAND ""
  BUILD_COMMAND make
  INSTALL_COMMAND ""
  )

LIST(APPEND INCDIRS ${CMAKE_BINARY_DIR}/Ext/src/prob3pp)

LIST(APPEND EXTRA_LINK_DIRS ${CMAKE_BINARY_DIR}/Ext/src/prob3pp)

LIST(APPEND EXTRA_LIBS ThreeProb_2.10)

cmessage(STATUS "Using Prob3++ 2.10")
