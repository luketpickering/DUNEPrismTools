SET(TEST_APPS
  MatrixToolsTester)

include_directories(${INCDIRS})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

foreach(appname ${TEST_APPS})

  add_executable(dpt_${appname} ${appname}.cxx)

  target_link_libraries(dpt_${appname} ${MODULETargets})
  target_link_libraries(dpt_${appname} ${ROOT_LIBS})

  if(NOT "${CMAKE_LINK_FLAGS}" STREQUAL "")
    set_target_properties(dpt_${appname} PROPERTIES LINK_FLAGS ${CMAKE_LINK_FLAGS})
  endif()

  if(NOT "${CMAKE_DEPENDLIB_FLAGS}" STREQUAL "")
    target_link_libraries(dpt_${appname} ${CMAKE_DEPENDLIB_FLAGS})
  endif()

  install(TARGETS dpt_${appname} DESTINATION bin)

endforeach()
