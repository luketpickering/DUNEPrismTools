SET(SIM_APPS
  Condenser
  StopProcessor
  CheckRunPlanCoverage)

include_directories(${INCDIRS})

foreach(appname ${SIM_APPS})

  add_executable(dp_${appname} ${appname}.cxx)

  target_link_libraries(dp_${appname} ${MODULETargets})
  target_link_libraries(dp_${appname} ${ROOT_LIBS})

  if(NOT "${CMAKE_LINK_FLAGS}" STREQUAL "")
    set_target_properties(dp_${appname} PROPERTIES LINK_FLAGS ${CMAKE_LINK_FLAGS})
  endif()

  if(NOT "${CMAKE_DEPENDLIB_FLAGS}" STREQUAL "")
    target_link_libraries(dp_${appname} ${CMAKE_DEPENDLIB_FLAGS})
  endif()

  install(TARGETS dp_${appname} DESTINATION bin)

endforeach()
