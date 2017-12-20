SET(APPS BuildFluxes
  FitFluxes
  OscillateFlux
  PredictRunPlanEventRates
  CombineBuiltFluxes
  BuildUncertaintyMatrix)

  include_directories(${INCLUDE_DIRECTORIES})

foreach(a ${APPS})

  add_executable(dp_${a} ${a}.cxx)

  target_link_libraries(dp_${a} ${MODULETargets})
  target_link_libraries(dp_${a} ${ROOT_LIBS})


  if(NOT "${CMAKE_LINK_FLAGS}" STREQUAL "")
    set_target_properties(dp_${a} PROPERTIES LINK_FLAGS ${CMAKE_LINK_FLAGS})
  endif()

  if(NOT "${CMAKE_DEPENDLIB_FLAGS}" STREQUAL "")
    target_link_libraries(dp_${a} ${CMAKE_DEPENDLIB_FLAGS})
  endif()

  install(TARGETS dp_${a} DESTINATION bin)


endforeach()
