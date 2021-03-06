SET(ANA_APPS
  # FluxVarFriendTreeBuilder
  # XSecVarFriendTreeBuilder
  # EfficiencyWeightFriendTreeBuilder
  # PRISMAnalysis
  # RunSelection
  # ProduceEfficiencyCorrector
  # ProduceDiagnosticPlots
  # CAFCombiner
  # ProcessedCAFDiagnosticPlots
  # ProcessedCAFDumpNominalEventRates
  # ProcessedCAFPRISMPrediction
  NDFDFluxRatioPlots
  # LinearCombinationUncertainties
  # gspl2root
)

include_directories(${INCDIRS})

foreach(appname ${ANA_APPS})

  add_executable(dp_${appname} ${appname}.cxx)

  target_link_libraries(dp_${appname} ${MODULETargets})
  target_link_libraries(dp_${appname} ${ROOT_LIBS})

  add_dependencies(dp_${appname} spectra)
  add_dependencies(dp_${appname} TH2Jagged_ext)

  if(NOT "${CMAKE_LINK_FLAGS}" STREQUAL "")
    set_target_properties(dp_${appname} PROPERTIES LINK_FLAGS ${CMAKE_LINK_FLAGS})
  endif()

  if(NOT "${CMAKE_DEPENDLIB_FLAGS}" STREQUAL "")
    target_link_libraries(dp_${appname} ${CMAKE_DEPENDLIB_FLAGS})
  endif()

  install(TARGETS dp_${appname} DESTINATION bin)

endforeach()
