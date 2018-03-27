add_library(dp_persistency
  CondensedDepositsTreeReader.cxx
  DepositsSummaryTreeReader.cxx
  FluxFitResultsTreeReader.cxx
  OscillationParametersTreeReader.cxx
  SimConfigTreeReader.cxx
  SliceConfigTreeReader.cxx
  StopConfigTreeReader.cxx
)

include_directories(${INCLUDE_DIRECTORIES})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

install(TARGETS dp_persistency DESTINATION lib)

target_link_libraries(dp_persistency ${MODULETargets})
target_link_libraries(dp_persistency ${ROOT_LIBS})

if(NOT "${CMAKE_LINK_FLAGS}" STREQUAL "")
  set_target_properties(dp_persistency PROPERTIES LINK_FLAGS ${CMAKE_LINK_FLAGS})
endif()

if(NOT "${CMAKE_DEPENDLIB_FLAGS}" STREQUAL "")
  target_link_libraries(dp_persistency ${CMAKE_DEPENDLIB_FLAGS})
endif()
