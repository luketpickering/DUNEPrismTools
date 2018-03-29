add_library(DUNE-PRISM_Persistency
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

install(TARGETS DUNE-PRISM_Persistency DESTINATION lib)

target_link_libraries(DUNE-PRISM_Persistency ${ROOT_LIBS})

if(NOT "${CMAKE_LINK_FLAGS}" STREQUAL "")
  set_target_properties(DUNE-PRISM_Persistency PROPERTIES LINK_FLAGS ${CMAKE_LINK_FLAGS})
endif()

if(NOT "${CMAKE_DEPENDLIB_FLAGS}" STREQUAL "")
  target_link_libraries(DUNE-PRISM_Persistency ${CMAKE_DEPENDLIB_FLAGS})
endif()

SET(MODULETargets ${MODULETargets};DUNE-PRISM_Persistency PARENT_SCOPE)
