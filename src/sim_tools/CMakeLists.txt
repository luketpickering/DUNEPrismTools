add_library(DUNE-PRISM_SimTools
  DepoParticle.cxx
  DepoEvent.cxx
  G4ArReader.cxx
  StopDimensions.cxx
)

include_directories(${INCDIRS})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

install(TARGETS DUNE-PRISM_SimTools DESTINATION lib)

target_link_libraries(DUNE-PRISM_SimTools ${ROOT_LIBS})

if(NOT "${CMAKE_LINK_FLAGS}" STREQUAL "")
  set_target_properties(DUNE-PRISM_SimTools PROPERTIES LINK_FLAGS ${CMAKE_LINK_FLAGS})
endif()

if(NOT "${CMAKE_DEPENDLIB_FLAGS}" STREQUAL "")
  target_link_libraries(DUNE-PRISM_SimTools ${CMAKE_DEPENDLIB_FLAGS})
endif()

SET(MODULETargets ${MODULETargets};DUNE-PRISM_SimTools PARENT_SCOPE)
