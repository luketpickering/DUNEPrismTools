add_library(DUNE-PRIMS_FluxTools
  dk2nu_TreeReader.cxx
  OscillationHelper.cxx
  FluxLikelihood.cxx
  FluxFitter.cxx
)

include_directories(${INCLUDE_DIRECTORIES})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

install(TARGETS DUNE-PRIMS_FluxTools DESTINATION lib)

target_link_libraries(DUNE-PRIMS_FluxTools ${ROOT_LIBS})

if(NOT "${CMAKE_LINK_FLAGS}" STREQUAL "")
  set_target_properties(DUNE-PRIMS_FluxTools PROPERTIES LINK_FLAGS ${CMAKE_LINK_FLAGS})
endif()

if(NOT "${CMAKE_DEPENDLIB_FLAGS}" STREQUAL "")
  target_link_libraries(DUNE-PRIMS_FluxTools ${CMAKE_DEPENDLIB_FLAGS})
endif()

SET(MODULETargets ${MODULETargets};DUNE-PRIMS_FluxTools PARENT_SCOPE)
