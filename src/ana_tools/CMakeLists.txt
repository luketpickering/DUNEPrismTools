add_library(DUNE-PRISM_AnaTools
  EffCorrector.cxx
  VALORModelClassifier.cxx
  CovarianceHelper.cxx
)

include_directories(${INCDIRS})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

install(TARGETS DUNE-PRISM_AnaTools DESTINATION lib)

target_link_libraries(DUNE-PRISM_AnaTools ${ROOT_LIBS})
target_link_libraries(DUNE-PRISM_AnaTools gomp)

if(NOT "${CMAKE_LINK_FLAGS}" STREQUAL "")
  set_target_properties(DUNE-PRISM_AnaTools PROPERTIES LINK_FLAGS ${CMAKE_LINK_FLAGS})
endif()

if(NOT "${CMAKE_DEPENDLIB_FLAGS}" STREQUAL "")
  target_link_libraries(DUNE-PRISM_AnaTools ${CMAKE_DEPENDLIB_FLAGS})
endif()

SET(MODULETargets ${MODULETargets};DUNE-PRISM_AnaTools PARENT_SCOPE)
