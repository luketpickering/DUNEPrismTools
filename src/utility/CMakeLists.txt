add_library(DUNE-PRISM_Utility
  BoundingBox.cxx
  DetectorStop.cxx
  GENIESplineReader.cxx
  ROOTUtility.cxx
  StringParserUtility.cxx
  PhysicsUtility.cxx
)

include_directories(${INCDIRS})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

install(TARGETS DUNE-PRISM_Utility DESTINATION lib)

target_link_libraries(DUNE-PRISM_Utility ${ROOT_LIBS})

if(DEFINED USE_EIGEN AND USE_EIGEN)
  add_dependencies(DUNE-PRISM_Utility eigen)
endif()

if(NOT "${CMAKE_LINK_FLAGS}" STREQUAL "")
  set_target_properties(DUNE-PRISM_Utility PROPERTIES LINK_FLAGS ${CMAKE_LINK_FLAGS})
endif()

if(NOT "${CMAKE_DEPENDLIB_FLAGS}" STREQUAL "")
  target_link_libraries(DUNE-PRISM_Utility ${CMAKE_DEPENDLIB_FLAGS})
endif()

SET(MODULETargets ${MODULETargets};DUNE-PRISM_Utility PARENT_SCOPE)
