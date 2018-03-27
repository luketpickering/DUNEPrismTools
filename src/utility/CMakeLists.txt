add_library(dp_utility
  DetectorStop.cxx
  GENIESplineReader.cxx
  OscillationHelper.cxx
  ROOTUtility.cxx
  StringParserUtility.cxx
  PhysicsUtility.cxx
)

include_directories(${INCLUDE_DIRECTORIES})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

install(TARGETS dp_utility DESTINATION lib)

target_link_libraries(dp_utility ${MODULETargets})
target_link_libraries(dp_utility ${ROOT_LIBS})

if(NOT "${CMAKE_LINK_FLAGS}" STREQUAL "")
  set_target_properties(dp_utility PROPERTIES LINK_FLAGS ${CMAKE_LINK_FLAGS})
endif()

if(NOT "${CMAKE_DEPENDLIB_FLAGS}" STREQUAL "")
  target_link_libraries(dp_utility ${CMAKE_DEPENDLIB_FLAGS})
endif()
