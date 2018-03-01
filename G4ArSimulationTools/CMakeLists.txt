SET(APPS DunePrismCondenser FullDetTreeStopProcessor
  FluxVarFriendTreeBuilder XSecVarFriendTreeBuilder G4ArEventDisplay)

include_directories(${INCLUDE_DIRECTORIES})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_library(G4ArReader G4ArReader.cxx)
install(TARGETS G4ArReader DESTINATION lib)

foreach(a ${APPS})

  add_executable(dp_${a} ${a}.cxx)

  target_link_libraries(dp_${a} ${MODULETargets})
  target_link_libraries(dp_${a} ${ROOT_LIBS})
  target_link_libraries(dp_${a} G4ArReader)


  if(NOT "${CMAKE_LINK_FLAGS}" STREQUAL "")
    set_target_properties(dp_${a} PROPERTIES LINK_FLAGS ${CMAKE_LINK_FLAGS})
  endif()

  if(NOT "${CMAKE_DEPENDLIB_FLAGS}" STREQUAL "")
    target_link_libraries(dp_${a} ${CMAKE_DEPENDLIB_FLAGS})
  endif()

  install(TARGETS dp_${a} DESTINATION bin)


endforeach()
