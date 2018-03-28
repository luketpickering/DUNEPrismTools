add_library(dp_flux_tools
  dk2nu_TreeReader.cxx
)

include_directories(${INCLUDE_DIRECTORIES})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

install(TARGETS dp_flux_tools DESTINATION lib)

target_link_libraries(dp_flux_tools ${ROOT_LIBS})

if(NOT "${CMAKE_LINK_FLAGS}" STREQUAL "")
  set_target_properties(dp_flux_tools PROPERTIES LINK_FLAGS ${CMAKE_LINK_FLAGS})
endif()

if(NOT "${CMAKE_DEPENDLIB_FLAGS}" STREQUAL "")
  target_link_libraries(dp_flux_tools ${CMAKE_DEPENDLIB_FLAGS})
endif()

SET(MODULETargets ${MODULETargets};dp_flux_tools PARENT_SCOPE)
