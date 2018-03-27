LIST(APPEND INCLUDE_DIRECTORIES ${CMAKE_CURRENT_SOURCE_DIR}/utility)
LIST(APPEND INCLUDE_DIRECTORIES ${CMAKE_CURRENT_SOURCE_DIR}/persistency)

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/utility)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/persistency)

#
# add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/flux_tools)
# add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/sim_tools)
# add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/ana_tools)
