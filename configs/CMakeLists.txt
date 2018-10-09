file(GLOB FLXWINXML ${CMAKE_CURRENT_SOURCE_DIR}/flux_windows/*.xml)
install(FILES ${FLXWINXML} DESTINATION ${CMAKE_INSTALL_PREFIX}/configs/flux_windows)

file(GLOB GDMLXML ${CMAKE_CURRENT_SOURCE_DIR}/gdml/*.gdml)
install(FILES ${GDMLXML} DESTINATION ${CMAKE_INSTALL_PREFIX}/configs/gdml)

file(GLOB GNSPLRDR ${CMAKE_CURRENT_SOURCE_DIR}/genie_spline_reader/*.xml)
install(FILES ${GNSPLRDR} DESTINATION ${CMAKE_INSTALL_PREFIX}/configs/genie_spline_reader)

file(GLOB RPXML ${CMAKE_CURRENT_SOURCE_DIR}/run_plans/*.xml)
install(FILES ${RPXML} DESTINATION ${CMAKE_INSTALL_PREFIX}/configs/run_plans)
