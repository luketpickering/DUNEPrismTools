SET(SRC_FILES)
LIST(APPEND SRC_FILES
  EffCorrector.cxx
  VALORModelClassifier.cxx)

if(DEFINED USE_EIGEN AND USE_EIGEN)
  LIST(APPEND SRC_FILES
    CovarianceHelper.cxx
    VariationBuilders.cxx)
endif()

add_library(DUNE-PRISM_AnaTools
  ${SRC_FILES}
)

if(DEFINED USE_FHICL AND USE_FHICL)
  add_dependencies(DUNE-PRISM_AnaTools fhiclcpp-simple)
endif()

include_directories(${INCDIRS})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# Requried for fhiclcpp
include_directories(${CMAKE_INSTALL_PREFIX}/include)

install(TARGETS DUNE-PRISM_AnaTools DESTINATION lib)

target_link_libraries(DUNE-PRISM_AnaTools ${ROOT_LIBS})
target_link_libraries(DUNE-PRISM_AnaTools gomp)

if(DEFINED USE_EIGEN AND USE_EIGEN)
  add_dependencies(DUNE-PRISM_AnaTools eigen)
endif()

if(NOT "${CMAKE_LINK_FLAGS}" STREQUAL "")
  set_target_properties(DUNE-PRISM_AnaTools PROPERTIES LINK_FLAGS ${CMAKE_LINK_FLAGS})
endif()

if(NOT "${CMAKE_DEPENDLIB_FLAGS}" STREQUAL "")
  target_link_libraries(DUNE-PRISM_AnaTools ${CMAKE_DEPENDLIB_FLAGS})
endif()

SET(MODULETargets ${MODULETargets};DUNE-PRISM_AnaTools PARENT_SCOPE)
