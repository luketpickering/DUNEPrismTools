find_package(Doxygen)

if(DOXYGEN_FOUND)

  file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/dox)

  configure_file(${CMAKE_SOURCE_DIR}/dox/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/dox/Doxyfile @ONLY)
  add_custom_target(doc_generate
  ${DOXYGEN_EXECUTABLE} ${CMAKE_BINARY_DIR}/dox/Doxyfile
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/dox
  COMMENT "Generating documentation with Doxygen... (this will take a while)" VERBATIM
  )

  add_custom_target(docs
  make
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/dox/latex
  COMMENT "Building latex documentation with Doxygen... (this will also take a while)" VERBATIM
  )
  add_dependencies(docs doc_generate)
  install(FILES ${CMAKE_BINARY_DIR}/dox/latex/refman.pdf
    DESTINATION ${CMAKE_BINARY_DIR}/dox
    RENAME DUNE-PRISM_Tools.${DUNEPrismTools_VERSION_STRING}.pdf OPTIONAL)
endif(DOXYGEN_FOUND)
