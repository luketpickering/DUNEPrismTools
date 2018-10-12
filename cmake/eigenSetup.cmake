if(DEFINED USE_EIGEN AND USE_EIGEN)
  ExternalProject_Add(eigen
    PREFIX "${PROJECT_BINARY_DIR}/Ext"
    GIT_REPOSITORY https://github.com/eigenteam/eigen-git-mirror.git
    GIT_TAG 3.3.5
    UPDATE_DISCONNECTED 1
    CONFIGURE_COMMAND ""
    UPDATE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND "")

  LIST(APPEND EXTRA_CXX_FLAGS -DUSE_EIGEN)

  LIST(APPEND INCDIRS ${PROJECT_BINARY_DIR}/Ext/src/eigen/)

  ExternalProject_Add(spectra
    PREFIX "${PROJECT_BINARY_DIR}/Ext"
    GIT_REPOSITORY https://github.com/yixuan/spectra.git
    GIT_TAG v0.6.2
    UPDATE_DISCONNECTED 1
    CONFIGURE_COMMAND ""
    UPDATE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND "")

  LIST(APPEND INCDIRS ${PROJECT_BINARY_DIR}/Ext/src/spectra/include)
endif()
