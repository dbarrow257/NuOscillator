if(NOT DEFINED CHICLinear_BRANCH)
  set(CHICLinear_BRANCH "v1.0.1")
endif()

if(${UseCHICLinear} EQUAL 1)
  if(CMAKE_CXX_STANDARD LESS 17)
    cmessage(WARNING "CHIC requires at least C++17 — changing from ${CMAKE_CXX_STANDARD} to 17")
    set(CMAKE_CXX_STANDARD 17)
  endif()

  CPMFindPackage(
    NAME Eigen
    VERSION ${Eigen_BRANCH}
    URL https://gitlab.com/libeigen/eigen/-/archive/${Eigen_BRANCH}/eigen-${Eigen_BRANCH}.tar.gz
    DOWNLOAD_ONLY YES
  )

  CPMAddPackage(
    NAME CHIC
    GIT_REPOSITORY https://github.com/pabloferm/CHIC.git
    GIT_TAG ${CHICLinear_BRANCH}
    DOWNLOAD_ONLY YES
  )

  set(CHIC_HEADERS
    ${CHIC_SOURCE_DIR}/src/CHIC.h
    ${CHIC_SOURCE_DIR}/src/opt_constants.h
  )

  add_library(CHIC STATIC
    ${CHIC_SOURCE_DIR}/src/CHIC.cpp
    ${CHIC_HEADERS}
  )

  target_include_directories(CHIC
  PUBLIC
    $<BUILD_INTERFACE:${CHIC_SOURCE_DIR}/src>
    $<INSTALL_INTERFACE:include/CHIC>
  )
  target_include_directories(CHIC PRIVATE ${Eigen_SOURCE_DIR})
  target_compile_options(CHIC PRIVATE -fPIC)
  target_link_libraries(CHIC PRIVATE NuOscillatorCompilerOptions)

  add_library(CHIC::CHIC ALIAS CHIC)
  install(TARGETS CHIC
    EXPORT NuOscillator-targets
    ARCHIVE DESTINATION lib
  )

  target_compile_definitions(NuOscillatorCompilerOptions INTERFACE UseCHICLinear=1)
endif()
