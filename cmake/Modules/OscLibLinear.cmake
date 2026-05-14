if(NOT DEFINED OscLib_BRANCH)
  set(OscLib_BRANCH "v1.2.1")
endif()

if(${UseOscLibLinear} EQUAL 1)
  # Check if C++ standard is less than 17
  if(CMAKE_CXX_STANDARD LESS 17)
    cmessage(FATAL_ERROR "OscLib requires at least C++17, but you are using C++${CMAKE_CXX_STANDARD}.")
  endif()

  if(${UseDoubles} EQUAL 0)
    cmessage(FATAL_ERROR "OscLib currently don't support float, please switch to double or change engine. Gomenasai")
  endif()

  find_package(Boost REQUIRED)
  find_package(GSL REQUIRED)

  CPMFindPackage(
    NAME Eigen
    VERSION ${Eigen_BRANCH}
    URL https://gitlab.com/libeigen/eigen/-/archive/${Eigen_BRANCH}/eigen-${Eigen_BRANCH}.tar.gz
    DOWNLOAD_ONLY YES
  )
  if(Eigen_ADDED)
    add_library(Eigen INTERFACE IMPORTED)
    target_include_directories(Eigen INTERFACE ${Eigen_SOURCE_DIR})
  endif()

  CPMAddPackage(
    NAME OscLib
    GITHUB_REPOSITORY cafana/OscLib
    GIT_TAG ${OscLib_BRANCH}
    DOWNLOAD_ONLY YES
  )
  # KS: Collect OscLib sources and header files
  file(GLOB OSCLIB_SOURCES
    ${OscLib_SOURCE_DIR}/OscLib/*.cxx
  )
  file(GLOB OSCLIB_HEADERS
    ${OscLib_SOURCE_DIR}/OscLib/*.h
  )
  add_library(OscLibLinear SHARED
    ${OSCLIB_SOURCES}
  )

  set_target_properties(OscLibLinear PROPERTIES
    PUBLIC_HEADER "${OSCLIB_HEADERS}"
    EXPORT_NAME OscLib
  )

  target_link_libraries(OscLibLinear PRIVATE ROOT::Core Eigen GSL::gsl)
  # Properly scoped include dirs
  target_include_directories(OscLibLinear PUBLIC
      $<BUILD_INTERFACE:${OscLib_SOURCE_DIR}>
      $<INSTALL_INTERFACE:include>
  )
  install(TARGETS OscLibLinear
    EXPORT NuOscillator-targets
    LIBRARY DESTINATION lib/
    PUBLIC_HEADER DESTINATION include/
  )
  target_compile_definitions(NuOscillatorCompilerOptions INTERFACE UseOscLibLinear=1)
endif()

if(OSCLIB_LIB_DIR)
    file(RELATIVE_PATH OSCLIB_REL_LIB_PATH
       "${CMAKE_INSTALL_PREFIX}"
       "${OSCLIB_LIB_DIR}")
endif()
