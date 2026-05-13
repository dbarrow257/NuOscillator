if(NOT DEFINED NuFASTEarth_BRANCH)
  set(NuFASTEarth_BRANCH "v1.2.0")
endif()

if(${UseNuFASTEarth} EQUAL 1)
  CPMFindPackage(
    NAME NUFAST_EARTH
    GIT_TAG ${NuFASTEarth_BRANCH}
    GITHUB_REPOSITORY PeterDenton/NuFast-Earth
    DOWNLOAD_ONLY YES
  )

  file(GLOB NUFAST_EARTH_SOURCES "${NUFAST_EARTH_SOURCE_DIR}/src/*.cpp")
  add_library(nufast_earth ${NUFAST_EARTH_SOURCES})
  set_target_properties(nufast_earth PROPERTIES POSITION_INDEPENDENT_CODE ON)

  target_include_directories(nufast_earth PUBLIC
    $<BUILD_INTERFACE:${NUFAST_EARTH_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
  )

  install(TARGETS nufast_earth
    EXPORT NuOscillator-targets
    LIBRARY DESTINATION lib/
  )

  target_compile_definitions(NuOscillatorCompilerOptions INTERFACE UseNuFASTEarth=1)
endif()
