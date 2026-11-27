if(NOT DEFINED nuOscillationLinear_BRANCH)
  set(nuOscillationLinear_BRANCH "54f55c10f2e4c03d6589b7bede4afe62e7523583")
endif()

if(${UsenuOscillationLinear} EQUAL 1)
  CPMAddPackage(
    NAME vOscillation
    GIT_REPOSITORY https://github.com/railroad2/vOscillation.git
    GIT_TAG ${nuOscillationLinear_BRANCH}
    DOWNLOAD_ONLY YES
  )

  # Terrible hack to fix compilation
  set(NUOSCILLATION_HEADER "${vOscillation_SOURCE_DIR}/header/vOscillation.hh")
  set(NUOSCILLATION_SOURCE "${vOscillation_SOURCE_DIR}/src/vOscillation.cc")

  add_library(vOscillation
    ${NUOSCILLATION_SOURCE}
    ${NUOSCILLATION_HEADERS}
  )
  set_target_properties(vOscillation PROPERTIES
    POSITION_INDEPENDENT_CODE ON
    PUBLIC_HEADER "${NUOSCILLATION_HEADERS}"
  )
  target_link_libraries(vOscillation PUBLIC ROOT::Core
                                     PRIVATE NuOscillatorCompilerOptions)
  target_include_directories(vOscillation PUBLIC
    $<BUILD_INTERFACE:${vOscillation_SOURCE_DIR}/>
    $<INSTALL_INTERFACE:include/vOscillation>
  )
  install(TARGETS vOscillation
    EXPORT NuOscillator-targets
    ARCHIVE DESTINATION lib
    PUBLIC_HEADER DESTINATION include
  )
  target_compile_definitions(NuOscillatorCompilerOptions INTERFACE UsenuOscillationLinear=1)
endif()
