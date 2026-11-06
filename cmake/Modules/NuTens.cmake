if(NOT DEFINED nuTens_BRANCH)
  set(nuTens_BRANCH "v0.7.0")
endif()

if(${UseNuTens} EQUAL 1)
  message(INFO "Adding nutens")
  CPMAddPackage(
    NAME nuTens
    GITHUB_REPOSITORY "ewanwm/nuTens"
    GIT_TAG ${nuTens_BRANCH}
    DOWNLOAD_ONLY
    OPTIONS
    "NT_PROFILING OFF"
  )

  target_compile_definitions(NuOscillatorCompilerOptions INTERFACE UseNuTens=1)
endif()
