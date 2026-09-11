if(NOT DEFINED Prob3pp_BRANCH)
  set(Prob3pp_BRANCH "v3r20")
endif()

if(${UseProb3pp} EQUAL 1)
  CPMAddPackage(
    NAME Prob3plusplus
    GITHUB_REPOSITORY rogerwendell/Prob3plusplus
    GIT_TAG ${Prob3pp_BRANCH}
  )

  if(NOT TARGET Prob3plusplus)
    cmessage(FATAL_ERROR "Could not find Prob3plusplus")
  endif()

  # KS: Add additional compilation flags to be used for Prob3++
  target_compile_options(Prob3plusplus PRIVATE ${NuOscillator_Compiler_Flags})
  target_compile_definitions(NuOscillatorCompilerOptions INTERFACE UseProb3pp=1)

  install(TARGETS Prob3plusplus
    EXPORT NuOscillator-targets
    LIBRARY DESTINATION lib/)
endif()
