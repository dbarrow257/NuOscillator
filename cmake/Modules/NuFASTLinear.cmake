if(NOT DEFINED NuFASTLinear_BRANCH)
  set(NuFASTLinear_BRANCH "v1.1")
endif()

if(${UseNuFASTLinear} EQUAL 1)
  CPMAddPackage(
    NAME NuFAST_LBL
    GIT_TAG ${NuFASTLinear_BRANCH}
    GITHUB_REPOSITORY PeterDenton/NuFast-LBL
    DOWNLOAD_ONLY
  )
  target_compile_definitions(NuOscillatorCompilerOptions INTERFACE UseNuFASTLinear=1)
endif()
