if(NOT DEFINED CUDAProb3_BRANCH)
  set(CUDAProb3_BRANCH "v1.2.1")
endif()

if(${UseCUDAProb3} EQUAL 1)
  CPMFindPackage(
    NAME CUDAProb3
    GITHUB_REPOSITORY "dbarrow257/CUDAProb3"
    GIT_TAG ${CUDAProb3_BRANCH}
    DOWNLOAD_ONLY
  )
  if(NOT CUDAProb3_ADDED)
    cmessage(FATAL_ERROR "Could not download CUDAProb3")
  endif()

  target_compile_definitions(NuOscillatorCompilerOptions INTERFACE UseCUDAProb3=1)
endif()

# KS: This extract relative path to CudaProb and others. This allow in setup to use ${NUOSCILLATOR_ROOT}/RELATIVE_PART
# this way even if we move directory we can still use NuOscillator
if(CUDAProb3_SOURCE_DIR)
  file(RELATIVE_PATH CUDAProb3_RELATIVE_PATH
       "${CMAKE_INSTALL_PREFIX}"
       "${CUDAProb3_SOURCE_DIR}")
endif()
