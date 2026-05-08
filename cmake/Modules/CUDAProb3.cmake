if(NOT DEFINED CUDAProb3_BRANCH)
  set(CUDAProb3_BRANCH "v1.2.1")
endif()

if(NOT DEFINED CUDAProb3Linear_BRANCH)
  set(CUDAProb3Linear_BRANCH "v1.1.0")
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

if(${UseCUDAProb3Linear} EQUAL 1)
  if(${UseGPU} EQUAL 1)
    CPMFindPackage(
      NAME CUDAProb3Linear
      GITHUB_REPOSITORY "mach3-software/CUDAProb3"
      GIT_TAG ${CUDAProb3Linear_BRANCH}
      OPTIONS
      "GPU_ON ON"
      "CMAKE_CUDA_ARCHITECTURES ${CMAKE_CUDA_ARCHITECTURES_STRING}"
    )
    # KS: Terrible hack because CUDAProb3Beam is mess :(
    install(TARGETS CUDAProb3Beam
      EXPORT NuOscillator-targets
      LIBRARY DESTINATION lib/
      )

   # KS: This is bit hacky but CUDAProb3Linear is a mess. Whole functions are written in header file
   # and someone wrote them to be hidden under GPU_ON and wrote shitty cmake which do not propagate this definition.
   # Therefore NuOscillator is fixing other people mess
   target_compile_definitions(NuOscillatorCompilerOptions INTERFACE GPU_ON)
  else()
    CPMFindPackage(
      NAME CUDAProb3Linear
      GITHUB_REPOSITORY "mach3-software/CUDAProb3"
      GIT_TAG ${CUDAProb3Linear_BRANCH}
      OPTIONS
      "CPU_ONLY ON"
      "GPU_ON OFF"
    )
    # KS: Terrible hack because CUDAProb3Beam is mess :(
    install(TARGETS CUDAProb3Beam
      EXPORT NuOscillator-targets
      LIBRARY DESTINATION lib/
      PUBLIC_HEADER DESTINATION include/CUDAProb3Beam)
  endif()
  if(NOT TARGET CUDAProb3Beam OR NOT TARGET CUDAProb3Atmos)
    cmessage(FATAL_ERROR "Could not find CUDAProb3Linear")
  endif()
  target_compile_definitions(NuOscillatorCompilerOptions INTERFACE UseCUDAProb3Linear=1)
endif()
