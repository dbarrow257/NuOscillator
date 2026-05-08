if(NOT DEFINED ProbGPULinear_BRANCH)
  set(ProbGPULinear_BRANCH "v1.0.0")
endif()

if(${UseProbGPULinear} EQUAL 1)
  if(${UseGPU} EQUAL 1)

    CPMFindPackage(
      NAME ProbGPU
      GIT_TAG ${ProbGPULinear_BRANCH}
      GITHUB_REPOSITORY dbarrow257/ProbGPU
      OPTIONS
      "CMAKE_CUDA_ARCHITECTURES ${CMAKE_CUDA_ARCHITECTURES_STRING}"
    )
    if(NOT TARGET ProbGPU)
      cmessage(FATAL_ERROR "Could not find ProbGPU")
    endif()

    target_compile_definitions(NuOscillatorCompilerOptions INTERFACE UseProbGPULinear=1)
    install(TARGETS ProbGPU
      EXPORT NuOscillator-targets
      LIBRARY DESTINATION lib/)
  else()
    cmessage(FATAL_ERROR "Enabled ProbGPU, however UseGPU is off. Disable ProbGPU or enable UseGPU option")
  endif()
endif()
