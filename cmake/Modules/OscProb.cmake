if(NOT DEFINED OscProb_BRANCH)
  set(OscProb_BRANCH "v2.3.0")
endif()

if(${UseOscProb} EQUAL 1)
  CPMAddPackage(
    NAME OscProb
    GITHUB_REPOSITORY joaoabcoelho/OscProb
    GIT_TAG ${OscProb_BRANCH}
    DOWNLOAD_ONLY YES
  )

  CPMFindPackage(
    NAME Eigen
    VERSION ${Eigen_BRANCH}
    URL https://gitlab.com/libeigen/eigen/-/archive/${Eigen_BRANCH}/eigen-${Eigen_BRANCH}.tar.gz
    DOWNLOAD_ONLY YES
  )

  if(CMAKE_CXX_STANDARD LESS 17)
    cmessage(WARNING "OscProb requires at least C++17 — changing from ${CMAKE_CXX_STANDARD} to 17")
    set(CMAKE_CXX_STANDARD 17)
  endif()

  if(Eigen_ADDED)
    add_library(Eigen INTERFACE IMPORTED)
    target_include_directories(Eigen INTERFACE ${Eigen_SOURCE_DIR})
  endif()

  execute_process(COMMAND make WORKING_DIRECTORY ${OscProb_SOURCE_DIR})
  set(OSCPROB_LIB_DIR ${OscProb_SOURCE_DIR}/lib)

  add_library(OSCPROB INTERFACE IMPORTED)
  target_include_directories(OSCPROB INTERFACE ${OscProb_SOURCE_DIR})
  target_compile_definitions(NuOscillatorCompilerOptions INTERFACE UseOscProb=1)
endif()
