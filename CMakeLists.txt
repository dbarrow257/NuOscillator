#================================================================================================
#Project details

cmake_minimum_required(VERSION 3.14 FATAL_ERROR)

find_program(CMAKE_C_COMPILER NAMES $ENV{CC} gcc PATHS ENV PATH NO_DEFAULT_PATH)
find_program(CMAKE_CXX_COMPILER NAMES $ENV{CXX} g++ PATHS ENV PATH NO_DEFAULT_PATH)

if(NOT DEFINED CMAKE_CXX_STANDARD)
  set(CMAKE_CXX_STANDARD 11)
endif()

#KS: Store some handy cmake functions
include(${CMAKE_CURRENT_LIST_DIR}/cmake/Modules/NuOscillatorUtils.cmake)

# By default GPU is off
DefineEnabledRequiredSwitch(UseGPU 0)

set(NUOSCILLATOR_VERSION 1.3.0)

if(${UseGPU} EQUAL 0)
  project(NUOSCILLATOR VERSION ${NUOSCILLATOR_VERSION} LANGUAGES CXX)
else()
  project(NUOSCILLATOR VERSION ${NUOSCILLATOR_VERSION} LANGUAGES CXX CUDA)
endif()

LIST(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/Modules)

#Changes default install path to be a subdirectory of the build dir.
#Can set build dir at configure time with -DCMAKE_INSTALL_PREFIX=/install/path
if(CMAKE_INSTALL_PREFIX STREQUAL "" OR CMAKE_INSTALL_PREFIX STREQUAL
    "/usr/local")
  set(CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}/${CMAKE_SYSTEM_NAME}")
elseif(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}/${CMAKE_SYSTEM_NAME}")
endif()

#================================================================================================
#Custom variable sanitisation
add_library(NuOscillatorCompilerOptions INTERFACE)

# download CPM.cmake
file(
  DOWNLOAD
  https://github.com/cpm-cmake/CPM.cmake/releases/download/v0.40.2/CPM.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/cmake/CPM.cmake
)
include(${CMAKE_CURRENT_BINARY_DIR}/cmake/CPM.cmake)

DefineEnabledRequiredSwitch(UseMultithreading 1)
DefineEnabledRequiredSwitch(UseDoubles 1)

DefineEnabledRequiredSwitch(UseCUDAProb3 0)
DefineEnabledRequiredSwitch(UseCUDAProb3Linear 0)
DefineEnabledRequiredSwitch(UseProb3ppLinear 0)
DefineEnabledRequiredSwitch(UseProbGPULinear 0)
DefineEnabledRequiredSwitch(UseNuFASTLinear 0)
DefineEnabledRequiredSwitch(UseNuSQUIDSLinear 0)
DefineEnabledRequiredSwitch(UseOscProb 0)

# KS: NuFastLinear by default
if(UseCUDAProb3        EQUAL 0 AND
    UseCUDAProb3Linear EQUAL 0 AND
    UseProb3ppLinear   EQUAL 0 AND
    UseProbGPULinear   EQUAL 0 AND
    UseNuFASTLinear    EQUAL 0 AND
    UseNuSQUIDSLinear  EQUAL 0 AND
    UseOscProb         EQUAL 0)
  set(UseNuFASTLinear 1)
endif()

LIST(APPEND ALL_Engines
  CUDAProb3
  CUDAProb3Linear
  Prob3ppLinear
  ProbGPULinear
  NuFASTLinear
  NuSQUIDSLinear
  OscProb
)

########################################################
#Default Branches

if(NOT DEFINED CUDAProb3_BRANCH)
  set(CUDAProb3_BRANCH "develop")
endif()

if(NOT DEFINED CUDAProb3Linear_BRANCH)
  set(CUDAProb3Linear_BRANCH "v1.0.1")
endif()

if(NOT DEFINED Prob3ppLinear_BRANCH)
  set(Prob3ppLinear_BRANCH "main")
endif()

if(NOT DEFINED ProbGPULinear_BRANCH)
  set(ProbGPULinear_BRANCH "main")
endif()

if(NOT DEFINED NuFASTLinear_BRANCH)
  set(NuFASTLinear_BRANCH "v1.0")
endif()

if(NOT DEFINED NuSQUIDSLinear_BRANCH)
  set(NuSQUIDSLinear_BRANCH "master")
endif()

if(NOT DEFINED OscProb_BRANCH)
  set(OscProb_BRANCH "v2.1.2")
endif()

########################################################

set(NuOscillator_Engines_Enabled "")
cmessage(STATUS "NuOscillator Features: ")
foreach(f ${ALL_Engines})
  set(engine_var "Use${f}")
  cmessage(STATUS "     ${f}: ${${engine_var}}")
  if(DEFINED ${engine_var} AND ${${engine_var}})
    list(APPEND NuOscillator_Engines_Enabled ${f})
  endif()
endforeach()

cmessage(STATUS "Required variables being used:")
if(${UseGPU} EQUAL 1)
  cmessage(STATUS "\tUsing GPU")
else()
  cmessage(STATUS "\tNot using GPU")
  
  if(${UseProbGPULinear} EQUAL 1)
    cmessage(FATAL_ERROR "ProbGPULinear calculator can only be used when using GPU")
  endif()
endif()

if(${UseMultithreading} EQUAL 1)
  if(DEFINED ENV{OMP_NUM_THREADS} AND NOT "$ENV{OMP_NUM_THREADS}" STREQUAL "")
    cmessage(STATUS "\tUsing Multithreading with nThreads=$ENV{OMP_NUM_THREADS}")
  else()
    execute_process(
      COMMAND lscpu
      COMMAND awk -F: "/^(CPU\\(s\\)|CPU):/ {gsub(/^[ \\t]+/, \"\", \$2); print \$2; exit}"
      OUTPUT_VARIABLE CPU_COUNT
      OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    cmessage(STATUS "\tUsing Multithreading with nThreads=${CPU_COUNT}")
  endif()
else()
  cmessage(STATUS "\tNot using Multithreading")
  if(${UseCUDAProb3} EQUAL 1 OR ${UseCUDAProb3Linear} EQUAL 1)
    cmessage(FATAL_ERROR "CUDAProb3 do not support Multithreading turn off. Please enable it or use different engine")
  endif()
endif()

if(${UseCUDAProb3} EQUAL 1 AND ${UseCUDAProb3Linear} EQUAL 1)
  cmessage(FATAL_ERROR "CUDAProb3 and CUDAProb3Linear and not able to be built at the same time")
endif()

if(${UseMultithreading} EQUAL 1)
  target_compile_options(NuOscillatorCompilerOptions INTERFACE -fopenmp)
  # KS: Not a clue why clang need different...
  if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    target_link_libraries(NuOscillatorCompilerOptions INTERFACE omp)
  else()
    target_link_libraries(NuOscillatorCompilerOptions INTERFACE gomp)
  endif()
endif()

if(NOT DEFINED NuOscillator_Compiler_Flags)
  set(NuOscillator_Compiler_Flags "")
endif()

separate_arguments(NuOscillator_Compiler_Flags)

cmessage(STATUS "Set compiler options: ${NuOscillator_Compiler_Flags}")
add_compile_options(${NuOscillator_Compiler_Flags})

if(${UseGPU} EQUAL 1)
  #KS: Allow user to define CMAKE_CUDA_ARCHITECTURES
  if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
    #KS: See this for more info https://cmake.org/cmake/help/latest/prop_tgt/CUDA_ARCHITECTURES.html
    if( ${CMAKE_VERSION} VERSION_GREATER_EQUAL 3.23 )
      set(CMAKE_CUDA_ARCHITECTURES all )
      #KS: Consider using native, requires cmake 3.24... will be terrible for containers but should results in more optimised code
      #set(CMAKE_CUDA_ARCHITECTURES native )
    else()
      #KS: Apparently with newer cmake and GPU
      set(CMAKE_CUDA_ARCHITECTURES 60 61 70 75 80 86)
    endif()
    #KS: Bit hacky but to properly pass cuda flags we need
    string(REPLACE ";" " " CMAKE_CUDA_ARCHITECTURES_STRING "${CMAKE_CUDA_ARCHITECTURES}")
  else()
    #KS this may look hacky however CPM isn't build for passing stuff like this. If CMAKE_CUDA_ARCHITECTURES is passed CPM it will be string not list. Thus we convert it to list
    set(CMAKE_CUDA_ARCHITECTURES_STRING ${CMAKE_CUDA_ARCHITECTURES})
    string(REPLACE " " ";" CMAKE_CUDA_ARCHITECTURES "${CMAKE_CUDA_ARCHITECTURES}")
  endif()
  
  cmessage(STATUS "Set CUDA ARCHITECTURES options: ${CMAKE_CUDA_ARCHITECTURES_STRING}")
endif()

#================================================================================================
#Setup Deps
set(YAML_CPP_VERSION 0.8.0)
CPMAddPackage(
    NAME yaml-cpp
    VERSION ${YAML_CPP_VERSION}
    GITHUB_REPOSITORY "jbeder/yaml-cpp"
    GIT_TAG "${YAML_CPP_VERSION}"
    OPTIONS
      "YAML_CPP_INSTALL ON"
      "YAML_CPP_BUILD_TESTS OFF"
      "YAML_CPP_BUILD_CONTRIB OFF"
      "YAML_BUILD_SHARED_LIBS ON"
)

if(NOT TARGET yaml-cpp::yaml-cpp)
  cmessage(FATAL_ERROR "NuOscillator Expected dependency target: yaml-cpp::yaml-cpp")
endif()

find_package(ROOT 6.18 REQUIRED)

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

if(${UseGPU} EQUAL 1)
  if(${UseProbGPULinear} EQUAL 1)
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
  endif()
endif()

if(${UseProb3ppLinear} EQUAL 1)
  CPMAddPackage(
    NAME Prob3plusplus
    VERSION 3.10.4
    GITHUB_REPOSITORY rogerwendell/Prob3plusplus
    GIT_TAG ${Prob3ppLinear_BRANCH}
  )

  if(NOT TARGET Prob3plusplus)
    cmessage(FATAL_ERROR "Could not find Prob3plusplus")
  endif()

  # KS: Add additional compilation flags to be used for Prob3++
  target_compile_options(Prob3plusplus PRIVATE ${NuOscillator_Compiler_Flags})
  target_compile_definitions(NuOscillatorCompilerOptions INTERFACE UseProb3ppLinear=1)

  install(TARGETS Prob3plusplus
    EXPORT NuOscillator-targets
    LIBRARY DESTINATION lib/)
endif()

if(${UseNuFASTLinear} EQUAL 1)
  CPMAddPackage(
    NAME NuFAST
    GIT_TAG ${NuFASTLinear_BRANCH}
    GITHUB_REPOSITORY PeterDenton/NuFast
    DOWNLOAD_ONLY
  )
  target_compile_definitions(NuOscillatorCompilerOptions INTERFACE UseNuFASTLinear=1)
endif()

if(${UseOscProb} EQUAL 1)
  CPMAddPackage(
    NAME ${OscProb_BRANCH}
    GITHUB_REPOSITORY joaoabcoelho/OscProb
    GIT_TAG v2.1.2
  )

  CPMAddPackage(
    NAME Eigen
    VERSION 3.4.0
    URL https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz
    DOWNLOAD_ONLY YES
  )
  if(Eigen_ADDED)
    add_library(Eigen INTERFACE IMPORTED)
    target_include_directories(Eigen INTERFACE ${Eigen_SOURCE_DIR})
  endif()

  execute_process(COMMAND make WORKING_DIRECTORY ${OscProb_SOURCE_DIR})

  add_library(OSCPROB INTERFACE IMPORTED)
  target_include_directories(OSCPROB INTERFACE  ${OscProb_SOURCE_DIR})

  target_compile_definitions(NuOscillatorCompilerOptions INTERFACE UseOscProb=1)

  set(OSCPROB_LIB_DIR ${OscProb_SOURCE_DIR}/lib)
endif()

if(${UseNuSQUIDSLinear} EQUAL 1)
  find_package(GSL REQUIRED)

  find_package(HDF5 QUIET COMPONENTS CXX)
  if(NOT ${HDF5_FOUND})
    find_package(HDF5 QUIET REQUIRED)
  endif()

  include(ExternalProject)

  ExternalProject_Add(SQUIDS
    GIT_REPOSITORY https://github.com/jsalvado/SQuIDS.git
    GIT_TAG master
    BUILD_IN_SOURCE = 1
    CONFIGURE_HANDLED_BY_BUILD false
    CONFIGURE_COMMAND <SOURCE_DIR>/configure --with-gsl=${GSL_ROOT_DIR} --prefix=<INSTALL_DIR>
    BUILD_COMMAND make
    INSTALL_COMMAND make install
  )

  SET(SQUIDS_INC_DIR ${CMAKE_CURRENT_BINARY_DIR}/SQUIDS-prefix/include)
  SET(SQUIDS_LIB_DIR ${CMAKE_CURRENT_BINARY_DIR}/SQUIDS-prefix/lib)

  add_library(SQUIDS::SQUIDS SHARED IMPORTED)
  set_target_properties(SQUIDS::SQUIDS PROPERTIES IMPORTED_LOCATION ${SQUIDS_LIB_DIR}/libSQuIDS.so)

  ExternalProject_Add(NUSQUIDS
    GIT_REPOSITORY https://github.com/arguelles/nuSQuIDS.git
    GIT_TAG ${NuSQUIDSLinear_BRANCH}
    BUILD_IN_SOURCE = 1
    CONFIGURE_HANDLED_BY_BUILD false
    CONFIGURE_COMMAND <SOURCE_DIR>/configure --with-gsl=${GSL_ROOT_DIR} --with-hdf5=${HDF5_INCLUDE_DIRS}/../ --with-squids=${CMAKE_CURRENT_BINARY_DIR}/SQUIDS-prefix --prefix=<INSTALL_DIR>
    BUILD_COMMAND make examples
    INSTALL_COMMAND make install
  )

  SET(NUSQUIDS_SRC_DIR ${CMAKE_CURRENT_BINARY_DIR}/NUSQUIDS-prefix/src/NUSQUIDS/)
  SET(NUSQUIDS_INC_DIR ${CMAKE_CURRENT_BINARY_DIR}/NUSQUIDS-prefix/include/)
  SET(NUSQUIDS_LIB_DIR ${CMAKE_CURRENT_BINARY_DIR}/NUSQUIDS-prefix/lib/)

  add_dependencies(NUSQUIDS SQUIDS)
  add_library(NUSQUIDS::NUSQUIDS SHARED IMPORTED)
  set_target_properties(NUSQUIDS::NUSQUIDS PROPERTIES IMPORTED_LOCATION ${NUSQUIDS_LIB_DIR}/libnuSQuIDS.so)

  target_compile_definitions(NuOscillatorCompilerOptions INTERFACE UseNuSQUIDSLinear=1)
endif()

if(${UseGPU} EQUAL 1)
  target_compile_definitions(NuOscillatorCompilerOptions INTERFACE UseGPU=1)
else()
  target_compile_definitions(NuOscillatorCompilerOptions INTERFACE UseGPU=0)
endif()

if(${UseMultithreading} EQUAL 1)
  target_compile_definitions(NuOscillatorCompilerOptions INTERFACE UseMultithreading=1)
else()
  target_compile_definitions(NuOscillatorCompilerOptions INTERFACE UseMultithreading=0)
endif()

if(${UseDoubles} EQUAL 1)
  target_compile_definitions(NuOscillatorCompilerOptions INTERFACE UseDoubles=1)
else()
  target_compile_definitions(NuOscillatorCompilerOptions INTERFACE UseDoubles=0)
endif()

#================================================================================================
#Subdirectories to add

add_subdirectory(OscProbCalcer)
add_subdirectory(Oscillator)
add_subdirectory(Apps)

#================================================================================================

set(LIBDEST ${CMAKE_INSTALL_PREFIX}/lib)

configure_file(cmake/Templates/setup.NuOscillator.sh.in "${PROJECT_BINARY_DIR}/setup.NuOscillator.sh" @ONLY)
install(FILES "${PROJECT_BINARY_DIR}/setup.NuOscillator.sh" DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)

install(DIRECTORY NuOscillatorConfigs DESTINATION ${CMAKE_INSTALL_PREFIX}/)
install(DIRECTORY Inputs DESTINATION ${CMAKE_INSTALL_PREFIX}/)
install(DIRECTORY Constants DESTINATION ${CMAKE_INSTALL_PREFIX}/)

set_target_properties(NuOscillatorCompilerOptions PROPERTIES EXPORT_NAME CompilerOptions)
install(TARGETS NuOscillatorCompilerOptions
  EXPORT NuOscillator-targets
  LIBRARY DESTINATION lib/)

#================================================================================================

#This is to export the target properties of NuOscillator
#Anything that links to "NuOscillator" will get all of these target properties
add_library(NuOscillator INTERFACE)
target_link_libraries(NuOscillator INTERFACE Oscillator OscProbCalcer NuOscillatorCompilerOptions)
set_target_properties(NuOscillator PROPERTIES EXPORT_NAME All)

install(TARGETS NuOscillator
  EXPORT NuOscillator-targets
  LIBRARY DESTINATION lib/)

install(EXPORT NuOscillator-targets
  FILE NuOscillatorTargets.cmake
  NAMESPACE NuOscillator::
  DESTINATION ${CMAKE_INSTALL_PREFIX}/
)

install(DIRECTORY
  NuOscillatorConfigs
  DESTINATION ${CMAKE_BINARY_DIR}/NuOscillatorConfigs)

add_library(NuOscillator::All ALIAS NuOscillator)


############################  Install  ####################################
include(CMakePackageConfigHelpers)
configure_package_config_file(
  ${CMAKE_CURRENT_LIST_DIR}/cmake/Templates/NuOscillatorConfig.cmake.in ${CMAKE_BINARY_DIR}/NuOscillatorConfig.cmake
  INSTALL_DESTINATION
  /this/is/ignored/for/some/reason/thanks/kitware
  NO_SET_AND_CHECK_MACRO
  NO_CHECK_REQUIRED_COMPONENTS_MACRO)

write_basic_package_version_file(${CMAKE_BINARY_DIR}/NuOscillatorConfigVersion.cmake
  VERSION ${NUOSCILLATOR_VERSION}
  COMPATIBILITY AnyNewerVersion)

install(FILES
  ${CMAKE_BINARY_DIR}/NuOscillatorConfig.cmake
  ${CMAKE_BINARY_DIR}/NuOscillatorConfigVersion.cmake
  DESTINATION ${CMAKE_INSTALL_PREFIX}/lib/cmake/NuOscillator)

include(nuoscillator-config)
