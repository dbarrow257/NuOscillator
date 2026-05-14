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

set(NUOSCILLATOR_VERSION 2.0.0)

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

DefineEnabledRequiredSwitch(UseMultithreading 1)
DefineEnabledRequiredSwitch(UseDoubles 1)

LIST(APPEND ALL_Engines
  CUDAProb3
  CUDAProb3Linear
  Prob3ppLinear
  ProbGPULinear
  NuFASTLinear
  NuFASTEarth
  NuSQUIDSLinear
  OscProb
  OscLibLinear
  GLoBESLinear
  CHICLinear
)
# KS: Tells whether all oscillators were disabled
set(ALL_OSCILLATORS_DISABLED TRUE)

foreach(option IN LISTS ALL_Engines)
  DefineEnabledRequiredSwitch(Use${option} 0)
  # If at least one oscillator is enabled then set ALL_OSCILLATORS_DISABLED to false
  if(Use${option} EQUAL 1)
    set(ALL_OSCILLATORS_DISABLED FALSE)
  endif()
endforeach()

# KS: NuFastLinear by default
if(ALL_OSCILLATORS_DISABLED)
  set(UseNuFASTLinear 1)
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

if(${UseMultithreading} EQUAL 1)
  if(CMAKE_CXX_COMPILER_ID STREQUAL "IntelLLVM")
    #ETA: intel icpx requires this I believe
    find_package(OpenMP REQUIRED)
    target_compile_options(NuOscillatorCompilerOptions INTERFACE -qopenmp)
    target_link_libraries(NuOscillatorCompilerOptions INTERFACE OpenMP::OpenMP_CXX)
  # KS: Not a clue why clang need different...
  elseif(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    target_compile_options(NuOscillatorCompilerOptions INTERFACE -fopenmp)
    target_link_libraries(NuOscillatorCompilerOptions INTERFACE omp)
  else()
    target_compile_options(NuOscillatorCompilerOptions INTERFACE -fopenmp)
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
########################################################
# Load all dependencies
include(${CMAKE_CURRENT_LIST_DIR}/cmake/Modules/NuOscillatorDependencies.cmake)
foreach(engine IN LISTS ALL_Engines)
  include(${CMAKE_CURRENT_LIST_DIR}/cmake/Modules/${engine}.cmake)
endforeach()
########################################################
#Subdirectories to add

add_subdirectory(OscProbCalcer)
add_subdirectory(Oscillator)
add_subdirectory(Apps)

#================================================================================================
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
install(DIRECTORY cmake DESTINATION ${CMAKE_BINARY_DIR})

include(nuoscillator-config)
