#================================================================================================
#Project details

cmake_minimum_required(VERSION 3.14 FATAL_ERROR)

find_program(CMAKE_C_COMPILER NAMES $ENV{CC} gcc PATHS ENV PATH NO_DEFAULT_PATH)
find_program(CMAKE_CXX_COMPILER NAMES $ENV{CXX} g++ PATHS ENV PATH NO_DEFAULT_PATH)

if(NOT DEFINED CMAKE_CXX_STANDARD)
  set(CMAKE_CXX_STANDARD 11)
endif()

if(NOT DEFINED UseGPU)
    message(FATAL_ERROR "UseGPU is not defined, add -DUseGPU=<1/0>")
endif()

set(NUOSCILLATOR_VERSION 0.0)

if(${UseGPU} EQUAL 0)
    project(NUOSCILLATOR VERSION ${NUOSCILLATOR_VERSION} LANGUAGES CXX)
else()
    project(NUOSCILLATOR VERSION ${NUOSCILLATOR_VERSION} LANGUAGES CXX CUDA)
endif()

LIST(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/Modules)

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

message(STATUS "NuOscillator Features: ")
foreach(VAR UseGPU UseMultithreading UseDoubles UseCUDAProb3 UseCUDAProb3Linear UseProb3ppLinear UseProbGPULinear UseNuFASTLinear UseNuSQUIDSLinear)
  if(NOT DEFINED ${VAR})
    message(FATAL_ERROR "${VAR} not defined")
  endif()
  message(STATUS "\t${VAR}: ${${VAR}}")
endforeach(VAR)

message(STATUS "Required variables being used:")

if(${UseGPU} EQUAL 1)
	message(STATUS "\tUsing GPU")
else()
	message(STATUS "\tNot using GPU")

	if(${UseProbGPULinear} EQUAL 1)
		message(FATAL_ERROR "ProbGPULinear calculator can only be used when using GPU")
	endif()
endif()

if(${UseMultithreading} EQUAL 1)
    message(STATUS "\tUsing Multithreading with nThreads=$ENV{OMP_NUM_THREADS}")
else()
    message(STATUS "\tNot using Multithreading")
endif()

if(${UseDoubles} EQUAL 1)
    message(STATUS "\tUsing FLOAT_T=double")
else()
    message(STATUS "\tUsing FLOAT_T=float")
endif()

if(${UseCUDAProb3} EQUAL 1 AND ${UseCUDAProb3Linear} EQUAL 1)
    message(FATAL_ERROR "CUDAProb3 and CUDAProb3Linear and not able to be built at the same time")
endif()

if(${UseMultithreading} EQUAL 1)
  target_compile_options(NuOscillatorCompilerOptions INTERFACE -fopenmp)
  target_link_libraries(NuOscillatorCompilerOptions INTERFACE gomp)
endif()

if(NOT DEFINED NuOscillator_Compiler_Flags)
    set(NuOscillator_Compiler_Flags "")
endif()

separate_arguments(NuOscillator_Compiler_Flags)

message(STATUS "Set compiler options: ${NuOscillator_Compiler_Flags}")
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
  
  message(STATUS "Set CUDA ARCHITECTURES options: ${CMAKE_CUDA_ARCHITECTURES_STRING}")
endif()

#================================================================================================
#Setup Deps

# download CPM.cmake
file(
  DOWNLOAD
  https://github.com/cpm-cmake/CPM.cmake/releases/download/v0.38.8/CPM.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/cmake/CPM.cmake
)
include(${CMAKE_CURRENT_BINARY_DIR}/cmake/CPM.cmake)

CPMFindPackage(
    NAME yaml-cpp
    VERSION 0.7.0
    GITHUB_REPOSITORY "jbeder/yaml-cpp"
    GIT_TAG "yaml-cpp-0.7.0"
    OPTIONS
      "YAML_BUILD_SHARED_LIBS ON"
)

CPMFindPackage(
      NAME CMakeModules
      GIT_TAG stable
      GITHUB_REPOSITORY NuHepMC/CMakeModules
      DOWNLOAD_ONLY
  )
include(${CMakeModules_SOURCE_DIR}/NuHepMCModules.cmake)

include(ROOT)
if(NOT TARGET ROOT::ROOT)
  cmessage(FATAL_ERROR "Expected dependency target: ROOT::ROOT")
endif()
if(DEFINED ROOT_CXX_STANDARD AND ROOT_CXX_STANDARD GREATER CMAKE_CXX_STANDARD)
  set(CMAKE_CXX_STANDARD ${ROOT_CXX_STANDARD})
endif()

if(${UseCUDAProb3} EQUAL 1)
	   CPMFindPackage(
		NAME CUDAProb3
		GITHUB_REPOSITORY "dbarrow257/CUDAProb3"
		GIT_TAG "develop"
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
			NAME CUDAProb3
        		GITHUB_REPOSITORY "mach3-software/CUDAProb3"
        		GIT_TAG "main"
			OPTIONS
			"GPU_ON ON"
			"CMAKE_CUDA_ARCHITECTURES ${CMAKE_CUDA_ARCHITECTURES_STRING}"
    			)
    else()
    CPMFindPackage(
        NAME CUDAProb3
        GITHUB_REPOSITORY "mach3-software/CUDAProb3"
        GIT_TAG "main"
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
			GIT_TAG main
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
        GIT_TAG main
    )

    if(NOT TARGET Prob3plusplus)
        cmessage(FATAL_ERROR "Could not find Prob3plusplus")
    endif()
    target_compile_definitions(NuOscillatorCompilerOptions INTERFACE UseProb3ppLinear=1)

    install(TARGETS Prob3plusplus
        EXPORT NuOscillator-targets
        LIBRARY DESTINATION lib/)
endif()

if(${UseNuFASTLinear} EQUAL 1)
    CPMAddPackage(
      	NAME NuFAST
      	GIT_TAG main
      	GITHUB_REPOSITORY PeterDenton/NuFast
      	DOWNLOAD_ONLY
      )
      target_compile_definitions(NuOscillatorCompilerOptions INTERFACE UseNuFASTLinear=1)
endif()

if(${UseNuSQUIDSLinear} EQUAL 1)
  find_package(GSL REQUIRED)
  find_package(HDF5 REQUIRED)

  include(ExternalProject)

  ExternalProject_Add(SQUIDS
    GIT_REPOSITORY https://github.com/jsalvado/SQuIDS.git
    GIT_TAG master
    BINARY_DIR = <SOURCE_DIR> 
    CONFIGURE_HANDLED_BY_BUILD true
    CONFIGURE_COMMAND <SOURCE_DIR>/configure --with-gsl=${GSL_ROOT_DIR} --prefix=<INSTALL_DIR>
    BUILD_COMMAND make
    INSTALL_COMMAND make install
  )
  
  ExternalProject_Add(NUSQUIDS
    GIT_REPOSITORY https://github.com/arguelles/nuSQuIDS.git
    GIT_TAG master
    BINARY_DIR = <SOURCE_DIR>
    CONFIGURE_HANDLED_BY_BUILD true
    CONFIGURE_COMMAND <SOURCE_DIR>/configure --with-gsl=${GSL_ROOT_DIR} --with-hdf5=${HDF5_INCLUDE_DIRS}/../ --with-squids=${CMAKE_CURRENT_BINARY_DIR}/SQUIDS-prefix --prefix=<INSTALL_DIR>
    BUILD_COMMAND make
    INSTALL_COMMAND make install
  )
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

install(DIRECTORY Configs DESTINATION ${CMAKE_INSTALL_PREFIX}/)
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
    Configs
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
