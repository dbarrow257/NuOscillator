#================================================================================================
#Project details

cmake_minimum_required(VERSION 3.14 FATAL_ERROR)

find_program(CMAKE_C_COMPILER NAMES $ENV{CC} gcc PATHS ENV PATH NO_DEFAULT_PATH)
find_program(CMAKE_CXX_COMPILER NAMES $ENV{CXX} g++ PATHS ENV PATH NO_DEFAULT_PATH)

if(NOT DEFINED CMAKE_CXX_STANDARD)
  SET(CMAKE_CXX_STANDARD 11)
endif()

SET(NUOSCILLATOR_VERSION 0.0)
project(NUOSCILLATOR VERSION ${NUOSCILLATOR_VERSION} LANGUAGES CXX)

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

foreach(VAR UseGPU UseMultithreading UseDoubles UseCUDAProb3 UseCUDAProb3Linear UseProb3ppLinear UseProbGPULinear)
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
	find_package(OpenMP)
	if (OPENMP_FOUND)
	    SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
	    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    	    SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
	endif()
endif()

#================================================================================================
#Setup Deps

include(CPM)

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
	   	  message(FATAL_ERROR "Could not download CUDAProb3")
	   endif()

	   add_compile_definitions(UseCUDAProb3=1)
endif()

if(${UseCUDAProb3Linear} EQUAL 1)
	if(${UseGPU} EQUAL 1)
    		   CPMFindPackage(
			NAME CUDAProb3
        		GITHUB_REPOSITORY "mach3-software/CUDAProb3"
        		GIT_TAG "feature_cleanup"
			OPTIONS "GPU_ON ON"
    			)
	else()
	           CPMFindPackage(
			NAME CUDAProb3
                        GITHUB_REPOSITORY "mach3-software/CUDAProb3"
                        GIT_TAG "feature_cleanup"
			OPTIONS "CPU_ONLY ON"
			)
	endif()
	if(NOT TARGET CUDAProb3Beam OR NOT TARGET CUDAProb3Atmos)
	       message(FATAL_ERROR "Could not find CUDAProb3Linear")
	endif()

	add_compile_definitions(UseCUDAProb3Linear=1)
endif()

if(${UseGPU} EQUAL 1)
	if(${UseProbGPULinear} EQUAL 1)
		   CPMFindPackage(
		   	NAME ProbGPU
			GIT_TAG main   
			GITHUB_REPOSITORY dbarrow257/ProbGPU
			)
		   if(NOT TARGET ProbGPU)
       	       	   	message(FATAL_ERROR "Could not find ProbGPU")
		   endif()

		   add_compile_definitions(UseProbGPULinear=1)
	endif()
endif()

if(${UseProb3ppLinear} EQUAL 1)
	   CPMFindPackage(
		NAME Prob3plusplus
		VERSION 3.10.3
		GIT_TAG v3.10.3
        	GITHUB_REPOSITORY mach3-software/Prob3plusplus
		)
	   if(NOT TARGET Prob3plusplus)
       	   	message(FATAL_ERROR "Could not find Prob3plusplus")
	   endif()

	   add_compile_definitions(UseProb3ppLinear=1)
endif()

#================================================================================================
#Subdirectories to add

add_subdirectory(OscProbCalcer)
add_subdirectory(Oscillator)
add_subdirectory(Apps)

#================================================================================================
#

set(LIBDEST ${CMAKE_INSTALL_PREFIX}/lib)

configure_file(cmake/Templates/setup.NuOscillator.sh.in
  "${PROJECT_BINARY_DIR}/setup.Oscillator.sh" @ONLY)
install(FILES
  "${PROJECT_BINARY_DIR}/setup.Oscillator.sh" DESTINATION
  ${CMAKE_INSTALL_PREFIX})

install(DIRECTORY Configs DESTINATION ${CMAKE_INSTALL_PREFIX}/)
install(DIRECTORY Inputs DESTINATION ${CMAKE_INSTALL_PREFIX}/)

#================================================================================================

#This is to export the target properties of NuOscillator
#Anything that links to "NuOscillator" will get all of these target properties
add_library(NuOscillator INTERFACE)
target_link_libraries(NuOscillator INTERFACE Oscillator::Oscillator OscProbCalcer::OscProbCalcer)
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
