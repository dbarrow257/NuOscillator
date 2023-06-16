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

#================================================================================================
#Setup Deps

include(CPM)

if(DEFINED ENV{UseCUDAProb3} AND $ENV{UseCUDAProb3} EQUAL 1)
	   CPMFindPackage(
		NAME CUDAProb3
		GITHUB_REPOSITORY "dbarrow257/CUDAProb3"
		GIT_TAG "develop"
		DOWNLOAD_ONLY
	   )
	   if(NOT CUDAProb3_ADDED)
	   	  message(FATAL_ERROR "Could not download CUDAProb3")
	   endif()
endif()

if(DEFINED ENV{UseCUDAProb3Linear} AND $ENV{UseCUDAProb3Linear} EQUAL 1)
	if(DEFINED ENV{UseGPU} AND $ENV{UseGPU} EQUAL 1)
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
endif()

if(DEFINED ENV{UseGPU} AND $ENV{UseGPU} EQUAL 1)
	if(DEFINED ENV{UseProbGPU} AND $ENV{UseProbGPU} EQUAL 1)
		   CPMFindPackage(
		   	NAME ProbGPU
			GIT_TAG main   
			GITHUB_REPOSITORY dbarrow257/ProbGPU
			)
		   if(NOT TARGET ProbGPU)
       	       	   	message(FATAL_ERROR "Could not find ProbGPU")
		   endif()
	endif()
endif()

if(DEFINED ENV{UseProb3ppLinear} AND $ENV{UseProb3ppLinear} EQUAL 1)
	   CPMFindPackage(
		NAME Prob3plusplus
		GIT_TAG main
        	GITHUB_REPOSITORY rogerwendell/Prob3plusplus
		)
	   if(NOT TARGET Prob3plusplus)
       	   	message(FATAL_ERROR "Could not find Prob3plusplus")
	   endif()
endif()

#================================================================================================
#Subdirectories to add

add_subdirectory(OscProbCalcer)
add_subdirectory(Oscillator)
add_subdirectory(Apps)
