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

message(STATUS "Environment variables used:")
foreach(VAR UseGPU UseMultithreading UseDoubles UseCUDAProb3 UseCUDAProb3Linear UseProb3ppLinear UseProbGPULinear)
	    if(NOT DEFINED ENV{${VAR}})
	    	   message(FATAL_ERROR "${VAR} not defined")
	    endif()
            SET(${VAR} $ENV{${VAR}})
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
	if(DEFINED ENV{UseProbGPULinear} AND $ENV{UseProbGPULinear} EQUAL 1)
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

#add_subdirectory(Constants)
add_subdirectory(OscProbCalcer)
add_subdirectory(Oscillator)
add_subdirectory(Apps)

#================================================================================================
#

configure_file(cmake/Templates/setup.NuOscillator.sh.in
  "${PROJECT_BINARY_DIR}/setup.Oscillator.sh" @ONLY)
install(FILES
  "${PROJECT_BINARY_DIR}/setup.Oscillator.sh" DESTINATION
  ${CMAKE_INSTALL_PREFIX})
