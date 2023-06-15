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

#================================================================================================
#Setup Deps

include(CPM)
CPMFindPackage(
	NAME ProbGPU
	GIT_TAG main
	GITHUB_REPOSITORY dbarrow257/ProbGPU
)

if(NOT TARGET ProbGPU::All)
       message(FATAL_ERROR "Could not find Prob3plusplus")
endif()

#================================================================================================
#Subdirectories to add

add_subdirectory(OscProbCalcer)
add_subdirectory(Oscillator)
add_subdirectory(Apps)
