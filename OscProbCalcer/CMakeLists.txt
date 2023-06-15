SET(OSCPROBCALCER_VERSION 0.0)
project(OSCPROBCALCER VERSION ${OSCPROBCALCER_VERSION} LANGUAGES CXX)

add_library(OscProbCalcer SHARED OscProbCalcerBase.cpp)
target_include_directories(OscProbCalcer PUBLIC ../Constants)

find_package(OpenMP)
if (OPENMP_FOUND)
    SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

if(DEFINED ENV{UseCUDAProb3} AND $ENV{UseCUDAProb3} EQUAL 1)
	   message(STATUS "Adding CUDAProb3")
	   target_sources(OscProbCalcer PUBLIC OscProbCalcer_CUDAProb3.cpp)
	   target_include_directories(OscProbCalcer PUBLIC ../CUDAProb3)
endif()

if(DEFINED ENV{UseCUDAProb3Linear} AND $ENV{UseCUDAProb3Linear} EQUAL 1)
	   message(STATUS "Adding CUDAProb3Linear")
	   target_sources(OscProbCalcer PUBLIC OscProbCalcer_CUDAProb3Linear.cpp)
	   target_include_directories(OscProbCalcer PUBLIC ../CUDAProb3Linear)
endif()

if(DEFINED ENV{UseProbGPULinear} AND $ENV{UseProbGPULinear} EQUAL 1)
	   message(STATUS "Adding ProbGPULinear")
	   target_link_libraries(OscProbCalcer -L/opt/ppd/t2k/users/barrowd/DUNE/OscillatorBase/Testbench/ProbGPU -lProbGPU)
	   target_sources(OscProbCalcer PUBLIC OscProbCalcer_ProbGPULinear.cpp)
	   target_include_directories(OscProbCalcer PUBLIC ../ProbGPU)
endif()

if(DEFINED ENV{UseProb3ppLinear} AND $ENV{UseProb3ppLinear} EQUAL 1)
	   message(STATUS "Adding Prob3ppLinear")
	   target_link_libraries(OscProbCalcer -L../Prob3plusplus -lThreeProb_3.20)
	   target_sources(OscProbCalcer PUBLIC OscProbCalcer_Prob3ppLinear.cpp)
	   target_include_directories(OscProbCalcer PUBLIC ../Prob3plusplus)
endif()

set_target_properties(OscProbCalcer PROPERTIES 
				    EXPORT_NAME OscProbCalcer)

target_compile_options(OscProbCalcer PUBLIC -DUseDoubles=1)

add_library(OscProbCalcer::OscProbCalcer ALIAS OscProbCalcer)
