SET(OSCPROBCALCER_VERSION 0.0)

if (DEFINED ENV{UseGPU} AND $ENV{UseGPU} EQUAL 1)
   project(OSCPROBCALCER VERSION ${OSCPROBCALCER_VERSION} LANGUAGES CXX C CUDA)
else()
   project(OSCPROBCALCER VERSION ${OSCPROBCALCER_VERSION} LANGUAGES CXX)
endif()

#============================================================================================================
#Interface library used for setting compiler options

add_library(NuOscillatorInterface INTERFACE)

foreach(VAR UseGPU UseMultithreading UseDoubles UseCUDAProb3 UseCUDAProb3Linear UseProb3ppLinear UseProbGPULinear)
	    target_compile_options(NuOscillatorInterface INTERFACE -D${VAR}=$ENV{${VAR}})
endforeach(VAR)

#============================================================================================================

add_library(OscProbCalcer SHARED OscProbCalcerBase.cpp)
target_link_libraries(OscProbCalcer NuOscillatorInterface)
target_include_directories(OscProbCalcer PUBLIC ../Constants)

if(DEFINED ENV{UseCUDAProb3} AND $ENV{UseCUDAProb3} EQUAL 1)
	   target_sources(OscProbCalcer PRIVATE OscProbCalcer_CUDAProb3.cpp)
	   if(DEFINED ENV{UseGPU} AND $ENV{UseGPU} EQUAL 1)
	   	      set_source_files_properties(OscProbCalcer_CUDAProb3.cpp PROPERTIES LANGUAGE CUDA)
	   endif()
	   target_include_directories(OscProbCalcer PRIVATE ${CUDAProb3_SOURCE_DIR})
endif()

if(DEFINED ENV{UseCUDAProb3Linear} AND $ENV{UseCUDAProb3Linear} EQUAL 1)
	   target_link_libraries(OscProbCalcer CUDAProb3Beam CUDAProb3Atmos)
	   target_sources(OscProbCalcer PRIVATE OscProbCalcer_CUDAProb3Linear.cpp)
endif()

if(DEFINED ENV{UseProbGPULinear} AND $ENV{UseProbGPULinear} EQUAL 1)
	   target_link_libraries(OscProbCalcer ProbGPU)
	   target_sources(OscProbCalcer PRIVATE OscProbCalcer_ProbGPULinear.cpp)
endif()

if(DEFINED ENV{UseProb3ppLinear} AND $ENV{UseProb3ppLinear} EQUAL 1)
	   target_link_libraries(OscProbCalcer Prob3plusplus)
	   target_sources(OscProbCalcer PRIVATE OscProbCalcer_Prob3ppLinear.cpp)
endif()

if (DEFINED ENV{UseGPU} AND $ENV{UseGPU} EQUAL 1)
   set_property(TARGET ProbGPU PROPERTY CUDA_ARCHITECTURES 35 52 60 61 70 75)

   set_target_properties(OscProbCalcer PROPERTIES 
				    CUDA_SEPARABLE_COMPILATION ON 
				    LINKER_LANGUAGE CUDA
				    EXPORT_NAME OscProbCalcer)
else()
   set_target_properties(OscProbCalcer PROPERTIES
                                    EXPORT_NAME OscProbCalcer)
endif()

add_library(OscProbCalcer::OscProbCalcer ALIAS OscProbCalcer)
