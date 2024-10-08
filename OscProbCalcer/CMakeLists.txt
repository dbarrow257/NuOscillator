set(HEADERS OscProbCalcerBase.h OscProbCalcerFactory.h)

add_library(OscProbCalcer SHARED OscProbCalcerBase.cpp OscProbCalcerFactory.cpp)

target_link_libraries(OscProbCalcer yaml-cpp NuOscillatorCompilerOptions)

target_include_directories(OscProbCalcer PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/../>
    $<INSTALL_INTERFACE:include>
)

install(TARGETS yaml-cpp
        EXPORT NuOscillator-targets
        LIBRARY DESTINATION lib/
        PUBLIC_HEADER DESTINATION include/yaml-cpp
        )

install(DIRECTORY ${yaml-cpp_SOURCE_DIR}/include/
    DESTINATION include/
)

if(${UseCUDAProb3} EQUAL 1)
	   target_sources(OscProbCalcer PRIVATE OscProbCalcer_CUDAProb3.cpp)
	   if(${UseGPU} EQUAL 1)
	     set_source_files_properties(OscProbCalcer_CUDAProb3.cpp PROPERTIES LANGUAGE CUDA)
	   endif()
	   target_include_directories(OscProbCalcer PRIVATE ${CUDAProb3_SOURCE_DIR})
	   list(APPEND HEADERS OscProbCalcer_CUDAProb3.h)
endif()

if(${UseCUDAProb3Linear} EQUAL 1)
	   target_link_libraries(OscProbCalcer CUDAProb3Beam)
	   target_sources(OscProbCalcer PRIVATE OscProbCalcer_CUDAProb3Linear.cpp)
	   if(${UseGPU} EQUAL 1)
            set_source_files_properties(OscProbCalcer_CUDAProb3Linear.cpp PROPERTIES LANGUAGE CUDA)
            target_include_directories(OscProbCalcer PRIVATE ${CUDAProb3_SOURCE_DIR})
	   endif()
	   list(APPEND HEADERS OscProbCalcer_CUDAProb3Linear.h)
endif()

if(${UseProbGPULinear} EQUAL 1)
	   target_link_libraries(OscProbCalcer ProbGPU)
	   set_property(TARGET ProbGPU PROPERTY CUDA_ARCHITECTURES ${CMAKE_CUDA_ARCHITECTURES})
	   target_sources(OscProbCalcer PRIVATE OscProbCalcer_ProbGPULinear.cpp)
	   list(APPEND HEADERS OscProbCalcer_ProbGPULinear.h)
endif()

if(${UseProb3ppLinear} EQUAL 1)
	   target_link_libraries(OscProbCalcer Prob3plusplus)
	   target_sources(OscProbCalcer PRIVATE OscProbCalcer_Prob3ppLinear.cpp)
	   list(APPEND HEADERS OscProbCalcer_Prob3ppLinear.h)
endif()

if(${UseNuFASTLinear} EQUAL 1)
	   include_directories(${NuFAST_SOURCE_DIR})
	   target_sources(OscProbCalcer PRIVATE OscProbCalcer_NuFASTLinear.cpp)
	   list(APPEND HEADERS OscProbCalcer_NuFASTLinear.h)
endif()

if (${UseGPU} EQUAL 1)
   set_target_properties(OscProbCalcer PROPERTIES
				    CUDA_SEPARABLE_COMPILATION ON
				    LINKER_LANGUAGE CUDA
				    EXPORT_NAME OscProbCalcer)
endif()

set_target_properties(OscProbCalcer PROPERTIES
    PUBLIC_HEADER "${HEADERS}"
    EXPORT_NAME OscProbCalcer)

install(TARGETS OscProbCalcer
        EXPORT NuOscillator-targets
        LIBRARY DESTINATION lib/
        PUBLIC_HEADER DESTINATION include/OscProbCalcer)

add_library(OscProbCalcer::OscProbCalcer ALIAS OscProbCalcer)
