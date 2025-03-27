set(HEADERS OscProbCalcerBase.h OscProbCalcerFactory.h)

add_library(OscProbCalcer SHARED OscProbCalcerBase.cpp OscProbCalcerFactory.cpp)

target_link_libraries(OscProbCalcer yaml-cpp NuOscillatorCompilerOptions ROOT::Core)

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
  target_link_libraries(OscProbCalcer ROOT::Core ROOT::RIO)
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

if(${UseNuSQUIDSLinear} EQUAL 1)
  include_directories(${SQUIDS_INC_DIR})
  include_directories(${NUSQUIDS_INC_DIR})
  include_directories(${NUSQUIDS_SRC_DIR})

  set(NUSQUIDSEXAMPLES nuSQUIDSDecoh.o)

  target_link_libraries(OscProbCalcer GSL::gsl GSL::gslcblas HDF5::HDF5 SQUIDS::SQUIDS NUSQUIDS::NUSQUIDS)
  foreach(VAR ${NUSQUIDSEXAMPLES})
    target_link_libraries(OscProbCalcer ${NUSQUIDS_SRC_DIR}/build/${VAR})
  endforeach(VAR)
  target_sources(OscProbCalcer PRIVATE OscProbCalcer_NuSQUIDSLinear.cpp)
  list(APPEND HEADERS OscProbCalcer_NuSQUIDSLinear.h)

  add_dependencies(OscProbCalcer SQUIDS)
  add_dependencies(OscProbCalcer NUSQUIDS)
endif()

if(${UseOscProb} EQUAL 1)
  target_sources(OscProbCalcer PRIVATE OscProbCalcer_OscProb.cpp)
  target_link_libraries(OscProbCalcer ${OscProb_SOURCE_DIR}/lib/libOscProb.so)
  target_link_libraries(OscProbCalcer Eigen OSCPROB ROOT::Core)
  list(APPEND HEADERS OscProbCalcer_OscProb.h)
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
