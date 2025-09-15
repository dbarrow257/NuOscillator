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

set(NUOSCILLATOR_VERSION 1.4.3)

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

# download CPM.cmake
file(
  DOWNLOAD
  https://github.com/cpm-cmake/CPM.cmake/releases/download/v0.42.0/CPM.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/cmake/CPM.cmake
)
include(${CMAKE_CURRENT_BINARY_DIR}/cmake/CPM.cmake)

DefineEnabledRequiredSwitch(UseMultithreading 1)
DefineEnabledRequiredSwitch(UseDoubles 1)

DefineEnabledRequiredSwitch(UseCUDAProb3 0)
DefineEnabledRequiredSwitch(UseCUDAProb3Linear 0)
DefineEnabledRequiredSwitch(UseProb3ppLinear 0)
DefineEnabledRequiredSwitch(UseProbGPULinear 0)
DefineEnabledRequiredSwitch(UseNuFASTLinear 0)
DefineEnabledRequiredSwitch(UseNuSQUIDSLinear 0)
DefineEnabledRequiredSwitch(UseOscProb 0)
DefineEnabledRequiredSwitch(UseOscLib 0)
DefineEnabledRequiredSwitch(UseGLoBESLinear  0)

# KS: NuFastLinear by default
if(UseCUDAProb3        EQUAL 0 AND
    UseCUDAProb3Linear EQUAL 0 AND
    UseProb3ppLinear   EQUAL 0 AND
    UseProbGPULinear   EQUAL 0 AND
    UseNuFASTLinear    EQUAL 0 AND
    UseNuSQUIDSLinear  EQUAL 0 AND
    UseOscProb         EQUAL 0 AND
    UseOscLib          EQUAL 0 AND
    UseGLoBESLinear    EQUAL 0 )
  set(UseNuFASTLinear 1)
endif()

LIST(APPEND ALL_Engines
  CUDAProb3
  CUDAProb3Linear
  Prob3ppLinear
  ProbGPULinear
  NuFASTLinear
  NuSQUIDSLinear
  OscProb
  OscLib
  GLoBESLinear
)

########################################################
#Default Branches

if(NOT DEFINED CUDAProb3_BRANCH)
  set(CUDAProb3_BRANCH "v1.2.0")
endif()

if(NOT DEFINED CUDAProb3Linear_BRANCH)
  set(CUDAProb3Linear_BRANCH "v1.1.0")
endif()

if(NOT DEFINED Prob3ppLinear_BRANCH)
  set(Prob3ppLinear_BRANCH "v3r20")
endif()

if(NOT DEFINED ProbGPULinear_BRANCH)
  set(ProbGPULinear_BRANCH "v1.0.0")
endif()

if(NOT DEFINED NuFASTLinear_BRANCH)
  set(NuFASTLinear_BRANCH "v1.0")
endif()

if(NOT DEFINED NuSQUIDSLinear_BRANCH)
  set(NuSQUIDSLinear_BRANCH "v1.12.0")
endif()

if(NOT DEFINED OscProb_BRANCH)
  set(OscProb_BRANCH "v2.1.2")
endif()

if(NOT DEFINED OscLib_BRANCH)
  set(OscLib_BRANCH "v1.0.1")
endif()

if(NOT DEFINED GLoBES_BRANCH)
  set(GLoBES_BRANCH "3.2.18")
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

#================================================================================================
#Setup Deps
set(YAML_CPP_VERSION 0.8.0)
CPMAddPackage(
    NAME yaml-cpp
    VERSION ${YAML_CPP_VERSION}
    GITHUB_REPOSITORY "jbeder/yaml-cpp"
    GIT_TAG "${YAML_CPP_VERSION}"
    OPTIONS
      "YAML_CPP_INSTALL ON"
      "YAML_CPP_BUILD_TESTS OFF"
      "YAML_CPP_BUILD_CONTRIB OFF"
      "YAML_BUILD_SHARED_LIBS ON"
)

if(NOT TARGET yaml-cpp::yaml-cpp)
  cmessage(FATAL_ERROR "NuOscillator Expected dependency target: yaml-cpp::yaml-cpp")
endif()

find_package(ROOT 6.18 REQUIRED)
STRING(STRIP "${ROOT_CXX_FLAGS}" ROOT_CXX_FLAGS_LIST)
STRING(REPLACE " " ";" ROOT_CXX_FLAGS_LIST ${ROOT_CXX_FLAGS_LIST})

list (FIND ROOT_CXX_FLAGS_LIST "-std=c++14" CPP14_INDEX)
list (FIND ROOT_CXX_FLAGS_LIST "-std=c++1y" CPP1Y_INDEX)
list (FIND ROOT_CXX_FLAGS_LIST "-std=c++17" CPP17_INDEX)
list (FIND ROOT_CXX_FLAGS_LIST "-std=c++1z" CPP1Z_INDEX)
list (FIND ROOT_CXX_FLAGS_LIST "-std=c++20" CPP20_INDEX)

if (CPP14_INDEX GREATER -1)
  SET(ROOT_CXX_STANDARD 14)
elseif (CPP1Y_INDEX GREATER -1)
  SET(ROOT_CXX_STANDARD 14)
elseif (CPP17_INDEX GREATER -1)
  SET(ROOT_CXX_STANDARD 17)
elseif (CPP1Z_INDEX GREATER -1)
  SET(ROOT_CXX_STANDARD 17)
elseif (CPP20_INDEX GREATER -1)
  SET(ROOT_CXX_STANDARD 20)
endif()

if(DEFINED ROOT_CXX_STANDARD AND NOT ROOT_CXX_STANDARD EQUAL CMAKE_CXX_STANDARD)
  set(CMAKE_CXX_STANDARD ${ROOT_CXX_STANDARD})
  cmessage(STATUS "Set CXX standard due to ROOT: \"${ROOT_CXX_STANDARD}\"")
endif()

if(${UseCUDAProb3} EQUAL 1)
  CPMFindPackage(
    NAME CUDAProb3
    GITHUB_REPOSITORY "dbarrow257/CUDAProb3"
    GIT_TAG ${CUDAProb3_BRANCH}
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
      NAME CUDAProb3Linear
      GITHUB_REPOSITORY "mach3-software/CUDAProb3"
      GIT_TAG ${CUDAProb3Linear_BRANCH}
      OPTIONS
      "GPU_ON ON"
      "CMAKE_CUDA_ARCHITECTURES ${CMAKE_CUDA_ARCHITECTURES_STRING}"
    )
   # KS: This is bit hacky but CUDAProb3Linear is a mess. Whole functions are written in header file
   # and someone wrote them to be hidden under GPU_ON and wrote shitty cmake which do not propagate this definition.
   # Therefore NuOscillator is fixing other people mess
   target_compile_definitions(NuOscillatorCompilerOptions INTERFACE GPU_ON)
  else()
    CPMFindPackage(
      NAME CUDAProb3Linear
      GITHUB_REPOSITORY "mach3-software/CUDAProb3"
      GIT_TAG ${CUDAProb3Linear_BRANCH}
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
      GIT_TAG ${ProbGPULinear_BRANCH}
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
    GITHUB_REPOSITORY rogerwendell/Prob3plusplus
    GIT_TAG ${Prob3ppLinear_BRANCH}
  )

  if(NOT TARGET Prob3plusplus)
    cmessage(FATAL_ERROR "Could not find Prob3plusplus")
  endif()

  # KS: Add additional compilation flags to be used for Prob3++
  target_compile_options(Prob3plusplus PRIVATE ${NuOscillator_Compiler_Flags})
  target_compile_definitions(NuOscillatorCompilerOptions INTERFACE UseProb3ppLinear=1)

  install(TARGETS Prob3plusplus
    EXPORT NuOscillator-targets
    LIBRARY DESTINATION lib/)
endif()

if(${UseNuFASTLinear} EQUAL 1)
  CPMAddPackage(
    NAME NuFAST
    GIT_TAG ${NuFASTLinear_BRANCH}
    GITHUB_REPOSITORY PeterDenton/NuFast-LBL
    DOWNLOAD_ONLY
  )
  target_compile_definitions(NuOscillatorCompilerOptions INTERFACE UseNuFASTLinear=1)
endif()

if(${UseOscProb} EQUAL 1)
  CPMAddPackage(
    NAME OscProb
    GITHUB_REPOSITORY joaoabcoelho/OscProb
    GIT_TAG ${OscProb_BRANCH}
  )

  CPMAddPackage(
    NAME Eigen
    VERSION 3.4.0
    URL https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz
    DOWNLOAD_ONLY YES
  )
  if(Eigen_ADDED)
    add_library(Eigen INTERFACE IMPORTED)
    target_include_directories(Eigen INTERFACE ${Eigen_SOURCE_DIR})
  endif()

  execute_process(COMMAND make WORKING_DIRECTORY ${OscProb_SOURCE_DIR})
  set(OSCPROB_LIB_DIR ${OscProb_SOURCE_DIR}/lib)

  add_library(OSCPROB INTERFACE IMPORTED)
  target_include_directories(OSCPROB INTERFACE ${OscProb_SOURCE_DIR})
  target_compile_definitions(NuOscillatorCompilerOptions INTERFACE UseOscProb=1)
endif()

if(${UseOscLib} EQUAL 1)
  # Check if C++ standard is less than 17
  if(CMAKE_CXX_STANDARD LESS 17)
    cmessage(FATAL_ERROR "OscLib requires at least C++17, but you are using C++${CMAKE_CXX_STANDARD}.")
  endif()

  find_package(Boost REQUIRED)
  find_package(GSL REQUIRED)

  CPMAddPackage(
    NAME Eigen
    VERSION 3.4.0
    URL https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz
    DOWNLOAD_ONLY YES
  )
  if(Eigen_ADDED)
    add_library(Eigen INTERFACE IMPORTED)
    target_include_directories(Eigen INTERFACE ${Eigen_SOURCE_DIR})
  endif()

  CPMAddPackage(
    NAME OscLib
    GITHUB_REPOSITORY cafana/OscLib
    GIT_TAG ${OscLib_BRANCH}
    DOWNLOAD_ONLY YES
  )
  # KS: Collect OscLib sources and header files
  file(GLOB OSCLIB_SOURCES
    ${OscLib_SOURCE_DIR}/OscLib/*.cxx
  )
  file(GLOB OSCLIB_HEADERS
    ${OscLib_SOURCE_DIR}/OscLib/*.h
  )
  add_library(OscLib SHARED
    ${OSCLIB_SOURCES}
  )

  set_target_properties(OscLib PROPERTIES
    PUBLIC_HEADER "${OSCLIB_HEADERS}"
    EXPORT_NAME OscLib
  )

  target_link_libraries(OscLib ROOT::Core Eigen GSL::gsl)
  # Properly scoped include dirs
  target_include_directories(OscLib PUBLIC
      $<BUILD_INTERFACE:${OscLib_SOURCE_DIR}>
      $<INSTALL_INTERFACE:include>
  )
  install(TARGETS OscLib
    EXPORT NuOscillator-targets
    LIBRARY DESTINATION lib/
    PUBLIC_HEADER DESTINATION include/
  )
  target_compile_definitions(NuOscillatorCompilerOptions INTERFACE UseOscLib=1)
endif()

if(${UseNuSQUIDSLinear} EQUAL 1)
  find_package(GSL REQUIRED)

  set(BUILTIN_HDF5 FALSE)
  find_package(HDF5 QUIET COMPONENTS CXX)
  include(ExternalProject)

  if(NOT TARGET HDF5::HDF5)
    SET(HDF5_DIR ${CMAKE_CURRENT_BINARY_DIR}/HDF5_project-prefix/)
    SET(HDF5_LIB_DIR ${CMAKE_CURRENT_BINARY_DIR}/HDF5_project-prefix/lib)
    ExternalProject_Add(HDF5_project
        GIT_REPOSITORY https://github.com/HDFGroup/hdf5.git
        GIT_TAG hdf5-1.14.6
        UPDATE_DISCONNECTED TRUE
        DOWNLOAD_NO_PROGRESS TRUE
        CMAKE_ARGS
            -DCMAKE_INSTALL_PREFIX=${HDF5_DIR}
            -DBUILD_SHARED_LIBS=ON
            -DHDF5_BUILD_EXAMPLES=OFF
            -DHDF5_ENABLE_PARALLEL=OFF
        BUILD_COMMAND make -j${NPROC}
        INSTALL_COMMAND make install
    )

    file(MAKE_DIRECTORY "${HDF5_DIR}/include")

    # Imported target pointing to the installed library
    add_library(HDF5::HDF5 SHARED IMPORTED)
    set_target_properties(HDF5::HDF5 PROPERTIES
        IMPORTED_LOCATION "${HDF5_LIB_DIR}/libhdf5.so"
        INTERFACE_INCLUDE_DIRECTORIES "${HDF5_DIR}/include"
    )
    set(BUILTIN_HDF5 TRUE)
  endif()

  ExternalProject_Add(SQUIDS
    GIT_REPOSITORY https://github.com/jsalvado/SQuIDS.git
    GIT_TAG v1.3.1
    BUILD_IN_SOURCE = 1
    CONFIGURE_HANDLED_BY_BUILD false
    CONFIGURE_COMMAND <SOURCE_DIR>/configure --with-gsl=${GSL_ROOT_DIR} --prefix=<INSTALL_DIR>
    BUILD_COMMAND make
    INSTALL_COMMAND make install
  )

  SET(SQUIDS_INC_DIR ${CMAKE_CURRENT_BINARY_DIR}/SQUIDS-prefix/include)
  SET(SQUIDS_LIB_DIR ${CMAKE_CURRENT_BINARY_DIR}/SQUIDS-prefix/lib)

  add_library(SQUIDS::SQUIDS SHARED IMPORTED)
  set_target_properties(SQUIDS::SQUIDS PROPERTIES IMPORTED_LOCATION ${SQUIDS_LIB_DIR}/libSQuIDS.so)

  if(BUILTIN_HDF5)
      set(NUSQUIDS_HDF5_DIR "${HDF5_DIR}/")
      set(NUSQUIDS_DEPENDS HDF5_project SQUIDS)
  else()
      set(NUSQUIDS_HDF5_DIR "${HDF5_INCLUDE_DIRS}/..")
      set(NUSQUIDS_DEPENDS SQUIDS)
  endif()

  ExternalProject_Add(NUSQUIDS
    GIT_REPOSITORY https://github.com/arguelles/nuSQuIDS.git
    GIT_TAG ${NuSQUIDSLinear_BRANCH}
    BUILD_IN_SOURCE = 1
    CONFIGURE_HANDLED_BY_BUILD false
    CONFIGURE_COMMAND <SOURCE_DIR>/configure
        --with-gsl=${GSL_ROOT_DIR}
        --with-hdf5=${NUSQUIDS_HDF5_DIR}
        --with-squids=${CMAKE_CURRENT_BINARY_DIR}/SQUIDS-prefix
        --prefix=<INSTALL_DIR>
    BUILD_COMMAND make examples
    INSTALL_COMMAND make install
    DEPENDS ${NUSQUIDS_DEPENDS}    #ensure both are finished first
   )

  SET(NUSQUIDS_SRC_DIR ${CMAKE_CURRENT_BINARY_DIR}/NUSQUIDS-prefix/src/NUSQUIDS/)
  SET(NUSQUIDS_INC_DIR ${CMAKE_CURRENT_BINARY_DIR}/NUSQUIDS-prefix/include/)
  SET(NUSQUIDS_LIB_DIR ${CMAKE_CURRENT_BINARY_DIR}/NUSQUIDS-prefix/lib/)

  add_dependencies(NUSQUIDS SQUIDS)
  add_library(NUSQUIDS::NUSQUIDS SHARED IMPORTED)
  set_target_properties(NUSQUIDS::NUSQUIDS PROPERTIES IMPORTED_LOCATION ${NUSQUIDS_LIB_DIR}/libnuSQuIDS.so)
  target_compile_definitions(NuOscillatorCompilerOptions INTERFACE UseNuSQUIDSLinear=1)
endif()

if(${UseGLoBESLinear} EQUAL 1)
  include(ExternalProject)
  set(GLoBES_URL "https://www.mpi-hd.mpg.de/personalhomes/globes/download/globes-${GLoBES_BRANCH}.tar.gz")
  set(GLOBES_PREFIX "${CMAKE_BINARY_DIR}/globes-prefix")
  set(GLOBES_LIB_DIR "${GLOBES_PREFIX}/lib/")
  set(GLOBES_INC_DIR "${GLOBES_PREFIX}/include/")

  ExternalProject_Add(
    globes_project
    URL ${GLoBES_URL}
    URL_HASH SHA256=7edd0fea67800fc1f2d5309a47aba93d4afb5710a10cb40ceddec07079041db3

    PREFIX ${GLOBES_PREFIX}
    CONFIGURE_COMMAND <SOURCE_DIR>/configure --prefix=${GLOBES_PREFIX} CXXFLAGS="-O3 -g"
    BUILD_COMMAND make -C <BINARY_DIR>
    INSTALL_COMMAND make -C <BINARY_DIR> install
  )
  file(MAKE_DIRECTORY "${GLOBES_INC_DIR}/globes")

  add_library(GLOBES::GLOBES SHARED IMPORTED)
  set_target_properties(GLOBES::GLOBES PROPERTIES
      IMPORTED_LOCATION "${GLOBES_LIB_DIR}/libglobes.so.8.3.2"
      INTERFACE_INCLUDE_DIRECTORIES "${GLOBES_INC_DIR}"
  )

  # --- SNU Integration ---
  set(SNU_URL "https://www.mpi-hd.mpg.de/personalhomes/globes/tools/snu-1.1.tar.gz")
  set(SNU_PREFIX "${CMAKE_BINARY_DIR}/snu-prefix")
  set(SNU_INC_DIR "${SNU_PREFIX}/include/")
  set(SNU_LIB_DIR "${SNU_PREFIX}/lib/")

  ExternalProject_Add(
    snu_project
    URL ${SNU_URL}
    PREFIX ${SNU_PREFIX}
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ${CMAKE_COMMAND} -E echo "Building SNU..."
      COMMAND ${CMAKE_C_COMPILER} -O3 -g -fPIC -shared -I${GLOBES_INC_DIR} -o <BINARY_DIR>/libsnu.so <SOURCE_DIR>/snu.c
    INSTALL_COMMAND
      ${CMAKE_COMMAND} -E make_directory ${SNU_LIB_DIR}
      COMMAND ${CMAKE_COMMAND} -E copy <BINARY_DIR>/libsnu.so ${SNU_LIB_DIR}/
      COMMAND ${CMAKE_COMMAND} -E make_directory ${SNU_INC_DIR}
      COMMAND ${CMAKE_COMMAND} -E copy <SOURCE_DIR>/snu.h ${SNU_INC_DIR}/
    BUILD_IN_SOURCE 1
    DEPENDS globes_project
  )
  file(MAKE_DIRECTORY "${SNU_INC_DIR}/")
  add_library(SNU::SNU SHARED IMPORTED)
  set_target_properties(SNU::SNU PROPERTIES
    IMPORTED_LOCATION "${SNU_LIB_DIR}/libsnu.so"
    INTERFACE_INCLUDE_DIRECTORIES "${SNU_INC_DIR}"
  )
  add_dependencies(SNU::SNU GLOBES::GLOBES)
  target_compile_definitions(NuOscillatorCompilerOptions INTERFACE UseGLoBESLinear=1)
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

# KS: This extract relative path to CudaProb and others. This allow in setup to use ${NUOSCILLATOR_ROOT}/RELATIVE_PART
# this way even if we move directory we can still use NuOscillator
if(CUDAProb3_SOURCE_DIR)
  file(RELATIVE_PATH CUDAProb3_RELATIVE_PATH
       "${CMAKE_INSTALL_PREFIX}"
       "${CUDAProb3_SOURCE_DIR}")
endif()

if(HDF5_LIB_DIR)
  file(RELATIVE_PATH HDF5_REL_LIB_PATH
       "${CMAKE_INSTALL_PREFIX}"
       "${HDF5_LIB_DIR}")
endif()

if(SQUIDS_LIB_DIR)
  file(RELATIVE_PATH SQUIDS_REL_LIB_PATH
       "${CMAKE_INSTALL_PREFIX}"
       "${SQUIDS_LIB_DIR}")
endif()

if(NUSQUIDS_LIB_DIR)
  file(RELATIVE_PATH NUSQUIDS_REL_LIB_PATH
       "${CMAKE_INSTALL_PREFIX}"
       "${NUSQUIDS_LIB_DIR}")
endif()

if(OSCPROB_LIB_DIR)
  file(RELATIVE_PATH OSCPROB_REL_LIB_PATH
       "${CMAKE_INSTALL_PREFIX}"
       "${OSCPROB_LIB_DIR}")
endif()

if(OSCLIB_LIB_DIR)
    file(RELATIVE_PATH OSCLIB_REL_LIB_PATH
       "${CMAKE_INSTALL_PREFIX}"
       "${OSCLIB_LIB_DIR}")
endif()

if(GLOBES_LIB_DIR)
  file(RELATIVE_PATH GLOBES_REL_LIB_PATH
       "${CMAKE_INSTALL_PREFIX}"
       "${GLOBES_LIB_DIR}")
endif()

if(SNU_LIB_DIR)
  file(RELATIVE_PATH SNU_REL_LIB_PATH
       "${CMAKE_INSTALL_PREFIX}"
       "${SNU_LIB_DIR}")
endif()

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

include(nuoscillator-config)
