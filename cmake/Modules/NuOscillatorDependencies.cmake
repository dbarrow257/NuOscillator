#Default Branches for Eigen. Some engine requires it and we want ensure consistent version is used
if(NOT DEFINED Eigen_BRANCH)
  set(Eigen_BRANCH "3.4.0")
endif()

# download CPM.cmake
file(
  DOWNLOAD
  https://github.com/cpm-cmake/CPM.cmake/releases/download/v0.42.3/CPM.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/cmake/CPM.cmake
)
include(${CMAKE_CURRENT_BINARY_DIR}/cmake/CPM.cmake)

#================================================================================================
#Setup Deps
set(YAML_CPP_VERSION 0.9.0)
CPMAddPackage(
    NAME yaml-cpp
    VERSION ${YAML_CPP_VERSION}
    GITHUB_REPOSITORY "jbeder/yaml-cpp"
    GIT_TAG "yaml-cpp-${YAML_CPP_VERSION}"
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
