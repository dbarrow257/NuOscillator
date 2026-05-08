if(NOT DEFINED NuSQUIDSLinear_BRANCH)
  set(NuSQUIDSLinear_BRANCH "v1.13.3")
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
