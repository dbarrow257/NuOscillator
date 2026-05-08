if(NOT DEFINED GLoBES_BRANCH)
  set(GLoBES_BRANCH "3.2.18")
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
