#To run use:  docker build -t nuoscillator .
# KS: Get glorious container from Luke which will work as a base
FROM picker24/root_v6_26_10:alma9 AS nuoscillator_build

# Add a label for the author
LABEL website="https://dbarrow257.github.io/NuOscillator/"
LABEL compiler="GNU 11.3.1"
LABEL root_version="v6.26.10"
LABEL org.opencontainers.image.description="Official NuOscillator container"

ARG CMAKE_OPTIONS
ENV CMAKE_OPTIONS=${CMAKE_OPTIONS:-DUseDoubles=1}

# Declare the build argument
ARG NUOSCILLATOR_VERSION
ENV NUOSCILLATOR_VERSION=${NUOSCILLATOR_VERSION:-develop}
ENV NuOscillator_WORK_DIR=/opt/NuOscillator/
ENV NuOscillator_INSTALL_DIR=/opt/NuOscillator/build/

RUN mkdir -p ${NuOscillator_WORK_DIR}

WORKDIR /opt/
# KS: Let's clone NuOscillator
RUN --mount=type=ssh git clone https://github.com/dbarrow257/NuOscillator ${NuOscillator_WORK_DIR}
WORKDIR ${NuOscillator_WORK_DIR}
RUN git checkout ${NUOSCILLATOR_VERSION}

RUN mkdir -p ${NuOscillator_INSTALL_DIR}
WORKDIR ${NuOscillator_INSTALL_DIR}
RUN cmake -DUseGPU=0 -DUseMultithreading=1 ${CMAKE_OPTIONS} -DUseCUDAProb3=0 -DUseCUDAProb3Linear=1 -DUseProb3ppLinear=1 -DUseNuFASTLinear=1 -DUseProbGPULinear=0 -DUseOscProb=1 ${NuOscillator_WORK_DIR}
RUN make VERBOSE=1 && make install

# KS: Need to set them here, otherwise container using this container will not be able to find NuOscillator
ENV NuOscillator_ROOT=${NuOscillator_INSTALL_DIR}
ENV PATH=${NuOscillator_ROOT}/bin:${PATH} \
    LD_LIBRARY_PATH=${NuOscillator_ROOT}/lib:${LD_LIBRARY_PATH}

# Start from NuOscillator install dir
WORKDIR ${NuOscillator_INSTALL_DIR}
