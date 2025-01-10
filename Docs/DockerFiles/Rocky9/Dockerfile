#To run use:  docker build -t nuoscillator .
# KS: Get glorious container from Kamil which will work as a base
FROM kamilskwarczynski/nukamil:latest AS nuoscillator_build

# Add a label for the author
LABEL website="https://dbarrow257.github.io/NuOscillator/"
LABEL org.opencontainers.image.description="Official NuOscillator container"

ARG CMAKE_OPTIONS
ENV CMAKE_OPTIONS=${CMAKE_OPTIONS:-DUseGPU=1}

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
RUN cmake -DUseGPU=1 -DUseMultithreading=1 -DUseDoubles=1 -DUseCUDAProb3=1 -DUseCUDAProb3Linear=0 -DUseProb3ppLinear=1 -DUseNuFASTLinear=1 -DUseProbGPULinear=1 -DUseOscProb=1 ${NuOscillator_WORK_DIR}
RUN make VERBOSE=1 && make install

# KS: Need to set them here, otherwise container using this container will not be able to find NuOscillator
ENV NuOscillator_ROOT=${NuOscillator_INSTALL_DIR}
ENV PATH=${NuOscillator_ROOT}/bin:${PATH} \
    LD_LIBRARY_PATH=${NuOscillator_ROOT}/lib:${LD_LIBRARY_PATH}

# Start from NuOscillator install dir
WORKDIR ${NuOscillator_INSTALL_DIR}
