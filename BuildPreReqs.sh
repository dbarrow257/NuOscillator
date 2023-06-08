BASEDIR=${PWD}

#ProbGPU
if [ ${UseProbGPULinear} == 1 ]; then
    git clone git@github.com:dbarrow257/ProbGPU.git
    cd ./ProbGPU
    make
    cd ${BASEDIR}
fi

#Prob3pp
if [ ${UseProb3ppLinear} == 1 ]; then
    git clone https://github.com/rogerwendell/Prob3plusplus.git
    cd ./Prob3plusplus
    make
    make shared
    cd ${BASEDIR}
fi

#CUDAProb3
if [ ${UseCUDAProb3} == 1 ]; then
    git clone git@github.com:/dbarrow257/CUDAProb3.git
    cd CUDAProb3
    git checkout develop
    cd ${BASEDIR}
fi

#CUDAProb3
if [ ${UseCUDAProb3Linear} == 1 ]; then
    git clone git@github.com:/mach3-software/CUDAProb3.git CUDAProb3Linear
    cd CUDAProb3Linear
    git checkout feature_cleanup
    mkdir build
    cd build
    #The CMAKE build system in this directory is pretty screwy and the library is only built if GPU options are turned on
    #May also need a 'set(CMAKE_CXX_STANDARD 11)' dropped into the CMAKELists.txt file if using older compilers
    cmake3 .. -DGPU_ON=1
    make
    cd ${BASEDIR}
fi
