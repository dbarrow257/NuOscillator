BASEDIR=${PWD}

#ProbGPU
if [ ${UseGPU} == 1 ]; then
    git clone git@github.com:dbarrow257/ProbGPU.git
    cd ./ProbGPU
    make
    cd ${BASEDIR}
fi

#Prob3pp
git clone https://github.com/rogerwendell/Prob3plusplus.git
cd ./Prob3plusplus
make
make shared
cd ${BASEDIR}

#CUDAProb3
git clone git@github.com:/dbarrow257/CUDAProb3.git
cd CUDAProb3
git checkout develop
cd ${BASEDIR}
