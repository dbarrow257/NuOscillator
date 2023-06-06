export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${PWD}/Prob3plusplus:${PWD}/ProbGPU

#Environment configs
export UseGPU=1
export UseMultithreading=1

#Which calculators to compile
export UseCUDAProb3=1
export UseProbGPULinear=1
export UseProb3ppLinear=1

#ProbGPU only supported when using GPU
if [ ${UseGPU} == 0 ]; then
    export UseProbGPULinear=0
fi

#Set CUDA environment variables
if [ ${UseGPU} == 0 ]; then
    export CUDAPATH=/usr/local/cuda/
    export PATH=${PATH}:${CUDAPATH}/bin/
fi 
