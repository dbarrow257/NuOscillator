export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${PWD}/Prob3plusplus:${PWD}/ProbGPU:

#Environment configs
export UseGPU=1
export UseMultithreading=1
export UseBinned=1
export OMP_NUM_THREADS=8
export UseDoubles=1

if [ ${UseBinned} == 1 ]; then
    echo "Using Binned Probabilities"
else
    echo "Using Unbinned Probabilities"
fi

if [ ${UseGPU} == 1 ]; then
    echo "Using GPU"
else
    echo "Not using GPU"
fi

if [ ${UseMultithreading} == 1 ]; then
    echo "Using Multithreading with:" ${OMP_NUM_THREADS} "threads"
else
    echo "Not using Multithreading"
fi

#Which calculators to compile
export UseCUDAProb3=1
export UseCUDAProb3Linear=0
export UseProbGPULinear=1
export UseProb3ppLinear=1

#ProbGPU only supported when using GPU
if [ ${UseGPU} == 0 ]; then
    export UseProbGPULinear=0
fi

#Set CUDA environment variables
if [ ${UseGPU} == 1 ]; then
    export CUDAPATH=/usr/local/cuda/
    export PATH=${PATH}:${CUDAPATH}/bin/
fi 

if [ ${UseCUDAProb3} == 1 ] && [ ${UseCUDAProb3Linear} == 1 ]; then
    echo "error: CUDAProb3 and CUDAProb3Linear and not able to be built at the same time"
fi

if [ ${UseCUDAProb3Linear} == 1 ]; then
    export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${PWD}/CUDAProb3Linear/build
fi
