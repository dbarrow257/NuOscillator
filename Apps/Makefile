INCS = -I. -I..
LIBS = -L/usr/lib64 -L/usr/lib -L../Lib
CXXFLAGS = -Wall -fPIC -fopenmp -O3 -std=c++11 -g

UseGPUFlag = 
ifeq (${UseGPU},1)
	UseGPUFlag = -DUseGPU=1 -DGPU_ON=1
	CUDA_LIBS = -L$(CUDAPATH)/lib64 -I$(CUDAPATH)/include -lcudart
endif

UseBinnedProbsInExeFlag = 
ifeq (${UseBinned},1)
	UseBinnedProbsInExeFlag = -DUseBinned=1
endif

UseMultithreadingFlag = 
ifeq (${UseMultithreading},1)
	UseMultithreadingFlag = -DUseMultithreading=1
endif

CUDAProb3Lib =
CUDAProb3Inc =
CUDAProb3Flags =
ifeq (${UseCUDAProb3},1)
        CUDAProb3Lib = 
        CUDAProb3Inc = -I../CUDAProb3
        CUDAProb3Flags = -DUseCUDAProb3=1
endif

CUDAProb3LinearLib =
CUDAProb3LinearInc =
CUDAProb3LinearFlags =
ifeq (${UseCUDAProb3Linear},1)
        CUDAProb3LinearLib = -L../CUDAProb3Linear/build -lCUDAProb3Beam
        CUDAProb3LinearInc = -I../CUDAProb3Linear
        CUDAProb3LinearFlags = -DUseCUDAProb3Linear=1
endif

ProbGPULinearLib =
ProbGPULinearInc =
ProbGPULinearFlags =
ifeq (${UseProbGPULinear},1)
	ProbGPULinearLib = -L../ProbGPU -lProbGPU
	ProbGPULinearInc = -I../ProbGPU
	ProbGPULinearFlags = -DUseProbGPULinear=1
endif

Prob3ppLinearLib =
Prob3ppLinearInc =
Prob3ppLinearFlags =
ifeq (${UseProb3ppLinear},1)
	Prob3ppLinearLib = -L../Prob3plusplus -lThreeProb_3.20
	Prob3ppLinearInc = -I../Prob3plusplus
        Prob3ppLinearFlags = -DUseProb3ppLinear=1
endif

TARLIBS = ${CUDAProb3Lib} ${CUDAProb3LinearLib} ${ProbGPULinearLib} ${Prob3ppLinearLib}
TARINCS = ${CUDAProb3Inc} ${CUDAProb3LinearInc} ${ProbGPULinearInc} ${Prob3ppLinearInc}
TARFLAGS = ${CUDAProb3Flags} ${CUDAProb3LinearFlags} ${ProbGPULinearFlags} ${Prob3ppLinearFlags}

FLOAT_TFLAGS = -DUseDoubles=1

OscProbCalcerLib = -lOscProbCalcer
OscProbCalcerInc = -I../OscProbCalcer

OscillatorLib = -lOscillator
OscillatorInc = -I../Oscillator

all: Analysis.exe DragRace.exe

Analysis.exe: Analysis.cpp
	g++ ${CXXFLAGS} ${LIBS} ${TARLIBS} ${CUDA_LIBS} ${OscProbCalcerLib} ${OscillatorLib} ${INCS} ${TARINCS} ${OscillatorInc} ${OscProbCalcerInc} Analysis.cpp -o Analysis.exe ${TARFLAGS} ${FLOAT_TFLAGS} ${UseBinnedProbsInExeFlag}
	mv Analysis.exe ../Bin

DragRace.exe: DragRace.cpp
	g++ ${CXXFLAGS} ${LIBS} ${TARLIBS} ${CUDA_LIBS} ${OscProbCalcerLib} ${OscillatorLib} ${INCS} ${TARINCS} ${OscillatorInc} ${OscProbCalcerInc} DragRace.cpp -o DragRace.exe ${TARFLAGS} ${FLOAT_TFLAGS} ${UseBinnedProbsInExeFlag}
	mv DragRace.exe ../Bin

clean:
	rm -f *.o
	rm -f *.exe
	rm -f *~
