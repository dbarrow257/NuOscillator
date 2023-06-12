INCS = -I. -I..
LIBS = -L/usr/lib64 -L/usr/lib 
CXXFLAGS = -Wall -fPIC -fopenmp -O3 -std=c++11 -g

OscProbCalcerInc = -I../OscProbCalcer
OscProbCalcerLib = -L../Lib -lOscProbCalcer

UseGPUFlag = 
ARCH = 
CUDA_LIBS =
ifeq (${UseGPU},1)
	ARCH= 	-gencode arch=compute_52,code=sm_52 \
		-gencode arch=compute_60,code=sm_60 \
		-gencode arch=compute_61,code=sm_61 \
		-gencode arch=compute_70,code=sm_70 \
		-gencode arch=compute_75,code=sm_75 \
		-gencode arch=compute_80,code=sm_80 \
		-gencode arch=compute_86,code=sm_86 \
		-gencode arch=compute_86,code=compute_86
	UseGPUFlag = -DUseGPU=1 -DGPU_ON=1
	CUDA_LIBS = -L$(CUDAPATH)/lib64 -I$(CUDAPATH)/include -lcudart
endif

UseMultithreadingFlag = 
ifeq (${UseMultithreading},1)
	UseMultithreadingFlag = -DUseMultithreading=1
endif

CUDAProb3Lib =
CUDAProb3Inc =
CUDAProb3Obj =
CUDAProb3Flags =
ifeq (${UseCUDAProb3},1)
        CUDAProb3Lib = 
        CUDAProb3Inc = -I../CUDAProb3
#        CUDAProb3Obj = OscProbCalcer_CUDAProb3.o
        CUDAProb3Flags = -DUseCUDAProb3=1
endif

CUDAProb3LinearLib =
CUDAProb3LinearInc =
CUDAProb3LinearObj =
CUDAProb3LinearFlags =
ifeq (${UseCUDAProb3Linear},1)
        CUDAProb3LinearLib = -L../CUDAProb3Linear/build -lCUDAProb3Beam
        CUDAProb3LinearInc = -I../CUDAProb3Linear
        CUDAProb3LinearObj = OscProbCalcer_CUDAProb3Linear.o
        CUDAProb3LinearFlags = -DUseCUDAProb3Linear=1
endif

ProbGPULinearLib =
ProbGPULinearInc =
ProbGPULinearObj =
ProbGPULinearFlags =
ifeq (${UseProbGPULinear},1)
	ProbGPULinearLib = -L../ProbGPU -lProbGPU
	ProbGPULinearInc = -I../ProbGPU
	ProbGPULinearObj = OscProbCalcer_ProbGPULinear.o
	ProbGPULinearFlags = -DUseProbGPULinear=1
endif

Prob3ppLinearLib =
Prob3ppLinearInc =
Prob3ppLinearObj = 
Prob3ppLinearFlags =
ifeq (${UseProb3ppLinear},1)
	Prob3ppLinearLib = -L../Prob3plusplus -lThreeProb_3.20
	Prob3ppLinearInc = -I../Prob3plusplus
        Prob3ppLinearObj = OscProbCalcer_Prob3ppLinear.o
        Prob3ppLinearFlags = -DUseProb3ppLinear=1
endif

TARLIBS = ${CUDAProb3Lib} ${CUDAProb3LinearLib} ${ProbGPULinearLib} ${Prob3ppLinearLib}
TARINCS = ${CUDAProb3Inc} ${CUDAProb3LinearInc} ${ProbGPULinearInc} ${Prob3ppLinearInc}
TAROBJS = ${CUDAProb3Obj} ${CUDAProb3LinearObj} ${ProbGPULinearObj} ${Prob3ppLinearObj}
TARFLAGS = ${CUDAProb3Flags} ${CUDAProb3LinearFlags} ${ProbGPULinearFlags} ${Prob3ppLinearFlags}

FLOAT_TFLAGS = -DUseDoubles=1

OSCOBJS = OscillatorUnbinned.o OscillatorBinned.o

all: libOscillator.so

OscillatorBase.o : OscillatorBase.cpp
	g++ ${CXXFLAGS} ${LIBS} ${TARLIBS} ${INCS} ${TARINCS} ${OscProbCalcerInc} ${OscProbCalcerLib} -o OscillatorBase.o -c OscillatorBase.cpp ${TARFLAGS} ${FLOAT_TFLAGS}

OscillatorUnbinned.o : OscillatorBase.o OscillatorUnbinned.cpp
	g++ ${CXXFLAGS} ${LIBS} ${TARLIBS} ${INCS} ${TARINCS} ${OscProbCalcerInc} ${OscProbCalcerLib} -o OscillatorUnbinned.o -c OscillatorUnbinned.cpp ${TARFLAGS} ${FLOAT_TFLAGS}

OscillatorBinned.o : OscillatorBase.o OscillatorBinned.cpp
	g++ ${CXXFLAGS} ${LIBS} ${TARLIBS} ${INCS} ${TARINCS} ${OscProbCalcerInc} ${OscProbCalcerLib} -o OscillatorBinned.o -c OscillatorBinned.cpp ${TARFLAGS} ${FLOAT_TFLAGS}

libOscillator.a: OscillatorBase.o ${OSCOBJS}
	ar rcsv libOscillator.a OscillatorBase.o ${OSCOBJS}

libOscillator.so: libOscillator.a OscillatorBase.o ${OSCOBJS}
	g++ -shared OscillatorBase.o ${OSCOBJS} -o libOscillator.so ${LIBS} ${OscProbCalcerLib}
	mv libOscillator.so ../Lib/.

clean:
	rm -f *.o
	rm -f *~
	rm -f ../Lib/libOscillator.so
