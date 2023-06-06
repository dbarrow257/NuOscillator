INCS = -I.
LIBS = -L/usr/lib64 -L/usr/lib
CXXFLAGS = -Wall -fPIC -fopenmp -O3 -std=c++11 -g

UseGPUFlag = 
ARCH = 
CUDA_LIBS =
ifdef UseGPU
	ARCH= 	-gencode arch=compute_52,code=sm_52 \
		-gencode arch=compute_60,code=sm_60 \
		-gencode arch=compute_61,code=sm_61 \
		-gencode arch=compute_70,code=sm_70 \
		-gencode arch=compute_75,code=sm_75 \
		-gencode arch=compute_80,code=sm_80 \
		-gencode arch=compute_86,code=sm_86 \
		-gencode arch=compute_86,code=compute_86
	UseGPUFlag = -DUseGPU=1
	CUDA_LIBS = -L$(CUDAPATH)/lib64 -I$(CUDAPATH)/include -lcudart
endif

MultiThreadFlags = 
ifdef MULTITHREAD
	MultiThreadFlags = -DUseMultithread=1
endif

CUDAProb3Lib =
CUDAProb3Inc =
CUDAProb3Obj =
CUDAProb3Flags =
ifdef UseCUDAProb3
        CUDAProb3Lib = 
        CUDAProb3Inc = -I../CUDAProb3
        CUDAProb3Obj = OscProbCalcer_CUDAProb3.o
        CUDAProb3Flags = -DUseCUDAProb3=1
endif

ProbGPULinearLib =
ProbGPULinearInc =
ProbGPULinearObj =
ProbGPULinearFlags =
ifdef UseProbGPULinear
	ProbGPULinearLib = -L../ProbGPU -lProbGPU
	ProbGPULinearInc = -I../ProbGPU
	ProbGPULinearObj = OscProbCalcer_ProbGPULinear.o
	ProbGPULinearFlags = -DUseProbGPULinear=1
endif

Prob3ppLinearLib =
Prob3ppLinearInc =
Prob3ppLinearObj = 
Prob3ppLinearFlags =
ifdef UseProb3ppLinear
	Prob3ppLinearLib = -L../Prob3 -lThreeProb_3.20
	Prob3ppLinearInc = -I../Prob3
        Prob3ppLinearObj = OscProbCalcer_Prob3ppLinear.o
        Prob3ppLinearFlags = -DUseProb3ppLinear=1
endif

TARLIBS = ${CUDAProb3Lib} ${ProbGPULinearLib} ${Prob3ppLinearLib}
TARINCS = ${CUDAProb3Inc} ${ProbGPULinearInc} ${Prob3ppLinearInc}
TAROBJS = ${CUDAProb3Obj} ${ProbGPULinearObj} ${Prob3ppLinearObj}
TARFLAGS = ${CUDAProb3Flags} ${ProbGPULinearFlags} ${Prob3ppLinearFlags}

FLOAT_TFLAGS = -DUsingDoubles=1

all: Analysis.exe

OscProbCalcerBase.o : 
	g++ $(CXXFLAGS) ${LIBS} ${INCS} -o OscProbCalcerBase.o -c OscProbCalcerBase.cpp ${FLOAT_TFLAGS}

ifdef CUDAPATH
OscProbCalcer_CUDAProb3.o : OscProbCalcerBase.o
	nvcc -g -O3 -x cu $(ARCH) -lineinfo -std=c++11 -Xcompiler="${CUDA_LIBS} ${INCS} ${CUDAProb3Inc} $(CXXFLAGS) ${CUDAProb3Flags} ${FLOAT_TFLAGS} ${UseGPUFlag}" -c OscProbCalcer_CUDAProb3.cpp -o OscProbCalcer_CUDAProb3.o
else
OscProbCalcer_CUDAProb3.o : OscProbCalcerBase.o
	g++ ${CXXFLAGS} ${LIBS} ${CUDAProb3Lib} ${INCS} ${CUDAProb3Inc} -o OscProbCalcer_CUDAProb3.o -c OscProbCalcer_CUDAProb3.cpp ${CUDAProb3Flags} ${FLOAT_TFLAGS} ${MultiThreadFlags} ${UseGPUFlag}
endif

OscProbCalcer_ProbGPULinear.o : OscProbCalcerBase.o
	g++ ${CXXFLAGS} ${LIBS} ${ProbGPULinearLib} ${INCS} ${ProbGPULinearInc} -o OscProbCalcer_ProbGPULinear.o -c OscProbCalcer_ProbGPULinear.cpp ${ProbGPULinearFlags} ${FLOAT_TFLAGS}

OscProbCalcer_Prob3ppLinear.o : OscProbCalcerBase.o
	g++ ${CXXFLAGS} ${LIBS} ${Prob3ppLinearLib} ${INCS} ${Prob3ppLinearInc} -o OscProbCalcer_Prob3ppLinear.o -c OscProbCalcer_Prob3ppLinear.cpp ${Prob3ppLinearFlags} ${FLOAT_TFLAGS}

OscillatorBase.o : ${TAROBJS}
	g++ ${CXXFLAGS} ${LIBS} ${TARLIBS} ${INCS} ${TARINCS} -o OscillatorBase.o -c OscillatorBase.cpp ${TARFLAGS} ${FLOAT_TFLAGS}

OscillatorUnbinned.o : OscillatorBase.o
	g++ ${CXXFLAGS} ${LIBS} ${TARLIBS} ${INCS} ${TARINCS} -o OscillatorUnbinned.o -c OscillatorUnbinned.cpp ${TARFLAGS} ${FLOAT_TFLAGS}

OscillatorBinned.o : OscillatorBase.o
	g++ ${CXXFLAGS} ${LIBS} ${TARLIBS} ${INCS} ${TARINCS} -o OscillatorBinned.o -c OscillatorBinned.cpp ${TARFLAGS} ${FLOAT_TFLAGS}

Analysis.exe: OscillatorUnbinned.o OscillatorBinned.o
	g++ ${CXXFLAGS} ${LIBS} ${TARLIBS} ${CUDA_LIBS} ${INCS} ${TARINCS} Analysis.cpp OscProbCalcerBase.o ${TAROBJS} OscillatorBase.o OscillatorUnbinned.o OscillatorBinned.o -o Analysis.exe ${TARFLAGS} ${FLOAT_TFLAGS}

clean:
	rm -f *.o
	rm -f *.exe
	rm -f *~
