INCS = -I. -I..
LIBS = -L/usr/lib64 -L/usr/lib
CXXFLAGS = -Wall -fPIC -fopenmp -O3 -std=c++11 -g

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
        CUDAProb3Obj = OscProbCalcer_CUDAProb3.o
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

all: libOscProbCalcer.so

OscProbCalcerBase.o : OscProbCalcerBase.cpp
	g++ $(CXXFLAGS) ${LIBS} ${INCS} -o OscProbCalcerBase.o -c OscProbCalcerBase.cpp ${FLOAT_TFLAGS}

ifeq (${UseGPU},1)
OscProbCalcer_CUDAProb3.o : OscProbCalcerBase.o OscProbCalcer_CUDAProb3.cpp
	nvcc -g -O0 -x cu $(ARCH) -lineinfo -std=c++11 -Xcompiler="${INCS} ${CUDAProb3Inc} $(CXXFLAGS) ${CUDAProb3Flags} ${FLOAT_TFLAGS} ${UseGPUFlag}" -c OscProbCalcer_CUDAProb3.cpp -o OscProbCalcer_CUDAProb3.o
else
OscProbCalcer_CUDAProb3.o : OscProbCalcerBase.o OscProbCalcer_CUDAProb3.cpp
	g++ ${CXXFLAGS} ${LIBS} ${CUDAProb3Lib} ${INCS} ${CUDAProb3Inc} -o OscProbCalcer_CUDAProb3.o -c OscProbCalcer_CUDAProb3.cpp ${CUDAProb3Flags} ${FLOAT_TFLAGS} ${UseMultithreadingFlag} ${UseGPUFlag}
endif

ifeq (${UseGPU},1)
OscProbCalcer_CUDAProb3Linear.o : OscProbCalcerBase.o OscProbCalcer_CUDAProb3Linear.cpp
	nvcc -g -O0 -x cu $(ARCH) -lineinfo -std=c++11 -Xcompiler="${INCS} ${CUDAProb3LinearInc} $(CXXFLAGS) ${CUDAProb3LinearFlags} ${FLOAT_TFLAGS} ${UseGPUFlag}" -c OscProbCalcer_CUDAProb3Linear.cpp -o OscProbCalcer_CUDAProb3Linear.o
else
OscProbCalcer_CUDAProb3Linear.o : OscProbCalcerBase.o OscProbCalcer_CUDAProb3Linear.cpp
	g++ ${CXXFLAGS} ${LIBS} ${CUDAProb3LinearLib} ${INCS} ${CUDAProb3LinearInc} -o OscProbCalcer_CUDAProb3Linear.o -c OscProbCalcer_CUDAProb3Linear.cpp ${CUDAProb3LinearFlags} ${FLOAT_TFLAGS} ${UseMultithreadingFlag} ${UseGPUFlag}
endif

OscProbCalcer_ProbGPULinear.o : OscProbCalcerBase.o OscProbCalcer_ProbGPULinear.cpp
	g++ ${CXXFLAGS} ${LIBS} ${ProbGPULinearLib} ${INCS} ${ProbGPULinearInc} -o OscProbCalcer_ProbGPULinear.o -c OscProbCalcer_ProbGPULinear.cpp ${ProbGPULinearFlags} ${FLOAT_TFLAGS}

OscProbCalcer_Prob3ppLinear.o : OscProbCalcerBase.o OscProbCalcer_Prob3ppLinear.cpp
	g++ ${CXXFLAGS} ${LIBS} ${Prob3ppLinearLib} ${INCS} ${Prob3ppLinearInc} -o OscProbCalcer_Prob3ppLinear.o -c OscProbCalcer_Prob3ppLinear.cpp ${Prob3ppLinearFlags} ${FLOAT_TFLAGS}

libOscProbCalcer.a: OscProbCalcerBase.o  ${TAROBJS}
	ar rcsv libOscProbCalcer.a OscProbCalcerBase.o  ${TAROBJS}

libOscProbCalcer.so: libOscProbCalcer.a OscProbCalcerBase.o ${TAROBJS}
	g++ -shared OscProbCalcerBase.o ${TAROBJS} -o libOscProbCalcer.so ${LIBS} ${TARLIBS}
	mv libOscProbCalcer.so ../Lib/.

clean:
	rm -f *.o
	rm -f *.a
	rm -f *~
	rm -f ../Lib/libOscProbCalcer.so
