INCS = -I. -I.. 
LIBS = -L/usr/lib64 -L/usr/lib
CXXFLAGS = -Wall -fPIC -fcommon -O3 -g -std=c++11

all: Analysis.exe

OscillatorBase.o : OscillatorBase.cpp
	g++ $(CXXFLAGS) ${LIBS} ${INCS} -o OscillatorBase.o -c OscillatorBase.cpp

OscillatorCUDAProb3.o : OscillatorBase.o OscillatorCUDAProb3.cpp
	g++ ${CXXFLAGS} ${LIBS} ${INCS} -o OscillatorCUDAProb3.o -c OscillatorCUDAProb3.cpp

Analysis.exe: OscillatorBase.o OscillatorCUDAProb3.o
	g++ ${CXXFLAGS} ${LIBS} ${INCS} Analysis.cpp OscillatorCUDAProb3.o OscillatorBase.o -o Analysis.exe

clean:
	rm *.o
	rm *.exe
	rm *~
