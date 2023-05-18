#include "OscillatorCUDAProb3.h"
#include "OscillatorBase.h"

int main() {
  OscillatorCUDAProb3* CUDAProb3 = new OscillatorCUDAProb3();
  CUDAProb3->WhoAmI();

  OscillatorBase* Base = (OscillatorBase*)CUDAProb3;
  Base->WhoAmI();
}
