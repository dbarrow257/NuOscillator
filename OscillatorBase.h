#ifndef __OSCILLATOR_BASE_H__
#define __OSCILLATOR_BASE_H__

using FLOAT_T = double;

class OscillatorBase {
 public:
  OscillatorBase();
  void WhoAmI();

  int nEnergyPoints;
  int nCosZPoints;
};

#endif
