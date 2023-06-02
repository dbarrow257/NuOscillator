#ifndef __OSCILLATOR_BASE_H__
#define __OSCILLATOR_BASE_H__

#include "OscProbCalcerBase.h"

class OscillatorBase {
 public:
  OscillatorBase();

  // ========================================================================================================================================================================
  // Public functions which are calculation implementation agnostic
  
  void ImplementationName();

  // ========================================================================================================================================================================
  // Public virtual functions which need calculater specific implementations

 protected:

  // ========================================================================================================================================================================
  // Protected functions which are calculation implementation agnostic  

  void InitialiseOscProbCalcer();

  // ========================================================================================================================================================================
  // Protected virtual functions which are calculation implementation agnostic

  // ========================================================================================================================================================================
  // Basic protected variables required for oscillation probability calculation

  OscProbCalcerBase* OPCalcer;

 private:
  // ========================================================================================================================================================================
  // Basic private variables required for oscillation probability calculation

};

#endif
