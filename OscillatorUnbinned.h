#ifndef __OSCILLATOR_UNBINNED_BASE_H__
#define __OSCILLATOR_UNBINNED_BASE_H__

#include "OscillatorBase.h"

class OscillatorUnbinned : public OscillatorBase {
 public:
  OscillatorUnbinned(std::vector<std::string> OscProbCalcerImplementationToCreate);

  // ========================================================================================================================================================================
  // Public functions which are calculation implementation agnostic
  
  // ========================================================================================================================================================================
  // Public virtual functions which need calculater specific implementations

 protected:

  // ========================================================================================================================================================================
  // Protected functions which are calculation implementation agnostic  

  // ========================================================================================================================================================================
  // Protected virtual functions which are calculation implementation agnostic

  // ========================================================================================================================================================================
  // Basic protected variables required for oscillation probability calculation

 private:

  // ========================================================================================================================================================================
  // Basic private variables required for oscillation probability calculation

};

#endif
