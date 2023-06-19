#ifndef __OSCILLATOR_CONSTANTS__
#define __OSCILLATOR_CONSTANTS__

#include <string>
#include <iostream>

/**
 * @file OscillatorConstants.h
 */

/**
 * @brief Different verbosity levels for console output
 */
enum Verbosity{NONE,INFO};

//DB
inline int Verbosity_StrToInt(std::string Verbosity) {
  if (Verbosity == "NONE") {
    return NONE;
  } else if (Verbosity == "INFO") {
    return INFO;
  } else {
    std::cerr << "Invalid verbosity provided:" << Verbosity << std::endl;
    throw;
  }
  
  return -1;
}
#endif
