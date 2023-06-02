#include "OscillatorBase.h"

#include <iostream>

int main() {
  std::cout << "========================================================" << std::endl;
  std::cout << "Starting setup in executable" << std::endl;

  OscillatorBase* Oscill = new OscillatorBase();
  Oscill->ImplementationName();
}
