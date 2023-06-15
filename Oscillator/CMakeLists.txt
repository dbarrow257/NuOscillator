SET(OSCILLATOR_VERSION 0.0)
project(OSCILLATOR VERSION ${OSCILLATOR_VERSION} LANGUAGES CXX)

add_library(Oscillator SHARED OscillatorBase.cpp OscillatorUnbinned.cpp OscillatorBinned.cpp)

target_include_directories(Oscillator PUBLIC ../Constants PRIVATE ../OscProbCalcer)

target_link_libraries(Oscillator OscProbCalcer::OscProbCalcer)

set_target_properties(Oscillator PROPERTIES 
				 EXPORT_NAME Oscillator)

add_library(Oscillator::Oscillator ALIAS Oscillator)
