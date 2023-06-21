SET(OSCILLATOR_VERSION 0.0)
project(OSCILLATOR VERSION ${OSCILLATOR_VERSION} LANGUAGES CXX)

set(HEADERS OscillatorBase.h OscillatorUnbinned.h OscillatorBinned.h)
add_library(Oscillator SHARED OscillatorBase.cpp OscillatorUnbinned.cpp OscillatorBinned.cpp)

target_include_directories(Oscillator PUBLIC ../Constants PRIVATE ../OscProbCalcer)

target_link_libraries(Oscillator OscProbCalcer yaml-cpp ROOT::ROOT)

set_target_properties(Oscillator PROPERTIES 
				 EXPORT_NAME Oscillator)

install(TARGETS Oscillator
		EXPORT Oscillator-targets
                LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib/
                PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_PREFIX}/include/Oscillator)

foreach(head ${HEADERS})
             install(FILES ${head} DESTINATION ${CMAKE_INSTALL_PREFIX}/include/Oscillator)
endforeach(head)

add_library(Oscillator::Oscillator ALIAS Oscillator)
