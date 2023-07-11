SET(APPS_VERSION 0.0)
project(Apps VERSION ${APPS_VERSION} LANGUAGES CXX)

foreach(app DragRace Analysis SingleOscillator SingleOscProbCalcer)
	    add_executable(${app} ${app}.cpp)
	    target_link_libraries(${app} Oscillator OscProbCalcer)
	    target_include_directories(${app} PRIVATE ../Oscillator ../OscProbCalcer)
	    install(TARGETS ${app} DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)
endforeach(app)
