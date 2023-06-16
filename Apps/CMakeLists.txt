SET(APPS_VERSION 0.0)
project(Apps VERSION ${APPS_VERSION} LANGUAGES CXX)

foreach(app DragRace Analysis)
	    add_executable(${app} ${app}.cpp)
	    target_link_libraries(${app} Oscillator OscProbCalcer)
	    target_include_directories(${app} PRIVATE ../Oscillator ../OscProbCalcer)
	    install(TARGETS ${ap} DESTINATION .)
endforeach(app)
