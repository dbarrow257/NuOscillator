add_custom_target(OscApps)

foreach(app
        DragRace
        Analysis
        SingleOscillator
        SingleOscProbCalcer)

        add_executable(${app} ${app}.cpp)
        add_dependencies(OscApps ${app})
        target_link_libraries(${app} NuOscillator::All)
        install(TARGETS ${app} DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)
endforeach(app)
