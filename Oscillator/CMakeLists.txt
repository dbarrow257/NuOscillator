set(HEADERS
        OscillatorBase.h
        OscillatorUnbinned.h
        OscillatorBinned.h
        OscillatorFactory.h)

add_library(Oscillator SHARED
        OscillatorBase.cpp
        OscillatorUnbinned.cpp
        OscillatorBinned.cpp
        OscillatorFactory.cpp)


target_include_directories(Oscillator PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/../>
    $<INSTALL_INTERFACE:include>
)

target_link_libraries(Oscillator OscProbCalcer yaml-cpp ROOT::ROOT NuOscillatorCompilerOptions)

set_target_properties(Oscillator PROPERTIES
    PUBLIC_HEADER "${HEADERS}"
    EXPORT_NAME Oscillator)

install(TARGETS Oscillator
        EXPORT NuOscillator-targets
        LIBRARY DESTINATION lib/
        PUBLIC_HEADER DESTINATION include/Oscillator)

add_library(Oscillator::Oscillator ALIAS Oscillator)
