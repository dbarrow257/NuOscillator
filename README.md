# NuOscillator

[![Container Image](https://img.shields.io/badge/Container-Image-brightgreen)](https://github.com/dbarrow257/NuOscillator/pkgs/container/nuoscillator)
[![Code - Doxygen](https://img.shields.io/badge/Code-Doxygen-2ea44f)](https://dbarrow257.github.io/NuOscillator/)
[![Build CI](https://github.com/mach3-software/MaCh3/actions/workflows/CIBuild.yml/badge.svg)](https://github.com/dbarrow257/NuOscillator/actions/workflows/CIBuild.yml)

## How to start
```
mkdir build;
cd build;
cmake ../ -DUseGPU=0 -DUseMultithreading=1 -DUseDoubles=0 -DUseCUDAProb3=0 -DUseCUDAProb3Linear=1 -DUseProb3ppLinear=1 -DUseNuFASTLinear=1 -DUseProbGPULinear=0
make -jN [Where N is number of threads]
make install
```

## Implemented Engines
`UseCUDAProb3` etc. refers to implemented engines. Engines are loaded via yaml config files. In principle you can compile all of them and select one you want to use via config.

Following neutrino oscillation calculators are available:
|Oscillator        | Hardware   | Source     | Reference  |
|------------------|------------|------------|------------|
| CUDAProb3Linear  | CPU/GPU    | Beam       |            |
| CUDAProb3        | CPU/GPU    | Atm        | [Ref](https://doi.org/10.1016/j.cpc.2018.07.022)        |
| ProbGPULinear    | GPU        | Beam       | [Ref](http://dx.doi.org/10.3204/DESY-PROC-2014-05/23)   |
| Prob3++Linear    | CPU        | Beam       |            |
| NuFastLinear     | CPU        | Beam       | [Ref](https://doi.org/10.48550/arXiv.2405.02400)        |

## GPU
Some engines requires gpu like `ProbGPULinear` other can use both CPU and GPU. To use GPU functionality remember about
```
cmake ../ -DUseGPU=1
```

## Other

If the output from running cmake appears to have picked up the wrong ROOT install, try manually adding the path to the FindROOT.cmake file in your ROOT install to CMAKE_MODULE_PATH like so -DCMAKE_MODULE_PATH=${ROOTSYS}/etc/cmake. This shouldn't be a problem for more recent versions of ROOT6 where the ROOT developers have taken care to expose modern, conventional CMake.
