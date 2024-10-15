# NuOscillator

[![Container Image](https://img.shields.io/badge/Container-Image-brightgreen)](https://github.com/dbarrow257/NuOscillator/pkgs/container/nuoscillator)
[![Code - Doxygen](https://img.shields.io/badge/Code-Doxygen-2ea44f)](https://dbarrow257.github.io/NuOscillator/)
[![Build CI](https://github.com/dbarrow257/NuOscillator/actions/workflows/CIBuild.yml/badge.svg)](https://github.com/dbarrow257/NuOscillator/actions/workflows/CIBuild.yml)

## How to start
```bash
mkdir build;
cd build;
cmake ../ -DUseGPU=0 -DUseMultithreading=1 -DUseDoubles=0 -DUseCUDAProb3=0 -DUseCUDAProb3Linear=1 -DUseProb3ppLinear=1 -DUseNuFASTLinear=1 -DUseProbGPULinear=0
make -jN [Where N is number of threads]
make install
```

don't forget about
```bash
source Linux/bin/setup.NuOscillator.sh
```

then you can check if everything runs correctly by
```bash
cd ../
./build/Linux/bin/DragRace
```

## How to Integrate in Framework
Recommended way is to use CPM within you CmakeList.txt
```Cmake
CPMAddPackage(
  NAME NuOscillator
    VERSION 0.0
    GITHUB_REPOSITORY "dbarrow257/NuOscillator"
    GIT_TAG "main"
    OPTIONS
    "UseGPU 0"
    "UseMultithreading 1"
    "UseDoubles 1"

    "UseCUDAProb3Linear 0"
    "UseCUDAProb3 0"
    "UseProbGPULinear 0"
    "UseProb3ppLinear 0"
    "UseNuFASTLinear  1"
)
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
```bash
cmake ../ -DUseGPU=1
```
You can also pass CUDA architectures via cmake to ensure consistent architectures in multiple packages:
```bash
cmake ../ -DCMAKE_CUDA_ARCHITECTURES="all"
```

## Compiler Flags
NuOscillator natively doesn't include many c++ compiler flags (this include optimalisation flags). Since package is intended to be used withing other (fitting) frameworks we allow to pass compiler flags via cmake.
```bash
cmake ../ -DNuOscillator_Compiler_Flags="-O3;-g;-Wextra;-Wall"
```
This give full freedom to users in how to configure NuOscillator.


## How to Use in Fitting Framework
First initialise factory to produce engine defined by config. See exmaples of configs [here](https://github.com/dbarrow257/NuOscillator/tree/main/Configs)
```cpp
OscillatorFactory* OscillFactory = new OscillatorFactory();

std::string NuOscillatorConfigFile = "configs/NuFASTLinear.yaml");
NuOscProbCalcers = OscillFactory->CreateOscillator(NuOscillatorConfigFile);
```

Now you need specify energy bins for which oscillation are going to be calculated.
Energy must be in GeV!
It is good practice to sort energy values as NuOscillator expect them in ascending order.
```cpp
std::vector<double> EnergyArray = {0, 0.6, 1, 10};
std::sort(EnergyArray.begin(),EnergyArray.end());
```

Lastly pass energy vector and finish setup
```cpp
if (!NuOscProbCalcers->EvalPointsSetInConstructor()) {
NuOscProbCalcers->SetEnergyArrayInCalcer(EnergyArray);
}
NuOscProbCalcers->Setup();
```
