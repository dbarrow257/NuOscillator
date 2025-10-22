#!/bin/bash
# This script is supposed to make process of remaking validation inputs easier and more automatic

# Set the base directories
if [ -z "$NUOSCILLATOR_ROOT" ]; then
    echo "Error: NUOSCILLATOR_ROOT is not set"
    exit 1
fi

NUOSC_DIR=${NUOSCILLATOR_ROOT}
OUTPUT_DIR="$NUOSC_DIR/.github/TestOutputs"

# Make sure output directory exists
mkdir -p "$OUTPUT_DIR"

# Define tests in a similar way to your workflow
tests=(
  #"SingleOscProbCalcer NuOscillatorConfigs/OscProbCalcerConfigs/NuFASTLinear.yaml NuFASTLinear_Stored.txt"
  #"SingleOscProbCalcer NuOscillatorConfigs/OscProbCalcerConfigs/OscProbLinear.yaml OscProbLinear_Stored.txt"
  #"SingleOscProbCalcer NuOscillatorConfigs/OscProbCalcerConfigs/Prob3ppLinear.yaml Prob3ppLinear_Stored.txt"
  #"SingleOscProbCalcer NuOscillatorConfigs/OscProbCalcerConfigs/CUDAProb3Linear.yaml CUDAProb3Linear_Stored.txt"
  #"SingleOscProbCalcer NuOscillatorConfigs/OscProbCalcerConfigs/NuSQUIDSLinear.yaml NuSQUIDSLinear_Stored.txt"
  "SingleOscProbCalcer NuOscillatorConfigs/Binned_NuSQUIDSLinear_NSI.yaml NuSQUIDSLinear_UnbinnedOscillator_NSI_Stored.txt"
  #"SingleOscProbCalcer NuOscillatorConfigs/OscProbCalcerConfigs/CUDAProb3.yaml CUDAProb3Atm_Stored.txt"
  #"SingleOscProbCalcer NuOscillatorConfigs/OscProbCalcerConfigs/OscProb.yaml OscProbAtm_Stored.txt"
  #"SingleOscillator NuOscillatorConfigs/Binned_OscProb.yaml OscProb_BinnedOscillator_Stored.txt"
  #"SingleOscillator NuOscillatorConfigs/Unbinned_OscProb.yaml OscProb_UnbinnedOscillator_Stored.txt"
  "SingleOscillator NuOscillatorConfigs/Binned_OscProb.yaml OscProb_UnbinnedOscillator_NSI_Stored.txt"
  #"SingleOscProbCalcer NuOscillatorConfigs/OscProbCalcerConfigs/GLoBESLinear.yaml GLoBESLinear_BinnedOscillator_Stored.txt"
)

# Loop over tests
for test in "${tests[@]}"; do
  set -- $test
  exec_name=$1
  config_file=$2
  output_file=$3

  echo "Generating $output_file..."
  cd $NUOSC_DIR
  "$exec_name" "$config_file" | tee "$OUTPUT_DIR/$output_file"
done

echo "All old outputs regenerated in $OUTPUT_DIR"
