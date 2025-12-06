# PBPKPredictionVariability
https://doi.org/###/zenodo.####

Modeling of tissue disposition under deep epistemic uncertainty

# Introduction
This repository contains code and data files used in the preparation of the manuscript "Prediction variability in physiologically based pharmacokinetic modeling of tissue disposition unde3r deep uncertainty" 
by Farahat et al. (2025)

The following authors contributed to this repository:

Mustafa Farahat (University of Tennessee)
DT Flaherty (University of Tennessee)
Belinda S Akpa (University of Tennessee and National Institute for Modeling Biological Systems)
Zachary Fox (Oak Ridge National Laboratory)

# Code contents
FarahatAkpaKpuModel.m
MathewKpuModel.m
PearceKpuModelUncalibrated.m
PearceKpuModelCalibrated.m
MonteCarloModeling.m
GlobalSensitivityAnalysis.m

# Data file contents
Obach862Molecules.csv
Synthetic10000Molecules.csv

# Usage
*KpuModel.m
Each code file, *KpuModel.m, ...

MonteCarloModeling.m
Contains a function MonteCarloModeling(Nsamples) with the following argument:

Nsamples - integer specifying the number of realizations to execute for each molecule.
