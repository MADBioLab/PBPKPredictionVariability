# PBPKPredictionVariability
https://doi.org/###/zenodo.####

Modeling of tissue disposition under deep epistemic uncertainty

# Introduction
This repository contains code and data files used in the preparation of the manuscript "Prediction variability in physiologically based pharmacokinetic modeling of tissue disposition under deep uncertainty" 
by Farahat et al. (2025)

The following authors contributed to this repository:

Mustafa Farahat (University of Tennessee)

DT Flaherty (University of Tennessee)

Belinda S Akpa (University of Tennessee and National Institute for Modeling Biological Systems)

Zachary Fox (Oak Ridge National Laboratory)

# Code contents
AkpaFarahatKpuModel.m

MathewKpuModel.m

PearceKpuModelUncalibrated.m

PearceKpuModelCalibrated.m

MonteCarloModeling.m

GlobalSensitivityAnalysis.m

# Data file contents
Obach862Molecules.csv

Synthetic10000Molecules.csv

tissueCompositionParamsHuman.csv

tissueCompositionParamsRat.csv

tissueCompositionPearceParamsHuman.csv

tissueCompositionPearceParamsRat.csv

# Usage
AkpaFarahatKpuModel.m
Contains a function AkpaFarahatKpuModel(Params) with the following argument:

Params - parameter struct containing molecule properties required to predict plasma-tissue partition coefficients. 

Function reads in tissueCompositionParamsHuman.csv or tissueCompositionParamsRat.csv, depending on the species specified in the Params struct.

MathewKpuModel.m
Contains a function MathewKpuModel(Params) with the following argument:

Params - parameter struct containing molecule properties required to predict plasma-tissue partition coefficients. 

Function reads in tissueCompositionParamsHuman.csv or tissueCompositionParamsRat.csv, depending on the species specified in the Params struct

PearceKpuModel.m
Contains a function AkpaFarahatKpuModel(Params) with the following argument:

Params - parameter struct containing molecule properties required to predict plasma-tissue partition coefficients. 

Function reads in tissueCompositionPearceParamsHuman.csv or tissueCompositionPearceParamsRat.csv, depending on the species specified in the Params struct. Also reads in PearceRegressionParams.csv, which contains slope and intercept values for the linear regression corrections employed in the calibrated version of the Pearce model.

MonteCarloModeling.m
Contains a function MonteCarloModeling(Nsamples) with the following argument:

Nsamples - integer specifying the number of realizations to execute for each molecule.
