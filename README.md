# PBPKPredictionVariability
https://doi.org/###/zenodo.####

Modeling of tissue disposition under deep epistemic uncertainty

# Introduction
This repository contains code and data files used in the preparation of the manuscript "<i>Prediction variability in physiologically based pharmacokinetic modeling of tissue disposition under deep uncertainty</i>" 
by Farahat et al. (2025)

The following authors contributed to this repository:

Mustafa Farahat (University of Tennessee)

DT Flaherty (University of Tennessee)

Belinda S Akpa (University of Tennessee and National Institute for Modeling Biological Systems)

Zachary Fox (Oak Ridge National Laboratory)

# Code contents
<ul>
<li>AkpaFarahatKpuModel.m</li>

<li>MathewKpuModel.m</li>

<li>PearceKpuModel.m</li>

<li>MonteCarloModeling.m</li>

<li>GlobalSensitivityAnalysis.m</li>

</ul>

# Data file contents
<ul>
  
<li>Obach862Molecules.csv</li>

<li>Synthetic10000Molecules.csv</li>

<li>tissueCompositionParamsHuman.csv</li>

<li>tissueCompositionParamsRat.csv</li>

<li>tissueCompositionPearceParamsHuman.csv</li>

<li>tissueCompositionPearceParamsRat.csv</li>

</ul>

# Usage
<ul>
  
<li><strong>AkpaFarahatKpuModel.m</strong>

Contains a function AkpaFarahatKpuModel(Params) with the following argument:

<i>Params - parameter struct containing molecule properties required to predict plasma-tissue partition coefficients.</i>

Function reads in tissueCompositionParamsHuman.csv or tissueCompositionParamsRat.csv, depending on the species specified in the Params struct.</li>

<li><strong>MathewKpuModel.m</strong>

Contains a function MathewKpuModel(Params) with the following argument:

<i>Params - parameter struct containing molecule properties required to predict plasma-tissue partition coefficients using the model reported by Mathew et al. (2023).</i>

Function reads in tissueCompositionParamsHuman.csv or tissueCompositionParamsRat.csv, depending on the species specified in the Params struct.</li>

<li><strong>PearceKpuModel.m</strong>

Contains a function PearceKpuModel(Params) with the following argument:

<i>Params - parameter struct containing molecule properties required to predict plasma-tissue partition coefficients using the model reported by Pearce et al. (2017).</i>

Function reads in tissueCompositionPearceParamsHuman.csv or tissueCompositionPearceParamsRat.csv, depending on the species specified in the Params struct. Also reads in PearceRegressionParams.csv, which contains slope and intercept values for the linear regression corrections employed in the calibrated version of the Pearce model.</li>

<li><strong>MonteCarloModeling.m</strong>

Contains a function MonteCarloModeling(Nsamples) with the following argument:

Nsamples - integer specifying the number of realizations to execute for each molecule.</li>

</ul>
