# PBPKPredictionVariability
https://doi.org/###/zenodo.####

Modeling of tissue disposition under deep epistemic uncertainty

# Introduction
This repository contains code and data files used in the preparation of the manuscript "<i>Prediction variability in physiologically based pharmacokinetic modeling of tissue disposition under deep uncertainty</i>" 
by Farahat et al. (2025)
https://doi.org/10.64898%2F2025.12.05.692437

The following authors contributed to this repository:

Mustafa Farahat (University of Tennessee)

DT Flaherty (University of Tennessee)

Zachary Fox (Oak Ridge National Laboratory)

Belinda S Akpa (University of Tennessee and National Institute for Modeling Biological Systems)


# Code contents
<ul>
<li>AkpaFarahatKpuModel.m</li>

<li>MathewKpuModel.m</li>

<li>PearceKpuModel.m</li>

<li>PrepareParameterStructure.m</li>

<li>PBPKodes.m</li>

<li>SimulatePBPK.m</li>

<li>MonteCarloSimulations.m</li>

<li>SobolForCorrelatedParamsPBPK.m</li>

</ul>

# Data file contents
<ul>
  
<li>Obach862Molecules.csv</li>

<li>Synthetic10000Molecules.csv</li>

<li>tissueCompositionParamsHuman.csv</li>

<li>tissueCompositionParamsRat.csv</li>

<li>tissueCompositionPearceParamsHuman.csv</li>

<li>tissueCompositionPearceParamsRat.csv</li>

<li>PhysiologicalParameters.csv</li>

<li>PhysiologicalParametersRat.csv</li>

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

<li><strong>PrepareParameterStructure.m</strong>

Contains a function PrepareParameterStructure(Params) with the following argument:

<i>Params - parameter struct </i> </li>

<li><strong>PBPKodes.m</strong>

Contains a function PBPKodes(time, m, P, Y, bodyWeight) with the following arguments:

<i>time - current simulation time [s]

m - array containing total mass of drug in each model compartment [mg]

P - parameter struct containing model parameters (drug and system-specific)

Y - struct containing organ subcompartment concentrations for all organs

bodyWeight - total bodyweight [kg] </i> </li>

<li><strong>SimulatePBPK.m</strong>

Contains a function SimulatePBPK(DrugP, Tmax) with the following arguments:

<i>DrugP - struct containing drug-specific parameters 

Tmax - desired simulation timespan [hrs] </i> </li>

<li><strong>MonteCarloSimulations.m</strong>

Contains a function MonteCarloSimulations(referenceProperties, selectedTissue, Nsamples, errorScale) with the following arguments:

<i> referenceProperties - drug-specific properties for molecules of interest

selectedTissue - numerical ID of tissue of interest

Nsamples - integer specifying the number of realizations to execute for each molecule 

errorScale - scaling factor for mean average errors in parameter distributions </i> </li>

<li><strong>SobolForCorrelatedParamsPBPK.m</strong>

Contains a function SobolForCorrelatedParamsPBPK(ReferenceData) with the following argument:

<i> ReferenceData - array of molecule properties representating the chemical space of interest. NReference rows (number of molecules) and 6 columns. 
    Columns should be (1) molar mass, (2) donor pKa, (3) acceptor pKa, (4) fraction unbound, (5) logCLint, (6) logD74. </i> </li>

</ul>
