# PBPKPredictionVariability


Code + data for: Farahat et al. (2025), *Prediction variability in physiologically based pharmacokinetic modeling of tissue disposition under deep epistemic uncertainty*  
Preprint DOI: 10.64898/2025.12.05.692437

Zenodo archive: TBD (will be added after release)
---

## Authors

- Mustafa Farahat (University of Tennessee)
- DT Flaherty (University of Tennessee)
- Zachary Fox (Oak Ridge National Laboratory)
- Belinda S Akpa (University of Tennessee; National Institute for Modeling Biological Systems)

---

## What this repository does

This repository provides MATLAB code that:

1. Computes tissue partition terms (`Kpu`) with multiple tissue-partition models.
2. Runs a whole-body PBPK simulation (ODE system) that uses those partition terms.
3. Propagates molecular-property uncertainty via Monte Carlo (`MonteCarloSimulations.m`).
4. Runs variance-based sensitivity analysis with correlated inputs (`SobolForCorrelatedParamsPBPK.m`).

---

## Requirements

- MATLAB R2023a (uses `ode15s`).
- Parallel Computing Toolbox is optional.
  - `MonteCarloSimulations.m` uses `parfor`.
  - Without the toolbox, replace `parfor` with `for`.

---

## Quickstart (Monte Carlo for all molecules in a CSV)

### Step 0 — Set MATLAB path

Set MATLAB Current Folder to the repo root (the folder that contains the `.m` and `.csv` files), then:

```matlab
addpath(genpath(pwd));
```

### Step 1 — Build `referenceProperties` (N × 6) from a CSV

`MonteCarloSimulations` expects a numeric matrix with **6 columns** in this order:

1. `RMM`
2. `pKa_donor`
3. `pKa_acceptor`
4. `fu`
5. `log10(human inferred CLh_int [L/s])` (converted to `CLint` via `10.^...`)
6. `logD74`

#### Obach (862 molecules)

```matlab
T = readtable('Obach862Molecules.csv','VariableNamingRule','preserve');

referenceProperties = [ ...
    T.RMM, ...
    T.pKa_donor, ...
    T.pKa_acceptor, ...
    T.fu, ...
    T.("log10(human inferred CLh_int [L/s])"), ...
    T.logD74 ...
];
```
#### Synthetic (10,000 molecules)

```matlab
T = readtable('Synthetic10000Molecules.csv','VariableNamingRule','preserve');

referenceProperties = [ ...
    T.RMM, ...
    T.pKa_donor, ...
    T.pKa_acceptor, ...
    T.fu, ...
    T.("log10(human inferred CLh_int [L/s])"), ...
    T.logD74 ...
];
```

### Step 2 — Run Monte Carlo

```matlab
selectedTissue = 1;   % numeric tissue index (see "Tissue index" section below)
Nsamples      = 1000; % realizations per molecule
errorScale    = 1;    % scales MAE values in parameter distributions

MonteCarloSimulations(referenceProperties, selectedTissue, Nsamples, errorScale);
```

---

## How the code pieces connect

At a high level:

- `MonteCarloSimulations` samples uncertain inputs and calls `SimulatePBPK` many times.
- `SimulatePBPK` builds the full PBPK parameter struct via `PrepareParameterStructure`, then solves the ODE system via `ode15s` using `PBPKodes`.

Inside `MonteCarloSimulations`, each Monte Carlo draw runs **four** tissue-partition models:

1. `AkpaFarahatKpuModel`
2. `MathewKpuModel`
3. `PearceKpuModel` (calibrated Pearce; reads regression parameters)
4. `PearceKpuModelNoCalibration` (uncalibrated Pearce; does not use regression corrections)

---

## Outputs from Monte Carlo

`MonteCarloSimulations.m` saves two `.mat` files to the current folder:

- `nominalResults.mat`
  - size: `(Noutputs × Nmodels) × Nmolecules`
- `MCResults.mat`
  - size: `Nsamples × (Noutputs × Nmodels) × Nmolecules`

In the current code:
- `Noutputs = 4`  (Vss, Cmaxu, AUCu, T0.1µM)
- `Nmodels  = 4`

So the “feature” dimension equals `4 × 4 = 16`.

Outcome blocks (each block has 4 columns, one per model, in the model order above):

- Columns 1–4: `Vss`
- Columns 5–8: `Cmaxu` (tissue-specific)
- Columns 9–12: `AUCu`  (tissue-specific; stored as µM·min)
- Columns 13–16: `T0.1µM` (tissue-specific; stored as hours)

---

## Tissue index (`selectedTissue`)

`selectedTissue` selects the tissue for `Cmaxu`, `AUCu`, and `T0.1µM`.

The tissue list is created inside `PrepareParameterStructure`. To print it, run a small “dummy” parameter struct:

```matlab
P = struct;
P.species  = 'hum';
P.KpuFunc  = @AkpaFarahatKpuModel;

P.RMM    = 300;
P.pKa_a  = 7;
P.pKa_b  = 7;
P.fu     = 0.1;
P.CLint  = 1e-3;
P.logD74 = 2;

[~, P2] = PrepareParameterStructure(P);
disp(P2.tissueNames);
```

Pick the index that matches your tissue of interest.

---

## Correlated Sobol GSA

`SobolForCorrelatedParamsPBPK` expects `ReferenceData` with **N × 6** columns:

(1) molar mass, (2) donor pKa, (3) acceptor pKa, (4) fu, (5) `log10(human inferred CLh_int [L/s])`, (6) logD74.


Example:

```matlab
ReferenceData = referenceProperties;
SobolForCorrelatedParamsPBPK(ReferenceData);
```

---

## Repository contents

### Code files

#### `AkpaFarahatKpuModel.m`
Function: `AkpaFarahatKpuModel(Params)`  
Reads: `tissueCompositionParamsHuman.csv` or `tissueCompositionParamsRat.csv` (species in `Params.species`).

#### `MathewKpuModel.m`
Function: `MathewKpuModel(Params)`  
Reads: `tissueCompositionParamsHuman.csv` or `tissueCompositionParamsRat.csv` (species in `Params.species`).

#### `PearceKpuModel.m`
Function: `PearceKpuModel(Params)`  
Reads:
- `tissueCompositionPearceParamsHuman.csv` or `tissueCompositionPearceParamsRat.csv` (species in `Params.species`)
- `PearceRegressionParams.csv` (slope/intercept for regression corrections used in the calibrated Pearce model)

#### `PearceKpuModelNoCalibration.m`
Function: `PearceKpuModelNoCalibration(Params)`  
Purpose: Pearce model run without regression corrections (uncalibrated variant).  
Reads: `tissueCompositionPearceParamsHuman.csv` or `tissueCompositionPearceParamsRat.csv`.

#### `PrepareParameterStructure.m`
Function: `[Y, Params] = PrepareParameterStructure(Params)`  
Reads: `PhysiologicalParameters.csv` (human) or `PhysiologicalParametersRat.csv` (rat).  
Creates: full physiology + compartment parameters and initial-condition struct `Y`.  
Calls: the Kpu model via the function handle in `Params.KpuFunc`.

#### `PBPKodes.m`
Function: `PBPKodes(time, m, P, Y, bodyWeight)`  
ODE right-hand side used by the PBPK solver.

#### `SimulatePBPK.m`
Function: `SimulatePBPK(DrugP, Tmax)`  
Calls: `PrepareParameterStructure(DrugP)` then solves the ODE system via `PBPKodes`.

#### `MonteCarloSimulations.m`
Function: `MonteCarloSimulations(referenceProperties, selectedTissue, Nsamples, errorScale)`  
Runs Monte Carlo over molecular properties, then saves `MCResults.mat` and `nominalResults.mat`.

#### `SobolForCorrelatedParamsPBPK.m`
Function: `SobolForCorrelatedParamsPBPK(ReferenceData)`  
Computes variance-based sensitivity indices for correlated inputs.

### Data files

#### Molecule property inputs
- `Obach862Molecules.csv` (862 rows)
- `Synthetic10000Molecules.csv` (10,000 rows)

#### Physiology + tissue composition inputs
- `PhysiologicalParameters.csv`
- `PhysiologicalParametersRat.csv`
- `tissueCompositionParamsHuman.csv`
- `tissueCompositionParamsRat.csv`
- `tissueCompositionPearceParamsHuman.csv`
- `tissueCompositionPearceParamsRat.csv`
- `PearceRegressionParams.csv`

---

## Troubleshooting

- **MATLAB cannot find a `.csv` file**  
  Set Current Folder to the repo root and run `addpath(genpath(pwd))`.

- **Error about `parfor`**  
  Replace `parfor` with `for` inside `MonteCarloSimulations.m`.

- **Column name mismatch after `readtable`**  
  Use `readtable(...,'VariableNamingRule','preserve')` and access columns via `T.("...")`.

---

## License

AGPL-3.0 (see `LICENSE`).
