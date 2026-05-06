# Balloon_Simulation

This folder contains simulation studies evaluating CDCM under model misspecification.

The simulations investigate the performance of CDCM when data are generated from Dynamic Causal Models using the nonlinear Balloon hemodynamic observation model and then analyzed using the simplified CDCM observation model. Comparisons with standard SPM DCM estimation procedures are also included.

The simulation study evaluates:
- parameter recovery,
- posterior uncertainty,
- predictive performance,
- and robustness of CDCM under misspecification.

All simulation settings are evaluated across 50 independent replicates.

---

## Folder Structure

### `Data/`

Directory used to store generated simulation datasets and intermediate files.

The repository keeps this directory mostly empty to avoid storing large simulation outputs in version control.

### `Output/`

Directory used to store:
- posterior summaries,
- predictive results,
- figures,
- and simulation comparison outputs.

---

## Simulation Structure

The simulation study considers multiple model configurations involving:
- canonical (CDCM) versus Balloon (SPM) observation models,
- models with or without modulatory effects in the diagonal elements of $\mathbf{B}$ to reflect the differences in diagonal parameterization between CDCM and SPM,
- and differences in the type of differential equation used (ODE vs DDE)

The naming convention of scripts reflects the corresponding simulation setting.

For example:
- `diagA` indicates simulations involving diagonal $\mathbf{A}$ effects,
- `diagAB` indicates simulations involving diagonal effects in both $\mathbf{A}$ and $\mathbf{B}$,
- `zero` and `nonzero` indicate, using SPM alone, whether the model uses an ODE (no signal delays between regions) or a DDE (TR/2 signal delays).

---

## Main Simulation Scripts

### Data Generation

#### `simulate_balloon_mdl_diagA.m`
#### `simulate_balloon_mdl_diagAB.m`

Generate datasets under Balloon-model-based DCM settings in MATLAB/SPM.

#### `simulate_balloon_mdl_diagA.R`
#### `simulate_balloon_mdl_diagAB.R`

Corresponding R scripts to add Gaussian noise at a specified SNR for the above MATLAB scripts. 

#### `simulate_canonical_mdl_diagA.R`
#### `simulate_canonical_mdl_diagAB.R`

Generate datasets under CDCM.

---

## CDCM Estimation Scripts

### `balloon_mdl_diagA_zero.R`
### `balloon_mdl_diagA_nonzero.R`
### `balloon_mdl_diagAB_zero.R`
### `balloon_mdl_diagAB_nonzero.R`

Fit CDCM models to datasets generated under Balloon-model simulation settings (model misspecification).

### `canonical_mdl_diagA.R`
### `canonical_mdl_diagAB.R`

Fit CDCM models to datasets generated under CDCM (correctly specified model).

---

## SPM Estimation Scripts

### `balloon_mdl_diagA_zero.m`
### `balloon_mdl_diagA_nonzero.m`
### `balloon_mdl_diagAB_zero.m`
### `balloon_mdl_diagAB_nonzero.m`

MATLAB/SPM estimation scripts corresponding to the Balloon-model simulation settings (correctly specified model).

### `canonical_mdl_diagA.m`
### `canonical_mdl_diagAB.m`

MATLAB/SPM estimation scripts corresponding to CDCM generated data (model misspecification).

---

## SLURM Scripts

### `*.sh`

SLURM submission scripts used to execute the simulation analyses on a high-performance computing cluster across simulation replicates.

---

## Results and Comparison Scripts

### `results_comparison_MCMC.R`

Processes and summarizes CDCM posterior estimation results across simulation replicates.

### `results_comparison_SPM.R`

Processes and summarizes SPM estimation results across simulation replicates.

### `generate_pred_signal_MCMC.R`

Generates posterior predictive trajectories from CDCM estimation output to use in `mse_comparison.R`.

### `mse_comparison.R`

Computes predictive and estimation error summaries across competing methods.

### `canonical_balloon_scales.R`

Compares the amplitude differences between the simulated CDCM and SPM signals under the same neural parameters. 

---

## Software Requirements

The analyses primarily rely on:
- R,
- CmdStanR / Stan,
- MATLAB,
- and SPM12.

Several scripts additionally require the companion `cdcm` R package:

https://github.com/kaitlyn-fales/cdcm

---

## General Workflow

A typical workflow is:

1. Generate simulation datasets.
2. Run CDCM estimation scripts.
3. Run comparator SPM estimation scripts.
4. Aggregate simulation results across the 50 replicates.
5. Generate summary statistics and visualization outputs.

Depending on the workflow, scripts may need to be executed sequentially.

---

## Notes

Several analyses were originally executed on high-performance computing systems using SLURM job scheduling.

Paths within scripts may need to be modified depending on the local computing environment.
