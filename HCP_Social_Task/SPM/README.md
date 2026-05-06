# SPM

This directory contains comparator analyses for the HCP social cognition task application using standard Dynamic Causal Modeling procedures implemented in SPM.

The workflows in this folder parallel the primary CDCM analyses and include:
- single-subject SPM DCM estimation,
- posterior diagnostics,
- subject-level result summaries,
- group-level analyses,
- and supplementary sensitivity simulations.

The SPM analyses are used throughout the project as a benchmark comparison against CDCM estimation results.

---

## Folder Structure

### `Data/`

Directory used to store:
- processed SPM inputs,
- intermediate `.mat` files,
- and subject-level DCM structures.

The repository keeps this directory mostly empty to avoid storing large derived datasets in version control.

### `Output/`

Directory used to store:
- subject-level outputs,
- estimated DCM structures,
- diagnostics,
- and generated analysis results.

### `Results/`

Directory used to store:
- processed posterior summaries,
- compiled group-level outputs,
- and visualization results.

### `Sensitivity_Sim/`

Contains supplementary simulation studies used to evaluate sensitivity of group-level inference under different combinations of:
- generating DCM framework,
- and group-level posterior parameter settings.

The `Sensitivity_Sim/` directory in this folder contains the SPM-generated simulation settings:
- `MCMCpar_SPMdgm/`
- `SPMpar_SPMdgm/`

The corresponding CDCM-generated simulation settings are located in the parallel `Analysis/Sensitivity_Sim/` directory.

Each simulation subdirectory follows the same general structure as the primary real-data workflow, including:
- preprocessing,
- single-subject estimation,
- posterior processing,
- and group-level analyses.

---

## Main Scripts

### `single_sub.m`

MATLAB/SPM script used for single-subject Dynamic Causal Model estimation.

### `single_sub.sh`

SLURM submission script used to execute single-subject SPM analyses on high-performance computing systems.

### `single_sub_results.R`

Processes and summarizes single-subject SPM estimation results for comparison with CDCM outputs.

### `diagnostics_compilation.R`

Compiles subject-level diagnostic summaries and convergence checks across estimated SPM DCM models.

### `group_results_SPM.R`

Processes posterior summaries and generates group-level result outputs from the SPM analyses.

### `group_analysis_SPM.R`

Performs downstream group-level analyses using SPM estimation results.

### `group_analysis_SPM.sh`

SLURM submission script used to execute group-level analyses on high-performance computing systems.

---

## Additional Files

### `subjects.txt`

Text file containing the subject list used during single-subject estimation workflows.

---

## Software Requirements

The analyses primarily rely on:
- MATLAB,
- SPM12,
- and R for downstream processing and visualization.

Several workflows additionally rely on outputs generated from the CDCM analyses for comparison purposes.

---

## General Workflow

A typical workflow is:

1. Prepare subject-level SPM inputs.
2. Run single-subject DCM estimation using SPM.
3. Compile diagnostic summaries and subject-level outputs.
4. Generate posterior summaries and group-level results.
5. Perform downstream group-level analyses.
6. Execute supplementary sensitivity simulations where applicable.

The sensitivity simulation workflows mirror the structure of the primary real-data analyses.

---

## Notes

Several analyses were originally executed on high-performance computing systems using SLURM job scheduling.

Paths within scripts may need to be modified depending on the local computing environment.
