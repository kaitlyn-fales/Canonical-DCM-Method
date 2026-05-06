# Analysis

This directory contains downstream CDCM analyses for the HCP social cognition task application.

The analyses in this folder include:
- posterior diagnostics,
- subject-level result summaries,
- group-level modeling,
- visualization and result compilation,
- and additional sensitivity simulations described in the supplementary materials.

The workflows in this directory primarily operate on posterior outputs generated from the single-subject CDCM analyses in the parent `HCP_Social_Task/` directory.

---

## Folder Structure

### `Results/`

Directory used to store:
- processed posterior summaries,
- compiled diagnostics,
- group-level outputs,
- and generated analysis results.

The repository keeps this directory mostly empty to avoid storing large derived outputs in version control.

### `Sensitivity_Sim/`

Contains additional downstream simulation studies used to evaluate sensitivity of group-level inference under different combinations of:
- generating DCM framework,
- and group-level posterior parameter settings.

The sensitivity analyses include four simulation configurations across CDCM and SPM settings:
- CDCM parameters with CDCM-generated data,
- SPM parameters with CDCM-generated data,
- CDCM parameters with SPM-generated data,
- and SPM parameters with SPM-generated data.

The `Sensitivity_Sim/` directory in this folder contains the CDCM-generated simulation settings:
- `MCMCpar_MCMCdgm/`
- `SPMpar_MCMCdgm/`

The corresponding SPM-generated simulation settings are located in the parallel `SPM/Sensitivity_Sim/` directory.

Each simulation subdirectory follows the same general structure as the primary real-data workflow, including:
- preprocessing,
- single-subject estimation,
- posterior processing,
- and group-level analyses.

---

## Main Scripts

### `single_sub_results.R`

Processes and summarizes single-subject CDCM estimation results.

### `diagnostics_compilation.R`

Compiles posterior diagnostics and model quality summaries across subjects.

This includes summaries related to:
- convergence,
- effective sample sizes,
- and posterior sampling diagnostics.

### `group_results_MCMC.R`

Generates group-level posterior summaries and visualization outputs from CDCM estimation results.

### `group_analysis_MCMC.R`

Performs downstream group-level analyses using CDCM posterior outputs.

### `group_analysis_MCMC.sh`

SLURM submission script used to execute group-level analyses on high-performance computing systems.

---

## Software Requirements

The analyses primarily rely on:
- R,
- CmdStanR / Stan,
- and the companion `cdcm` R package:

https://github.com/kaitlyn-fales/cdcm

Several workflows additionally rely on MATLAB/SPM-generated intermediate outputs from the comparator analyses.

---

## General Workflow

A typical workflow is:

1. Process single-subject CDCM outputs.
2. Compile posterior diagnostics across subjects.
3. Generate posterior summaries and visualization outputs.
4. Perform group-level analyses.
5. Execute supplementary sensitivity simulations where applicable.

The sensitivity simulation workflows mirror the structure of the primary real-data analyses.

---

## Notes

Several analyses were originally executed on high-performance computing systems using SLURM job scheduling.

Paths within scripts may need to be modified depending on the local computing environment.
