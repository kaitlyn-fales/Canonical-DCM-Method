# Attention_Motion

This folder contains the attention-to-visual-motion real data application for the Canonical Dynamic Causal Model (CDCM) framework.

The analyses use the publicly available SPM attention-to-visual-motion dataset and compare CDCM with standard Dynamic Causal Modeling approaches implemented in:
- SPM,
- and Variational Bayesian Analysis (VBA) toolbox.

The application is used to evaluate:
- posterior estimation behavior,
- uncertainty quantification,
- agreement between CDCM and existing DCM implementations,
- and computational characteristics of CDCM in a realistic task-based fMRI setting.

---

## Dataset

The analyses are based on the publicly available SPM attention-to-visual-motion dataset distributed with SPM12. The tutorial users should follow to recreate the SPM results from scratch is in Chapter 36 of the SPM12 manual. The key outputs are included in this repository.

SPM12 and the associated tutorial datasets are available here:

- https://www.fil.ion.ucl.ac.uk/spm/software/spm12/
- https://www.fil.ion.ucl.ac.uk/spm/data/attention/

The dataset contains task-based fMRI recordings from a visual attention experiment involving regions including:
- V1,
- V5,
- and SPC.

This repository does not redistribute the original neuroimaging dataset. Users should obtain the data directly through SPM12.

---

## Main Scripts

### `prepare_dat.R`

Prepares the attention-to-motion data for CDCM estimation and exports the required inputs for MCMC sampling.

### `attention_motion.R`

Main CDCM estimation script for the attention-to-motion application.

### `attention_motion.sh`

SLURM submission script used to execute CDCM estimation on a high-performance computing cluster.

### `results_visualization.R`

Generates posterior summaries, figures, and visualization outputs from CDCM estimation results.

### `est_vs_obs.R`

Compares observed and estimated BOLD trajectories and computes uncertainty summaries using block bootstrap procedures.

### `est_vs_obs_z.R`

Additional observed-versus-estimated trajectory analyses on the latent neural signal scale.

### `spline_regression.R`

Estimates signal-to-noise characteristics from the attention-to-motion dataset using spline-based signal estimation. The result of this script is the SNR used in all simulation studies. 

---

## `SPM_VBA_Comparison/`

This directory contains scripts and outputs used for comparison with:
- SPM DCM,
- and VBA implementations.

Key contents include:
- MATLAB scripts for VBA estimation,
- processed SPM outputs,
- posterior summaries,
- and comparison tables.

### Important Files

#### `VBA_motion_model.m`

MATLAB script used to fit the VBA model.

#### `SPM_results.R`

Processes and summarizes SPM DCM estimation results.

#### `VBA_results.R`

Processes and summarizes VBA estimation results.

#### `DCM_mod_fwd.mat`

Modified forward DCM structure used for model fitting.

#### `VOI_*.mat`

Region-specific VOI files used in the DCM analyses.

---

## `Output/`

Directory used to store generated outputs and analysis results.

The repository keeps this directory mostly empty so users can generate outputs locally without storing large intermediate files in version control.

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

1. Obtain the SPM attention-to-visual-motion dataset.
2. Run `prepare_dat.R`.
3. Execute `attention_motion.R` (optionally through `attention_motion.sh`).
4. Run post-processing and visualization scripts.
5. Generate comparison summaries using the `SPM_VBA_Comparison/` directory.

---

## Notes

Some scripts contain paths that may need to be modified depending on the local computing environment.

Several analyses were originally executed on high-performance computing systems using SLURM job scheduling.
