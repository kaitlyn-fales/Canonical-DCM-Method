# MCMC-Based Dynamic Causal Modeling for fMRI Connectivity Repository

This repository contains the code, simulation studies, and real data applications accompanying the Canonical Dynamic Causal Model (CDCM) framework for effective connectivity analysis in task-based fMRI.

CDCM combines:
- a simplified observation model for the BOLD response,
- a piecewise analytic solution to the latent neural ODE,
- Bayesian inference using MCMC,
- and theoretical identifiability results under block-designed experiments.

The repository includes:
- simulation studies evaluating computational performance and robustness to observation-model misspecification,
- applications to publicly available task-based fMRI datasets,
- and scripts for reproducing figures, tables, and model estimation results.

---

**Title:** An MCMC-Based Method for Dynamic Causal Modeling of Effective Connectivity in Functional MRI

**Authors:** Kaitlyn R. Fales, Hyebin Song, Nicole A. Lazar

**Abstract:** Effective connectivity analysis in functional magnetic resonance imaging (fMRI) studies directional interactions among brain regions and experimental stimuli. Dynamic causal modeling (DCM) is a widely used method to estimate effective connectivity, based on a state-space representation consisting of a latent neural signal model and an observation model transforming this signal into the observed blood-oxygen–level-dependent (BOLD) response. A standard DCM combines ordinary differential equation (ODE) dynamics for the latent signal with a complex neural-hemodynamic system for the observation model, and typically uses variational Bayes for parameter estimation. While physically well-motivated, this approach can lead to practical challenges such as inexact solutions and underestimated uncertainty. We introduce Canonical DCM (CDCM), a Markov chain Monte Carlo (MCMC)-based method that adopts a simpler observation model and the No-U-Turn Sampler for posterior sampling. The simpler observation model admits a piecewise analytic solution to the neural ODE, increasing computational efficiency and enabling explicit derivation of sufficient conditions for parameter identifiability. Our results indicate that this approach provides reliable uncertainty quantification and consistent estimation of parameters related to experimental inputs for simulated and real data. We use publicly available data from the Wellcome Centre for Human Neuroimaging and the Human Connectome Project (HCP) to benchmark CDCM against standard DCM methods and examine replicability of estimated connectivity patterns in small- and large-scale neuroimaging settings.

---

## Companion `cdcm` R Package

Many analyses in this repository rely on the companion R package [`cdcm`](https://github.com/kaitlyn-fales/cdcm), which provides functions for:
- preparing CDCM model inputs,
- compiling Stan models,
- running MCMC estimation,
- posterior processing,
- and reproducible CDCM workflows.

Please see the package repository for installation instructions and package documentation.

---

## Repository Structure

### `Attention_Motion/`

Real data application using the SPM attention-to-visual-motion dataset. This folder contains scripts for fitting CDCM and comparator DCM approaches to the attention-to-motion task data, along with code for posterior summaries, figures, and tables.

### `Balloon_Simulation/`

Simulation study evaluating CDCM under three types of misspecification, including misspecification of the observation model. Data are generated under a Balloon-model-based DCM framework and analyzed using CDCM to assess robustness when the simplified CDCM observation model is fit to data generated from a more complex hemodynamic model.

### `Computational_Efficiency_Sim/`

Simulation study assessing computational performance of CDCM. This folder contains scripts used to compare CDCM's piecewise analytic solution to the neural ODE and numerical ODE evaluations and to evaluate how computation time changes with model dimension and simulation settings.

### `HCP_Social_Task/`

Real data application using the Human Connectome Project (HCP) social cognition task. This folder contains scripts for applying CDCM to task-based fMRI data from the HCP social task, including model fitting, posterior summaries, and downstream analyses.

---

## Software Requirements

Analyses are primarily conducted using:
- R
- Stan / CmdStanR
- MATLAB
- SPM12

Some workflows were executed on high-performance computing systems using SLURM job scheduling.

Folder-specific README files provide additional details regarding:
- required software,
- execution order,
- data access,
- and reproduction of manuscript results.

---

## Data Availability

The repository does not redistribute external neuroimaging datasets.

Information regarding data access for each real data application is provided in the corresponding folder-specific README file.

---

## Reproducibility

The repository is organized so that each top-level folder corresponds to a distinct simulation study or real data application. Scripts are generally structured to move from:
1. data preparation,
2. model fitting,
3. posterior processing,
4. and figure/table generation.

Intermediate outputs may need to be generated sequentially depending on the workflow.

---

## Citation

If you use this repository or the CDCM framework in your work, please cite the associated manuscript once available.

 
