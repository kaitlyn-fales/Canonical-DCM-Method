# Computational_Efficiency_Sim

This folder contains simulation studies evaluating the computational performance of CDCM under varying model structures and estimation strategies.

The simulations compare:
- analytic evaluation of the latent neural ODE solution,
- and numerical ODE evaluation using the CKRK solver within `cmdstanr`.

The study examines how computational behavior changes across:
- simple versus complex model structures,
- fixed versus varying diagonal self-connectivity settings,
- and different numerical solution approaches.

All simulation settings are evaluated across 50 simulation replicates.

---

## Folder Structure

### `Data/`

Directory used to store generated simulation datasets and intermediate files.

The repository keeps this directory mostly empty to avoid storing large simulation outputs in version control.

### `Output/`

Directory used to store:
- posterior summaries,
- runtime outputs,
- comparison statistics,
- and visualization results.

---

## Simulation Structure

The simulations compare computational performance across:
- analytic versus numerical ODE evaluation,
- simple versus complex connectivity structures,
- and fixed versus varying diagonal self-connectivity configurations.

Script naming conventions reflect the simulation settings.

For example:
- `simple_mdl` and `complex_mdl` indicate the underlying simulation model complexity,
- `analytic` indicates use of the closed-form CDCM neural ODE solution,
- `ckrk` indicates numerical ODE evaluation using Runge-Kutta methods,
- `diag_fix` indicates fixed diagonal self-connectivity,
- and `diag_var` indicates varying diagonal self-connectivity.

---

## Data Generation Scripts

### `simulate_simple_mdl.R`

Generates datasets under the simple connectivity model configuration.

### `simulate_complex_mdl.R`

Generates datasets under the complex connectivity model configuration.

---

## CDCM Estimation Scripts

### Analytic ODE Evaluation

#### `simple_mdl_analytic_diag_fix.R`
#### `simple_mdl_analytic_diag_var.R`

Fit CDCM models using the analytic neural ODE solution under simple model configurations.

#### `complex_mdl_analytic_diag_fix.R`
#### `complex_mdl_analytic_diag_var.R`

Fit CDCM models using the analytic neural ODE solution under complex model configurations.

---

### Numerical ODE Evaluation (CKRK in Stan)

#### `simple_mdl_ckrk_diag_fix.R`
#### `simple_mdl_ckrk_diag_var.R`

Fit CDCM models using numerical ODE evaluation under simple model configurations.

#### `complex_mdl_ckrk_diag_fix.R`
#### `complex_mdl_ckrk_diag_var.R`

Fit CDCM models using numerical ODE evaluation under complex model configurations.

---

## SLURM Scripts

### `*.sh`

SLURM submission scripts used to execute the simulation analyses on high-performance computing systems.

Each estimation script has a corresponding SLURM submission script for parallel execution across simulation replicates.

---

## Results and Comparison Scripts

### `results_comparison_diag_fix.R`

Processes and summarizes simulation results for fixed diagonal self-connectivity settings.

### `results_comparison_diag_var.R`

Processes and summarizes simulation results for varying diagonal self-connectivity settings.

The comparison scripts aggregate:
- runtime summaries,
- and posterior estimation results.

---

## Software Requirements

The analyses primarily rely on:
- R,
- and CmdStanR / Stan.

Several scripts additionally require the companion `cdcm` R package:

https://github.com/kaitlyn-fales/cdcm

---

## General Workflow

A typical workflow is:

1. Generate simulation datasets.
2. Run CDCM estimation using analytic ODE evaluation.
3. Run CDCM estimation using numerical ODE evaluation.
4. Aggregate runtime and estimation results across simulation replicates.
5. Generate summary statistics and comparison outputs.

Depending on the workflow, scripts may need to be executed sequentially.

---

## Notes

Several analyses were originally executed on high-performance computing systems using SLURM job scheduling.

Paths within scripts may need to be modified depending on the local computing environment.
