# HCP_Social_Task

This folder contains the Human Connectome Project (HCP) social cognition task application for the CDCM framework.

The analyses apply CDCM to task-based fMRI data from the HCP social cognition task to evaluate effective connectivity during social processing. The workflow includes:
- preprocessing and ROI preparation,
- single-subject CDCM estimation,
- posterior diagnostics and quality checks,
- and downstream group-level analyses.

We also compare CDCM with standard SPM Dynamic Causal Modeling procedures, additionally included in the `SPM/` directory.

---

## Dataset

The analyses are based on the Human Connectome Project Young Adult (HCP-YA) dataset. We use the 100 Unrelated Subjects Dataset from the HCP-YA 2025 release, available from BALSA.

Information regarding the HCP-YA project and data access is available here:

- https://www.humanconnectome.org/study/hcp-young-adult

Access to HCP data requires registration and agreement to the HCP data use terms.

This repository does not redistribute the original neuroimaging data.

---

## Folder Structure

### `Analysis/`

Contains downstream CDCM analyses, including:
- posterior diagnostics,
- subject-level quality checks,
- posterior summary processing,
- visualization scripts,
- and group-level analyses.

A folder-specific README provides additional details regarding the analysis workflow.

### `SPM/`

Contains comparator analyses using SPM Dynamic Causal Modeling procedures, including:
- single-subject estimation,
- posterior summaries,
- diagnostics,
- and group-level analyses.

The `SPM/` directory additionally contains its own `Data/` and `Output/` subdirectories for storing intermediate `.mat` files and generated results.

A folder-specific README provides additional details regarding the SPM workflow.

### `Data/`

Directory used to store processed CDCM inputs and intermediate data objects.

The repository keeps this directory mostly empty to avoid storing large derived datasets in version control.

### `Output/`

Directory used to store:
- posterior summaries,
- subject-level outputs,
- diagnostics,
- and generated analysis results for CDCM.

---

## Main Scripts

### `preprocess_data.R`

Preprocesses and organizes HCP social task data for CDCM estimation.

### `VOI_mask_generate.R`

Generates ROI mask files used in the CDCM analyses.

### `single-sub_phaseLR_maskL.R`
### `single-sub_phaseLR_maskR.R`
### `single-sub_phaseRL_maskL.R`
### `single-sub_phaseRL_maskR.R`

Single-subject CDCM estimation scripts corresponding to:
- left/right hemisphere ROI masks,
- and LR/RL HCP acquisition phase encodings.

### `single-sub_phaseLR_maskL.sh`
### `single-sub_phaseLR_maskR.sh`
### `single-sub_phaseRL_maskL.sh`
### `single-sub_phaseRL_maskR.sh`

SLURM submission scripts used to execute the single-subject CDCM analyses on high-performance computing systems.

### `group_results.R`

Processes posterior outputs and generates group-level CDCM summaries and visualizations.

---

## Additional Files

### `mask_L.nii`
### `mask_R.nii`

ROI mask files used for extracting regional time series and constructing CDCM inputs.

### `file_list_*.txt`

Text files containing subject-level file paths and input lists used during preprocessing and model estimation workflows.

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

A typical CDCM workflow is:

1. Obtain and organize the HCP social cognition task data.
2. Generate ROI masks and preprocess subject-level inputs.
3. Run single-subject CDCM estimation.
4. Perform posterior diagnostics and quality checks.
5. Aggregate posterior summaries across subjects.
6. Generate group-level analyses and visualizations.

Comparator SPM analyses are organized separately within the `SPM/` directory.

---

## Notes

Several analyses were originally executed on high-performance computing systems using SLURM job scheduling.

Paths within scripts may need to be modified depending on the local computing environment.
