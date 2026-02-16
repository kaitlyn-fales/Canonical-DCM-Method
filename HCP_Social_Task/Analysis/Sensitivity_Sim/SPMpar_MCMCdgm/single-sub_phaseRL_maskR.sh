#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1 
#SBATCH --mem-per-cpu=4gb
#SBATCH --time=48:00:00
#SBATCH --account=open
#SBATCH --output=Output/output_phaseRL_maskR_%A_%a.out
#SBATCH --array=1-100

# Get started
echo " "
echo "Job started on `hostname` at `date`"
echo " "

# Environment setup
module purge
module use /gpfs/group/RISE/sw7/modules
module load r/4.2.1

export R_LIBS="~/R_packages"

# Path to file list
FILE_LIST="/storage/work/krf5429/Canonical-DCM-Method/HCP_Social_Task/Analysis/Sensitivity_Sim/SPMpar_MCMCdgm/file_list_phaseRL_maskR.txt"

# Pick the file corresponding to this array task
FILE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $FILE_LIST)

echo "Processing $FILE"

Rscript single-sub_phaseRL_maskR.R "$FILE"

# Finish up
echo " "
echo "Job Ended at `date`"
echo " "






