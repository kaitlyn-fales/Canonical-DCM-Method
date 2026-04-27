#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1 
#SBATCH --mem-per-cpu=4gb
#SBATCH --time=24:00:00
#SBATCH --account=statsresearch_sc_default
#SBATCH --output=Output/canonical_mdl_diagA_rep%03a.out
#SBATCH --array=1-50

# Get started
echo " "
echo "Job started on `hostname` at `date`"
echo " "

# Environment setup
module purge
module use /gpfs/group/RISE/sw7/modules
module load r/4.2.1

export R_LIBS="~/R_packages"

Rscript canonical_mdl_diagA.R "index=$SLURM_ARRAY_TASK_ID"

# Finish up
echo " "
echo "Job Ended at `date`"
echo " "