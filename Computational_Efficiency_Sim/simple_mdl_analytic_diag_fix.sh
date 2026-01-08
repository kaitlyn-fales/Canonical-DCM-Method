#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1 
#SBATCH --mem-per-cpu=4gb
#SBATCH --time=5:00:00
#SBATCH --account=open
#SBATCH --partition=open
#SBATCH --output=Output/simple_mdl_analytic_diag_fix_%a.out
#SBATCH --array=1-5

# Get started
echo " "
echo "Job started on `hostname` at `date`"
echo " "

# Environment setup
module purge
module use /gpfs/group/RISE/sw7/modules
module load r/4.2.1

export R_LIBS="~/R_packages"

Rscript simple_mdl_analytic_diag_fix.R "index=$SLURM_ARRAY_TASK_ID"

# Finish up
echo " "
echo "Job Ended at `date`"
echo " "