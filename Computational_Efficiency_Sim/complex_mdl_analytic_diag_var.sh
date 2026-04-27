#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1 
#SBATCH --mem-per-cpu=4gb
#SBATCH --time=48:00:00
#SBATCH --account=open
#SBATCH --output=Output/complex_mdl_analytic_diag_var_%a.out
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

Rscript complex_mdl_analytic_diag_var.R "index=$SLURM_ARRAY_TASK_ID"

# Finish up
echo " "
echo "Job Ended at `date`"
echo " "