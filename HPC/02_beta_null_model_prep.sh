#!/bin/bash

#SBATCH --job-name=beta_null_prep
#SBATCH --array=1-2
#SBATCH --cpus-per-task=1
#SBATCH --mem=40gb
#SBATCH --time=00:10:00
#SBATCH --mail-type=ALL

# Load in modules
module load anaconda3/2023.09-0
module load r/4.4.0
module load gdal/3.8.3
module load geos/3.12.1
module load proj/9.2.1
module load sqlite/3.43.2

# Run the task for each index in the job array
Rscript HPC/02_beta_null_model_prep.R ${SLURM_ARRAY_TASK_ID}
