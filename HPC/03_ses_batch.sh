#!/bin/bash

#SBATCH --job-name=ses_batch
#SBATCH --array=1-28
#SBATCH --cpus-per-task=1
#SBATCH --mem=22gb
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
Rscript HPC/03_hpc_batch_ses.R ${SLURM_ARRAY_TASK_ID}
