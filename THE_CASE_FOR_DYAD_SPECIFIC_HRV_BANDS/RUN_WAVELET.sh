#!/bin/bash

#SBATCH --job-name=wavelet_test

#SBATCH --output=logs/wavelet_test_%A_%a.out

#SBATCH --error=logs/wavelet_test_%A_%a.err

#SBATCH --array=1-180%30

#SBATCH --time=04:00:00              # adjust based on expected runtime

#SBATCH --partition=batch            # Sapelo partition
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

#SBATCH --mem=128G                     # adjust as needed

#SBATCH --mail-type=END,FAIL

#SBATCH --mail-user=cjh37695@uga.edu



# ---------------------------

#   Load R module

# ---------------------------

module purge

module load R/4.3.2-gfbf-2023a



# ---------------------------

#   Set personal R library

# ---------------------------

export R_LIBS_USER=/home/cjh37695/R/x86_64-pc-linux-gnu-library/4.3



# ---------------------------

#   Go to working directory

# ---------------------------

cd /home/cjh37695/WAVELET_WORKING/TEST_RESP



# ---------------------------

#   Prepare file list

# ---------------------------

FILELIST=filelist.txt

INPUT_CSV=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $FILELIST)



# ---------------------------

#   Run the R script

# ---------------------------

Rscript TR_STEP5_CROSS_WAVELET_ANALYSIS.R  "$INPUT_CSV"

