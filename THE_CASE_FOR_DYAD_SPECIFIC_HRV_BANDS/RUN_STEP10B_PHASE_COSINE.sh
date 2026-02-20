#!/bin/bash

#SBATCH --job-name=Phase_Cosine_Analysis

#SBATCH --output=PHS_COSINE_logs/PHS_COSINE_%A_%a.out

#SBATCH --error=PHS_COSINE_logs/PHS_COSINE_%A_%a.err

#SBATCH --array=1-29%15

#SBATCH --time=3:00:00              # adjust based on expected runtime

#SBATCH --partition=batch            # Sapelo partition
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

#SBATCH --mem=16G                     # adjust as needed

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

#   Run the R script

# ---------------------------
export N_TASKS=29

Rscript TR_STEP10B_DYAD_LEVEL_COSINE_PATTERNS_PHASE_SLURM_FRIENDLY.R

