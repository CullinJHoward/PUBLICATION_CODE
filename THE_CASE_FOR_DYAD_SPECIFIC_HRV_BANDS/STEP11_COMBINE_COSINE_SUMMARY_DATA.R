##LIBRARY STUFF 
library(dplyr)


## SET WORKING DIRECTORY
work_dir <- "/home/cjh37695/WAVELET_WORKING/TEST_RESP//"
setwd(work_dir)


################################################################################
################## LOAD AND BIND ALL COHERENCE SUMMARY DATA
COHERENCE_FILES <- "COHERENCE_COSINE_SUMMARY_PARTS/"
COH_sum_files  <- list.files(COHERENCE_FILES, pattern = "^DYAD_COHERENCE_COSINE_SUMMARY_part_\\d+\\.csv$", full.names =  TRUE)
COH_full_files <- list.files(COHERENCE_FILES, pattern = "^DYAD_COHERENCE_COSINE_FULL_OUTPUT_part_\\d+\\.csv$", full.names =  TRUE)

COH_DYAD_SUM  <- bind_rows(lapply(COH_sum_files, read.csv, stringsAsFactors = FALSE))
COH_DYAD_FULL <- bind_rows(lapply(COH_full_files, read.csv, stringsAsFactors = FALSE))

## SAVE THE FULL SUMMARY TO MAIN WORKING DIRECTORY 
write.csv(COH_DYAD_SUM,  "DYAD_COHERENCE_COSINE_SUMMARY.csv", row.names = FALSE)
write.csv(COH_DYAD_FULL, "DYAD_COHERNECE_COSINE_FULL_OUTPUT.csv", row.names = FALSE)



################################################################################
################## LOAD AND BIND ALL PHASE SUMMARY DATA
PHASE_FILES <- "PHASE_COSINE_SUMMARY_PARTS/"
PSH_sum_files  <- list.files(PHASE_FILES, pattern = "^DYAD_PHASE_CIRCULAR_COSINE_SUMMARY_part_\\d+\\.csv$", full.names =  TRUE)
PSH_full_files <- list.files(PHASE_FILES, pattern = "^DYAD_PHASE_CIRCULAR_COSINE_FULL_OUTPUT_part_\\d+\\.csv$", full.names =  TRUE)

PSH_DYAD_SUM  <- bind_rows(lapply(PSH_sum_files, read.csv, stringsAsFactors = FALSE))
PSH_DYAD_FULL <- bind_rows(lapply(PSH_full_files, read.csv, stringsAsFactors = FALSE))

## SAVE THE FULL SUMMARY TO MAIN WORKING DIRECTORY 
write.csv(PSH_DYAD_SUM,  "DYAD_PHASE_COSINE_SUMMARY.csv", row.names = FALSE)
write.csv(PSH_DYAD_FULL, "DYAD_PHASE_COSINE_FULL_OUTPUT.csv", row.names = FALSE)



