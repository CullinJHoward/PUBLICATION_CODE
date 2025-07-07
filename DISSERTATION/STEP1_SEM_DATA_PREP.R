
### LIBRARY
library(dplyr)
library(ggplot2)
library(MplusAutomation)

POUNDTOWN <- 0

# set working directory
if (POUNDTOWN == 1) {
  work_dir <- 'D:\\DISSERTATION\\ANALYSIS\\STATISTICAL_MODELS\\DORRY_ONLY\\'
} else {
  work_dir <- 'F:\\DISSERTATION\\ANALYSIS\\STATISTICAL_MODELS\\DORRY_ONLY\\'
}

setwd(work_dir)


################################################################################
################ PREP FOR FACTOR ANALYSIS - DONE 2.23.25

df <- read.csv("SURVEY_DATA\\ONLY_DORRY_DATA.csv")
names(df)

## SAVE .dat file 

prepareMplusData(df,"DISSERTATION_DORRY_ONLY_MM_2.23.25.dat")


################################################################################
################ MERGE IN NEW VARIABLES SAVED P-FACTOR SCORES - 2.25.25


OG_DF <- read.csv("SURVEY_DATA\\ONLY_DORRY_DATA.csv")
NEW_VARS <- read.csv("SURVEY_DATA\\EXTRA_DORRY_VARS.csv")


# JOINING OLD AND NEW VARIABLE PULL 

MERGED_DF <- full_join(OG_DF, NEW_VARS, by = "ID") 

write.csv(MERGED_DF, "DISSERTATION_DORRY_ONLY_ADD_VARS_MM_2.25.25.csv", row.names=FALSE, na="")
prepareMplusData(MERGED_DF,"DISSERTATION_DORRY_ONLY_ADD_VARS_MM_2.25.25.dat")


