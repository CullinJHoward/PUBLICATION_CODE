### LIBRARY
library(dplyr)
library(tidyr)
library(MplusAutomation)

################################################################################
############################ SET DIRECTORY & LOAD PATTERN PARAMETERS

work_dir <- "/home/cjh37695/WAVELET_WORKING/TEST_RESP"
setwd(work_dir)

FOB_SUM <- read.csv("FINAL_FREQUENCY_OPTIMIZED_BAND_SELECTION.csv")


## LOAD IN COSINE ESTIMATES 

## PAIRED DYADS 
COH_EST_PAIRED <- read.csv("/home/cjh37695/WAVELET_WORKING/TEST_RESP/DYAD_COHERENCE_COSINE_SUMMARY.csv")
PHS_EST_PAIRED <- read.csv("/home/cjh37695/WAVELET_WORKING/TEST_RESP/DYAD_PHASE_COSINE_SUMMARY.csv")

## RANDOM DYADS 
COH_EST_RAND <- read.csv("/home/cjh37695/WAVELET_WORKING/TEST_RESP/RESAMPLING_VALIDATION/CV_RANDOM_DYAD_COHERENCE_COSINE_SUMMARY.csv")
PHS_EST_RAND <- read.csv("/home/cjh37695/WAVELET_WORKING/TEST_RESP/RESAMPLING_VALIDATION/CV_RANDOM_DYAD_PHASE_COSINE_SUMMARY.csv")

## TRIM DOWN

COH_EST_PAIRED_RED <- COH_EST_PAIRED %>%
  select(-c("COH_C_CANON", "R2", "AIC", "BIC", "SSE", "converged"))
PHS_EST_PAIRED_RED <- PHS_EST_PAIRED %>%
  select(-c("PHS_C_CANON", "R2_circ", "AIC_circ", "BIC_circ", "SSE_circ",  "converged"))
COH_EST_RAND_RED <- COH_EST_RAND %>%
  select(-c("COH_C_CANON", "R2", "AIC", "BIC", "SSE", "converged"))
PHS_EST_RAND_RED <- PHS_EST_RAND %>%
  select(-c("PHS_C_CANON", "R2_circ", "AIC_circ", "BIC_circ", "SSE_circ",  "converged"))

## ADD PAIRING VARIABLE 

COH_EST_PAIRED_RED$PAIR <- "FAM"
PHS_EST_PAIRED_RED$PAIR <- "FAM"
COH_EST_RAND_RED$PAIR <- "RAND"
PHS_EST_RAND_RED$PAIR <- "RAND"

## MERGE TOGETHER 
COHERENCE <- rbind(COH_EST_PAIRED_RED, COH_EST_RAND_RED)
PHASE <- rbind(PHS_EST_PAIRED_RED, PHS_EST_RAND_RED)

## REMOVE DYADS FAILING QC CHECKS
BAD_IDS <- c(1035)  # VISUAL INSPECTION

COHERENCE$FAIL_QC <- ifelse(COHERENCE$ID %in% BAD_IDS, 1, 0)
COHERENCE <- subset(COHERENCE, FAIL_QC == 0)

PHASE$FAIL_QC <- ifelse(PHASE$ID %in% BAD_IDS, 1, 0)
PHASE <- subset(PHASE, FAIL_QC == 0)

COHERENCE <- COHERENCE %>%
  select(-FAIL_QC)
PHASE <- PHASE %>%
  select(-FAIL_QC)

#### ENSURE WE ONLY LOOK AT DYADS IN THE ORIGINAL PAIRING 
VALID_IDS <- setdiff(unique(COH_EST_PAIRED$ID), "1035")

COHERENCE_RED <- COHERENCE[COHERENCE$ID %in% VALID_IDS, ]
PHASE_RED <- PHASE[PHASE$ID %in% VALID_IDS, ]


## COMPUTE BAND AVERAGE/MIDPOINT 

COHERENCE_RED <- COHERENCE_RED %>%
  mutate(MID_BAND = (period_start + period_end) / 2)

PHASE_RED <- PHASE_RED %>%
  mutate(MID_BAND = (period_start + period_end) / 2)

################################################################################
############################ ADD IN ADDED SURVEY MEASURES 

df <- read.csv("SURVEY_MEASURES/JOINT_ALL_VARIABLES.csv")

## MERGE IN JUST THE PARTICIPANTS NEEDED 

COH_FULL_DF <- left_join(COHERENCE_RED, df, by = "ID")
PHS_FULL_DF <- left_join(PHASE_RED, df, by = "ID")

names(COH_FULL_DF)


################################################################################
############################ PREP VARIABLES FOR MULTI-LEVEL MODEL 

## MUST BE IN LONG FORMAT 

# DROP OUT RSB VARIABLES 

COH_FULL_DF_PAIRED <- subset(COH_FULL_DF, BAND != "RSB")
PHS_FULL_DF_PAIRED <- subset(PHS_FULL_DF, BAND != "RSB")

## GRAND MEAN CENTER FIXED EFFECTS 

#HELPER FUCTION
gmc <- function(df, vars) {
  for (v in vars) {
    df[[paste0("GMC_", v)]] <- df[[v]] - mean(df[[v]], na.rm = TRUE)
  }
  df
}

# APPLIED 
COH_FULL_DF_PAIRED <- gmc(
  df   = COH_FULL_DF_PAIRED,
  vars = c("C_AGE", "INCOME", 
           "ACCp1", "PCONp1", "FCONp1", 
           "ACCc1", "PCONc1", "FCONc1",
           "SOM_P_1", "ANX_P_1", "WTH_P_1", "SOC_P_1", "THO_P_1", 
           "ATT_P_1", "RBK_P_1", "AGG_P_1", "INT_P_1", "EXT_P_1", 
           "TOT_P_1", "FOOD_INS", "SUBSEScm", "SUBSESus",
           "COH_MU",  "COH_A", "COH_B", "COH_C"))

PHS_FULL_DF_PAIRED <- gmc(
  df   = PHS_FULL_DF_PAIRED,
  vars = c("C_AGE", "INCOME", 
           "ACCp1", "PCONp1", "FCONp1", 
           "ACCc1", "PCONc1", "FCONc1",
           "SOM_P_1", "ANX_P_1", "WTH_P_1", "SOC_P_1", "THO_P_1", 
           "ATT_P_1", "RBK_P_1", "AGG_P_1", "INT_P_1", "EXT_P_1", 
           "TOT_P_1", "FOOD_INS", "SUBSEScm", "SUBSESus",
           "PHS_MU", "PHS_A", "PHS_B", "PHS_C"))

COH_FULL_DF_PAIRED_RED <- COH_FULL_DF_PAIRED %>%
  select(-c("NBC3PW1", "NBC5PW1", "NBC9PW1", "CTS30PW1", "NER1PW1", "NER2PW1",     
            "NER3PW1", "NER4PW1", "NER5PW1", "NER6PW1",   "NER9PW1",     
            "NER10PW1", "NER11PW1", "NER13PW1",   "NER14PW1" ))

PHS_FULL_DF_PAIRED_RED <- PHS_FULL_DF_PAIRED %>%
  select(-c("NBC3PW1", "NBC5PW1", "NBC9PW1", "CTS30PW1", "NER1PW1", "NER2PW1",     
            "NER3PW1", "NER4PW1", "NER5PW1", "NER6PW1",   "NER9PW1",     
            "NER10PW1", "NER11PW1", "NER13PW1",   "NER14PW1" ))

## PERSON(DYAD)-MEAN CENTER RANDOM EFFECTS 

pmc <- function(df, vars, id) {
  for (v in vars) {
    df[[paste0("PMC_", v)]] <-
      df[[v]] - ave(df[[v]], df[[id]], FUN = function(x) mean(x, na.rm = TRUE))
  }
  df
}

#APPLIED 

COH_FULL_DF_PAIRED_RED <- pmc(
  df   = COH_FULL_DF_PAIRED_RED,
  vars = c("COH_MU",  "COH_A", "COH_B", "COH_C"),
  id   = "ID"
)

PHS_FULL_DF_PAIRED_RED <- pmc(
  df   = PHS_FULL_DF_PAIRED_RED,
  vars = c("PHS_MU", "PHS_A", "PHS_B", "PHS_C"),
  id   = "ID"
)



## SAVE DATA 

write.csv(COH_FULL_DF_PAIRED_RED, "COHERENCE_FULL_DF_LONG_RAW.csv", row.names = F)
write.csv(PHS_FULL_DF_PAIRED_RED, "PHASE_FULL_DF_LONG_RAW.csv", row.names = F)


################################################################################
############################ PREP VARIABLES FOR SEM PATH MODEL - WIDE FORMAT 

## MUST BE IN WIDE FORMAT

COH_PAIR_FOB <- subset(COH_FULL_DF, BAND == "FOB" & PAIR == "FAM")
COH_PAIR_RSB <- subset(COH_FULL_DF, BAND == "RSB" & PAIR == "FAM")
COH_PAIR_HF <- subset(COH_FULL_DF, BAND == "HF" & PAIR == "FAM")
COH_PAIR_LF <- subset(COH_FULL_DF, BAND == "LF" & PAIR == "FAM")
COH_RAND_FOB <- subset(COH_FULL_DF, BAND == "FOB" & PAIR == "RAND")
COH_RAND_HF <- subset(COH_FULL_DF, BAND == "HF" & PAIR == "RAND")
COH_RAND_LF <- subset(COH_FULL_DF, BAND == "LF" & PAIR == "RAND")

PHS_PAIR_FOB <- subset(PHS_FULL_DF, BAND == "FOB" & PAIR == "FAM")
PHS_PAIR_RSB <- subset(PHS_FULL_DF, BAND == "RSB" & PAIR == "FAM")
PHS_PAIR_HF <- subset(PHS_FULL_DF, BAND == "HF" & PAIR == "FAM")
PHS_PAIR_LF <- subset(PHS_FULL_DF, BAND == "LF" & PAIR == "FAM")
PHS_RAND_FOB <- subset(PHS_FULL_DF, BAND == "FOB" & PAIR == "RAND")
PHS_RAND_HF <- subset(PHS_FULL_DF, BAND == "HF" & PAIR == "RAND")
PHS_RAND_LF <- subset(PHS_FULL_DF, BAND == "LF" & PAIR == "RAND")

## REDUCE VARS 
COH_PAIR_FOB <- COH_PAIR_FOB %>%
  select(c("ID", "COH_MU", "COH_A", "COH_B",  "COH_C", "MID_BAND"))
COH_PAIR_RSB <- COH_PAIR_RSB %>%
  select(c("ID", "COH_MU", "COH_A", "COH_B",  "COH_C", "MID_BAND"))
COH_PAIR_HF <- COH_PAIR_HF %>%
  select(c("ID", "COH_MU", "COH_A", "COH_B",  "COH_C"))
COH_PAIR_LF <- COH_PAIR_LF %>%
  select(c("ID", "COH_MU", "COH_A", "COH_B",  "COH_C"))
COH_RAND_FOB <- COH_RAND_FOB %>%
  select(c("ID", "COH_MU", "COH_A", "COH_B",  "COH_C", "MID_BAND"))
COH_RAND_HF <- COH_RAND_HF %>%
  select(c("ID", "COH_MU", "COH_A", "COH_B",  "COH_C"))
COH_RAND_LF <- COH_RAND_LF %>%
  select(c("ID", "COH_MU", "COH_A", "COH_B",  "COH_C"))

PHS_PAIR_FOB <- PHS_PAIR_FOB %>%
  select(c("ID", "PHS_MU", "PHS_A", "PHS_B",  "PHS_C"))
PHS_PAIR_RSB <- PHS_PAIR_RSB %>%
  select(c("ID", "PHS_MU", "PHS_A", "PHS_B",  "PHS_C"))
PHS_PAIR_HF <- PHS_PAIR_HF %>%
  select(c("ID", "PHS_MU", "PHS_A", "PHS_B",  "PHS_C"))
PHS_PAIR_LF <- PHS_PAIR_LF %>%
  select(c("ID", "PHS_MU", "PHS_A", "PHS_B",  "PHS_C"))
PHS_RAND_FOB <- PHS_RAND_FOB %>%
  select(c("ID", "PHS_MU", "PHS_A", "PHS_B",  "PHS_C"))
PHS_RAND_HF <- PHS_RAND_HF %>%
  select(c("ID", "PHS_MU", "PHS_A", "PHS_B",  "PHS_C"))
PHS_RAND_LF <- PHS_RAND_LF %>%
  select(c("ID", "PHS_MU", "PHS_A", "PHS_B",  "PHS_C"))

## ADD SUFFIX  
names(COH_PAIR_FOB)[names(COH_PAIR_FOB) != "ID"] <- paste0(names(COH_PAIR_FOB)[names(COH_PAIR_FOB) != "ID"], "_FM_FOB")
names(COH_PAIR_RSB)[names(COH_PAIR_RSB) != "ID"] <- paste0(names(COH_PAIR_RSB)[names(COH_PAIR_RSB) != "ID"], "_FM_RSB")
names(COH_PAIR_HF)[names(COH_PAIR_HF) != "ID"] <- paste0(names(COH_PAIR_HF)[names(COH_PAIR_HF) != "ID"], "_FM_HF")
names(COH_PAIR_LF)[names(COH_PAIR_LF) != "ID"] <- paste0(names(COH_PAIR_LF)[names(COH_PAIR_LF) != "ID"], "_FM_LF")
names(COH_RAND_FOB)[names(COH_RAND_FOB) != "ID"] <- paste0(names(COH_RAND_FOB)[names(COH_RAND_FOB) != "ID"], "_RD_FOB")
names(COH_RAND_HF)[names(COH_RAND_HF) != "ID"] <- paste0(names(COH_RAND_HF)[names(COH_RAND_HF) != "ID"], "_RD_HF")
names(COH_RAND_LF)[names(COH_RAND_LF) != "ID"] <- paste0(names(COH_RAND_LF)[names(COH_RAND_LF) != "ID"], "_RD_LF")

names(PHS_PAIR_FOB)[names(PHS_PAIR_FOB) != "ID"] <- paste0(names(PHS_PAIR_FOB)[names(PHS_PAIR_FOB) != "ID"], "_FM_FOB")
names(PHS_PAIR_RSB)[names(PHS_PAIR_RSB) != "ID"] <- paste0(names(PHS_PAIR_RSB)[names(PHS_PAIR_RSB) != "ID"], "_FM_RSB")
names(PHS_PAIR_HF)[names(PHS_PAIR_HF) != "ID"] <- paste0(names(PHS_PAIR_HF)[names(PHS_PAIR_HF) != "ID"], "_FM_HF")
names(PHS_PAIR_LF)[names(PHS_PAIR_LF) != "ID"] <- paste0(names(PHS_PAIR_LF)[names(PHS_PAIR_LF) != "ID"], "_FM_LF")
names(PHS_RAND_FOB)[names(PHS_RAND_FOB) != "ID"] <- paste0(names(PHS_RAND_FOB)[names(PHS_RAND_FOB) != "ID"], "_RD_FOB")
names(PHS_RAND_HF)[names(PHS_RAND_HF) != "ID"] <- paste0(names(PHS_RAND_HF)[names(PHS_RAND_HF) != "ID"], "_RD_HF")
names(PHS_RAND_LF)[names(PHS_RAND_LF) != "ID"] <- paste0(names(PHS_RAND_LF)[names(PHS_RAND_LF) != "ID"], "_RD_LF")

## MERGE TOGETHER COSINE DATA 

ALL_COS <- Reduce(
  function(x, y) merge(x, y, by = "ID", all = TRUE),
  list(
    COH_PAIR_FOB,
    COH_PAIR_RSB,
    COH_PAIR_HF,
    COH_PAIR_LF,
    COH_RAND_FOB,
    COH_RAND_HF,
    COH_RAND_LF,
    PHS_PAIR_FOB,
    PHS_PAIR_RSB,
    PHS_PAIR_HF,
    PHS_PAIR_LF,
    PHS_RAND_FOB,
    PHS_RAND_HF,
    PHS_RAND_LF
  )
)

## ADD IN SURVEY MEASURES 

ALL_COS_VARS <- left_join(ALL_COS, df, by = "ID")

## CENTER THE VARIABLES 

gmc <- function(df, vars) {
  for (v in vars) {
    df[[paste0("GMC_", v)]] <- df[[v]] - mean(df[[v]], na.rm = TRUE)
  }
  df
}


# APPLIED 
CENT_ALL_COS_VARS <- gmc(
  df   = ALL_COS_VARS,
  vars = c("C_AGE", "INCOME", 
           "ACCp1", "PCONp1", "FCONp1", 
           "ACCc1", "PCONc1", "FCONc1",
           "SOM_P_1", "ANX_P_1", "WTH_P_1", "SOC_P_1", "THO_P_1", 
           "ATT_P_1", "RBK_P_1", "AGG_P_1", "INT_P_1", "EXT_P_1", 
           "TOT_P_1", "FOOD_INS", "SUBSEScm", "SUBSESus",
           "COH_MU_FM_FOB",   "COH_A_FM_FOB",    "COH_B_FM_FOB",    "COH_C_FM_FOB",    "MID_BAND_FM_FOB",
           "COH_MU_FM_RSB",   "COH_A_FM_RSB",    "COH_B_FM_RSB",    "COH_C_FM_RSB",    "MID_BAND_FM_RSB", "COH_MU_FM_HF",   
           "COH_A_FM_HF",     "COH_B_FM_HF",     "COH_C_FM_HF",     "COH_MU_FM_LF",    "COH_A_FM_LF",     "COH_B_FM_LF",    
           "COH_C_FM_LF",     "COH_MU_RD_FOB",   "COH_A_RD_FOB",    "COH_B_RD_FOB",    "COH_C_RD_FOB",    "MID_BAND_RD_FOB",
           "COH_MU_RD_HF",    "COH_A_RD_HF",     "COH_B_RD_HF",     "COH_C_RD_HF",     "COH_MU_RD_LF",    "COH_A_RD_LF",    
           "COH_B_RD_LF",     "COH_C_RD_LF",     "PHS_MU_FM_FOB",   "PHS_A_FM_FOB",    "PHS_B_FM_FOB",    "PHS_C_FM_FOB",   
           "PHS_MU_FM_RSB",   "PHS_A_FM_RSB",    "PHS_B_FM_RSB",    "PHS_C_FM_RSB",    "PHS_MU_FM_HF",    "PHS_A_FM_HF",    
           "PHS_B_FM_HF",     "PHS_C_FM_HF",     "PHS_MU_FM_LF",    "PHS_A_FM_LF",     "PHS_B_FM_LF",     "PHS_C_FM_LF",    
           "PHS_MU_RD_FOB",   "PHS_A_RD_FOB",    "PHS_B_RD_FOB",    "PHS_C_RD_FOB",    "PHS_MU_RD_HF",    "PHS_A_RD_HF",    
           "PHS_B_RD_HF",     "PHS_C_RD_HF",     "PHS_MU_RD_LF",    "PHS_A_RD_LF",     "PHS_B_RD_LF",     "PHS_C_RD_LF"  ))

## REMOVE SOME SHITS 

CENT_ALL_COS_VARS_RED <- CENT_ALL_COS_VARS %>%
  select(-c("NBC3PW1", "NBC5PW1", "NBC9PW1", "CTS30PW1", "NER1PW1", "NER2PW1",     
            "NER3PW1", "NER4PW1", "NER5PW1", "NER6PW1",   "NER9PW1",     
            "NER10PW1", "NER11PW1", "NER13PW1",   "NER14PW1" ))

## SAVE WIDE DF FOR CSV and DAT 

write.csv(CENT_ALL_COS_VARS_RED, "COHERENCE_FULL_DF_WIDE_RAW.csv", row.names = F)
prepareMplusData(CENT_ALL_COS_VARS_RED,"COHERENCE_FULL_DF_WIDE_RAW.dat", inpfile =T)
