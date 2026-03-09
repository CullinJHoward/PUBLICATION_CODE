################################################################################
############### Load Packages


library("dplyr")
library("psych")
library("MplusAutomation")


################################################################################
############### Set working directory 

POUNDTOWN <- 0

# set working directory
if (POUNDTOWN == 1) {
  work_dir <- 'C:\\Users\\cjh37695\\OneDrive - University of Georgia\\Suicide sleep paper\\ANALYSIS\\ANALYSIS_R&R\\DATA_WRANGLE\\'
} else {
  work_dir <- 'C:\\Users\\0910h\\OneDrive - University of Georgia\\Suicide sleep paper\\ANALYSIS\\ANALYSIS_R&R\\DATA_WRANGLE\\'
}

setwd(work_dir)


################################################################################
############## Load Data Frame 


OG <- read.csv("ABCD_SUD_LINEAR_MODEL_3.10.25.csv")

OG <- OG %>% 
  rename(subID = subid)

SUP <- read.csv("ABCD_SES_DMN_SUI_6.11.25.csv")

## LOAD DATA FRAME 
FULL <- left_join(OG, SUP, by = "subID")


################################################################################
############### COMPUTE LAG  

FULL$LAG_1_3 <- (FULL$CORE_VD_3 - FULL$CORE_VD_1)
FULL$LAG_1_5 <- (FULL$CORE_VD_5 - FULL$CORE_VD_1)
FULL$LAG_3_5 <- (FULL$CORE_VD_3 - FULL$CORE_VD_1)

## AGE 

mean(FULL$Y_AGE_1, na.rm = TRUE)/12 # 9.94 (SD = .63)
sd(FULL$Y_AGE_1, na.rm = TRUE)/12

mean(FULL$Y_AGE_3, na.rm = TRUE)/12 # 10.95 (SD = .65)
sd(FULL$Y_AGE_3, na.rm = TRUE)/12

mean(FULL$Y_AGE_5, na.rm = TRUE)/12 # 12.05 (SD = .67)
sd(FULL$Y_AGE_5, na.rm = TRUE)/12


##### INCIDENTAL FINDINGS 
table(FULL$INC_FND_1)

# Grouping variable, 1 = some image/incidental finding issue, 0 = no issues
FULL$INC_GRP <- ifelse(FULL$INC_FND_1 == 0 | FULL$INC_FND_1 == 2, 1, 0)

## EARLY LIFE ADVERSITY 

table(FULL$TRA_HISb_1)

# LOSE THE 888 VALUES 

FULL$TRA_HISb_1[FULL$TRA_HISb_1 == 888] <- NA

table(FULL$TRA_HISb_1)


# LOSE "" in SCANNER MODEL & TYPE 

FULL$SCANman_1[FULL$SCANman_1  == ""] <- NA
FULL$SCANman_5[FULL$SCANman_5  == ""] <- NA
FULL$SCANmod_1[FULL$SCANmod_1  == ""] <- NA
FULL$SCANmod_5[FULL$SCANmod_5  == ""] <- NA



################################################################################
############### CENTER NEW VARIABLES 

## AGE
FULL$C_Y_AGE1 <- FULL$Y_AGE_1 - mean(FULL$Y_AGE_1, na.rm = TRUE)
FULL$C_Y_AGE3 <- FULL$Y_AGE_3 - mean(FULL$Y_AGE_3, na.rm = TRUE)
FULL$C_Y_AGE5 <- FULL$Y_AGE_5 - mean(FULL$Y_AGE_5, na.rm = TRUE)

## INTERNALIZING 

FULL$C_CIntSRP_1 <- FULL$CIntSRP_1 - mean(FULL$CIntSRP_1, na.rm = TRUE)
FULL$C_CIntSRY_2 <- FULL$CIntSRY_2 - mean(FULL$CIntSRY_2, na.rm = TRUE)

## LAG 

FULL$C_LAG_1_3 <- FULL$LAG_1_3 - mean(FULL$LAG_1_3, na.rm = TRUE)
FULL$C_LAG_1_5 <- FULL$LAG_1_5 - mean(FULL$LAG_1_5, na.rm = TRUE)
FULL$C_LAG_3_5 <- FULL$LAG_3_5 - mean(FULL$LAG_3_5, na.rm = TRUE)

################################################################################
############### REDUCE DF 

FULL_RED <- FULL %>%
  select(-c(Yagey_W1, FamilyID_W1, Ppensity_W1, siteid,
            "CIntSRP_7", "CIntSRP_9", "CIntSRY_6", "CIntSRY_7", "CIntSRY_9",
            "CORE_VD_6", "CORE_VD_7", "CORE_VD_9", "SCANman_9", "SCANmod_9",
            "TRA_HISf_7", "TRA_HISf_9", "Y_AGE_6", "Y_AGE_7", "Y_AGE_9"))

################################################################################
############### SAVE DAT SHIT 

write.csv(FULL_RED, "ABCD_SLP_DMN_SUI_RR_6.11.25.csv", row.names=FALSE, na="")

prepareMplusData(FULL_RED,"ABCD_SLP_DMN_SUI_RR_6.11.25.dat")


