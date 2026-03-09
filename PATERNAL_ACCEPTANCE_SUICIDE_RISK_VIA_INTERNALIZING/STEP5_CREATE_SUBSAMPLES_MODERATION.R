## LIBRARY
library(MplusAutomation)
library(dplyr)
## SETWD


setwd("/home/cjh37695/ABCD_PROJECTS/SAPPENFIELD_FATHERS/")


## LOAD DF 

df <- read.csv("ABCD_SAPPENFIELD_FATHERS_RR1_FINAL_STRUCTURAL.csv")


## SPLIT FOT T1T9 and T1T5

T1T9 <- subset(df, IMG_PT1T9 == 1)

T1T5 <- subset(df, IMG_PT1T5 == 1)

## SAVE DATA 

write.csv(T1T9, "SAPPENFIELD_T1T9_SUBSAMPLE.csv", row.names = F)
prepareMplusData(T1T9,"SAPPENFIELD_T1T9_SUBSAMPLE.dat", inpfile =T)

write.csv(T1T5, "SAPPENFIELD_T1T5_SUBSAMPLE.csv", row.names = F)
prepareMplusData(T1T5,"SAPPENFIELD_T1T5_SUBSAMPLE.dat", inpfile =T)


######################################################################
############# ADD IN MORE LATENT VARIABLES (MOTHER FATHER INTERCEPTS)

## LOAD IN SAVED FACTORS
PAR <- read.table("PARENTING_LGC.txt", header = FALSE, na.strings = "*")

# Subset the data frame to keep only NUMID, and the latent scores you want

PAR_FACTORS <- PAR[, c("V9", "V11", "V13", "V15", "V18")]


# Give names to the variables
colnames(PAR_FACTORS) <- c("FACC_I", "FACC_S", "MACC_I", "MACC_S", "NUMID")

## CENTER INTERCEPT 

PAR_FACTORS$C_FACC_I <- (PAR_FACTORS$FACC_I - mean(PAR_FACTORS$FACC_I, na.rm = TRUE))
PAR_FACTORS$C_MACC_I <- (PAR_FACTORS$MACC_I - mean(PAR_FACTORS$MACC_I, na.rm = TRUE))
                                                  
######################################
## Merge the latent variables together 

#LOAD IN OG DF 
df <- read.csv("ABCD_SAPPENFIELD_FATHERS_RR1_FINAL_STRUCTURAL.csv")


df_FULL <- full_join(df, PAR_FACTORS, by = "NUMID") 


## SPLIT FOT T1T9 and T1T5

T1T9 <- subset(df_FULL, IMG_PT1T9 == 1)

T1T5 <- subset(df_FULL, IMG_PT1T5 == 1)

## SAVE DATA 

write.csv(df_FULL, "SAPPENFIELD_STRUCTURAL_LAST_FULL.csv", row.names = F)
prepareMplusData(T1T9,"SAPPENFIELD_STRUCTURAL_LAST_FULL.dat", inpfile =T)

write.csv(T1T9, "SAPPENFIELD_STRUCTURAL_LAST_T1T9.csv", row.names = F)
prepareMplusData(T1T9,"SAPPENFIELD_STRUCTURAL_LAST_T1T9.dat", inpfile =T)

write.csv(T1T5, "SAPPENFIELD_STRUCTURAL_LAST_T1T5.csv", row.names = F)
prepareMplusData(T1T5,"SAPPENFIELD_STRUCTURAL_LAST_T1T5.dat", inpfile =T)


################################################################################
###################### ADD IN MORE LATENT VARIABLES 

## LOAD IN SAVED FACTORS
SALSAL <- read.table("AVG_SAL_SAL_T5T9_LCS.txt", header = FALSE, na.strings = "*")
SALAMY <- read.table("AVG_SAL_AMY_T5T9_LCS.txt", header = FALSE, na.strings = "*")
SALNAC <- read.table("AVG_SAL_NAC_T5T9_LCS.txt", header = FALSE, na.strings = "*")
SAL_T1T5 <- read.table("AVG_SAL_SAL_T1T5_LCS.txt", header = FALSE, na.strings = "*")
AMY_T1T5 <- read.table("AVG_SAL_AMY_T1T5_LCS.txt", header = FALSE, na.strings = "*")
NACC_T1T5 <- read.table("AVG_SAL_NAC_T1T5_LCS.txt", header = FALSE, na.strings = "*")
SAL_T1T9 <- read.table("AVG_SAL_SAL_T1T9_LCS.txt", header = FALSE, na.strings = "*")
AMY_T1T9 <- read.table("AVG_SAL_AMY_T1T9_LCS.txt", header = FALSE, na.strings = "*")
NACC_T1T9 <- read.table("AVG_SAL_NAC_T1T9_LCS.txt", header = FALSE, na.strings = "*")

# Subset the data frame to keep only NUMID, and the latent scores you want

SALSALFACTORS <- SALSAL[, c("V7", "V9", "V12")]
SALAMYFACTORS <- SALAMY[, c("V7", "V9", "V12")]
SALNACFACTORS <- SALNAC[, c("V7", "V9", "V12")]
SAL_T1T5_FACTORS <- SAL_T1T5[, c("V7", "V9", "V12")]
AMY_T1T5_FACTORS <- AMY_T1T5[, c("V7", "V9", "V12")]
NACC_T1T5_FACTORS <- NACC_T1T5[, c("V7", "V9", "V12")]
SAL_T1T9_FACTORS <- SAL_T1T9[, c("V7", "V9", "V12")]
AMY_T1T9_FACTORS <- AMY_T1T9[, c("V7", "V9", "V12")]
NACC_T1T9_FACTORS <- NACC_T1T9[, c("V7", "V9", "V12")]

# Give names to the variables
colnames(SALSALFACTORS) <- c("fSALSALt5t9_I", "fSALSALt5t9_C", "NUMID")
colnames(SALAMYFACTORS) <- c("fSALAMYt5t9_I", "fSALAMYt5t9_C", "NUMID")
colnames(SALNACFACTORS) <- c("fSALNACt5t9_I", "fSALNACt5t9_C", "NUMID")
colnames(SAL_T1T5_FACTORS) <- c("fSALSALt1t5_I", "fSALSALt1t5_C", "NUMID")
colnames(AMY_T1T5_FACTORS) <- c("fSALAMYt1t5_I", "fSALAMYt1t5_C", "NUMID")
colnames(NACC_T1T5_FACTORS) <- c("fSALNACt1t5_I", "fSALNACt1t5_C", "NUMID")
colnames(SAL_T1T9_FACTORS) <- c("fSALSALt1t9_I", "fSALSALt1t9_C", "NUMID")
colnames(AMY_T1T9_FACTORS) <- c("fSALAMYt1t9_I", "fSALAMYt1t9_C", "NUMID")
colnames(NACC_T1T9_FACTORS) <- c("fSALNACt1t9_I", "fSALNACt1t9_C", "NUMID")

######################################
## Merge the latent variables together 

ALL_LAT_VARS <- full_join(SALSALFACTORS, SALAMYFACTORS, by = "NUMID") %>%
  full_join(SALNACFACTORS, by = "NUMID")  %>%
  full_join(SAL_T1T5_FACTORS, by = "NUMID")  %>%
  full_join(AMY_T1T5_FACTORS, by = "NUMID")  %>%
  full_join(NACC_T1T5_FACTORS, by = "NUMID")  %>%
  full_join(SAL_T1T9_FACTORS, by = "NUMID")  %>%
  full_join(AMY_T1T9_FACTORS, by = "NUMID")  %>%
  full_join(NACC_T1T9_FACTORS, by = "NUMID") 



names(ALL_LAT_VARS)
## VISUALIZE DISTRIBUTIONS PRIOR TO SCALING 
library(dplyr)
library(tidyr)
library(ggplot2)


vars <- c("fSALSALt5t9_I", "fSALSALt5t9_C",
          "fSALAMYt5t9_I", "fSALAMYt5t9_C",
          "fSALNACt5t9_I", "fSALNACt5t9_C",
          "fSALSALt1t5_I", "fSALSALt1t5_C",
          "fSALAMYt1t5_I", "fSALAMYt1t5_C",
          "fSALNACt1t5_I", "fSALNACt1t5_C",
          "fSALSALt1t9_I", "fSALSALt1t9_C",
          "fSALAMYt1t9_I", "fSALAMYt1t9_C",
          "fSALNACt1t9_I", "fSALNACt1t9_C")

df_long <- ALL_LAT_VARS %>%
  select(all_of(vars)) %>%
  pivot_longer(cols = everything(),
               names_to = "Variable",
               values_to = "Value")

ggplot(df_long, aes(x = Value)) +
  geom_histogram(bins = 30) +
  facet_wrap(~ Variable, scales = "free") +
  theme_minimal()


################################################################################
#################### WINSORIZE AND CENTER IMAGING VARIABLES 

library(psych)

names(ALL_LAT_VARS)

ALL_LAT_VARS$W_fSALSALt5t9_I <- winsor(ALL_LAT_VARS$fSALSALt5t9_I, trim = 0.01)
ALL_LAT_VARS$W_fSALSALt5t9_C <- winsor(ALL_LAT_VARS$fSALSALt5t9_C, trim = 0.01)
ALL_LAT_VARS$W_fSALAMYt5t9_I <- winsor(ALL_LAT_VARS$fSALAMYt5t9_I, trim = 0.01)
ALL_LAT_VARS$W_fSALAMYt5t9_C <- winsor(ALL_LAT_VARS$fSALAMYt5t9_C, trim = 0.01)
ALL_LAT_VARS$W_fSALNACt5t9_I <- winsor(ALL_LAT_VARS$fSALNACt5t9_I, trim = 0.01)
ALL_LAT_VARS$W_fSALNACt5t9_C <- winsor(ALL_LAT_VARS$fSALNACt5t9_C, trim = 0.01)
ALL_LAT_VARS$W_fSALSALt1t5_I <- winsor(ALL_LAT_VARS$fSALSALt1t5_I, trim = 0.01)
ALL_LAT_VARS$W_fSALSALt1t5_C <- winsor(ALL_LAT_VARS$fSALSALt1t5_C, trim = 0.01)
ALL_LAT_VARS$W_fSALAMYt1t5_I <- winsor(ALL_LAT_VARS$fSALAMYt1t5_I, trim = 0.01)
ALL_LAT_VARS$W_fSALAMYt1t5_C <- winsor(ALL_LAT_VARS$fSALAMYt1t5_C, trim = 0.01)
ALL_LAT_VARS$W_fSALNACt1t5_I<- winsor(ALL_LAT_VARS$fSALNACt1t5_I, trim = 0.01)
ALL_LAT_VARS$W_fSALNACt1t5_C <- winsor(ALL_LAT_VARS$fSALNACt1t5_C, trim = 0.01)
ALL_LAT_VARS$W_fSALSALt1t9_I <- winsor(ALL_LAT_VARS$fSALSALt1t9_I, trim = 0.01)
ALL_LAT_VARS$W_fSALSALt1t9_C <- winsor(ALL_LAT_VARS$fSALSALt1t9_C, trim = 0.01)
ALL_LAT_VARS$W_fSALAMYt1t9_I <- winsor(ALL_LAT_VARS$fSALAMYt1t9_I, trim = 0.01)
ALL_LAT_VARS$W_fSALAMYt1t9_C <- winsor(ALL_LAT_VARS$fSALAMYt1t9_C, trim = 0.01)
ALL_LAT_VARS$W_fSALNACt1t9_I<- winsor(ALL_LAT_VARS$fSALNACt1t9_I, trim = 0.01)
ALL_LAT_VARS$W_fSALNACt1t9_C <- winsor(ALL_LAT_VARS$fSALNACt1t9_C, trim = 0.01)

## CENTER INTERCEPTS 

ALL_LAT_VARS$CW_fSALSALt5t9_I <- (ALL_LAT_VARS$W_fSALSALt5t9_I - mean(ALL_LAT_VARS$W_fSALSALt5t9_I, na.rm = TRUE))
ALL_LAT_VARS$CW_fSALAMYt5t9_I <- (ALL_LAT_VARS$W_fSALAMYt5t9_I - mean(ALL_LAT_VARS$W_fSALAMYt5t9_I, na.rm = TRUE))
ALL_LAT_VARS$CW_fSALNACt5t9_I <- (ALL_LAT_VARS$W_fSALNACt5t9_I - mean(ALL_LAT_VARS$W_fSALNACt5t9_I, na.rm = TRUE))
ALL_LAT_VARS$CW_fSALSALt1t5_I <- (ALL_LAT_VARS$W_fSALSALt1t5_I - mean(ALL_LAT_VARS$W_fSALSALt1t5_I, na.rm = TRUE))
ALL_LAT_VARS$CW_fSALAMYt1t5_I <- (ALL_LAT_VARS$W_fSALAMYt1t5_I - mean(ALL_LAT_VARS$W_fSALAMYt1t5_I, na.rm = TRUE))
ALL_LAT_VARS$CW_fSALNACt1t5_I <- (ALL_LAT_VARS$W_fSALNACt1t5_I - mean(ALL_LAT_VARS$W_fSALNACt1t5_I, na.rm = TRUE))
ALL_LAT_VARS$CW_fSALSALt1t9_I <- (ALL_LAT_VARS$W_fSALSALt1t9_I - mean(ALL_LAT_VARS$W_fSALSALt1t9_I, na.rm = TRUE))
ALL_LAT_VARS$CW_fSALAMYt1t9_I <- (ALL_LAT_VARS$W_fSALAMYt1t9_I - mean(ALL_LAT_VARS$W_fSALAMYt1t9_I, na.rm = TRUE))
ALL_LAT_VARS$CW_fSALNACt1t9_I <- (ALL_LAT_VARS$W_fSALNACt1t9_I - mean(ALL_LAT_VARS$W_fSALNACt1t9_I, na.rm = TRUE))


vars <- c("CW_fSALSALt5t9_I", "W_fSALSALt5t9_C",
          "CW_fSALAMYt5t9_I", "W_fSALAMYt5t9_C",
          "CW_fSALNACt5t9_I", "W_fSALNACt5t9_C",
          "CW_fSALSALt1t5_I", "W_fSALSALt1t5_C",
          "CW_fSALAMYt1t5_I", "W_fSALAMYt1t5_C",
          "CW_fSALNACt1t5_I", "W_fSALNACt1t5_C",
          "CW_fSALSALt1t9_I", "W_fSALSALt1t9_C",
          "CW_fSALAMYt1t9_I", "W_fSALAMYt1t9_C",
          "CW_fSALNACt1t9_I", "W_fSALNACt1t9_C")

df_long <- ALL_LAT_VARS %>%
  select(all_of(vars)) %>%
  pivot_longer(cols = everything(),
               names_to = "Variable",
               values_to = "Value")

ggplot(df_long, aes(x = Value)) +
  geom_histogram(bins = 30) +
  facet_wrap(~ Variable, scales = "free") +
  theme_minimal()

#LOAD IN OG DF 
df <- read.csv("SAPPENFIELD_STRUCTURAL_LAST_FULL.csv")

## DROP OLD LATENT CHANGE SCORES 
names(df)

df_RED <- df %>%
  select(-c("SALSALt1t5_I",    "SALSALt1t5_C",
           "SALAMYt1t5_I",    "SALAMYt1t5_C",
           "SALNACt1t5_I",    "SALNACt1t5_C",
           "SALSALt1t9_I",    "SALSALt1t9_C",   
           "SALAMYt1t5_I",    "SALAMYt1t5_C",    
           "SALAMYt1t9_I",   
           "SALAMYt1t9_C",    "SALNACt1t5_I",    
           "SALNACt1t5_C",   
           "SALNACt1t9_I",    "SALNACt1t9_C"))

df_FULL <- full_join(df_RED, ALL_LAT_VARS, by = "NUMID") 

## CREATE MOTHER AND FATHER CHUNK AVERAGES 
df_FULL$DaccCH1 <- (df_FULL$D_acc_M_1 + df_FULL$D_acc_M_3)/2
df_FULL$DaccCH2 <- (df_FULL$D_acc_M_7 + df_FULL$D_acc_M_9)/2
df_FULL$MaccCH1 <- (df_FULL$M_acc_M_1 + df_FULL$M_acc_M_3)/2
df_FULL$MaccCH2 <- (df_FULL$M_acc_M_7 + df_FULL$M_acc_M_9)/2

df_FULL$C_DaccCH1 <- (df_FULL$DaccCH1 - mean(df_FULL$DaccCH1, na.rm = TRUE))*10
df_FULL$C_DaccCH2 <- (df_FULL$DaccCH2 - mean(df_FULL$DaccCH2, na.rm = TRUE))*10
df_FULL$C_MaccCH1 <- (df_FULL$MaccCH1 - mean(df_FULL$MaccCH1, na.rm = TRUE))*10
df_FULL$C_MaccCH2 <- (df_FULL$MaccCH2 - mean(df_FULL$MaccCH2, na.rm = TRUE))*10


## SPLIT FOT T1T9 and T1T5

T1T9 <- subset(df_FULL, IMG_PT1T9 == 1)

T1T5 <- subset(df_FULL, IMG_PT1T5 == 1)

## MORE LIBERAL INCLUSION 
vars <- c("CW_fSALNACt5t9_I", "CW_fSALSALt1t9_I", "CW_fSALAMYt1t9_I")

df_subset <- df %>%
  filter(!if_all(all_of(vars), is.na))

## SAVE DATA 

write.csv(df_FULL, "SAPPENFIELD_STRUCTURAL_LAST_FULL_V3.csv", row.names = F)
prepareMplusData(df_FULL,"SAPPENFIELD_STRUCTURAL_LAST_FULL_V3.dat", inpfile =T)

write.csv(T1T9, "SAPPENFIELD_STRUCTURAL_LAST_T1T9_V3.csv", row.names = F)
prepareMplusData(T1T9,"SAPPENFIELD_STRUCTURAL_LAST_T1T9_V3.dat", inpfile =T)

write.csv(T1T5, "SAPPENFIELD_STRUCTURAL_LAST_T1T5_V3.csv", row.names = F)
prepareMplusData(T1T5,"SAPPENFIELD_STRUCTURAL_LAST_T1T5_V3.dat", inpfile =T)
