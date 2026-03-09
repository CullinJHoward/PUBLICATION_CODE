## LIBRARY
library(MplusAutomation)
library(dplyr)
## SETWD


setwd("/home/cjh37695/ABCD_PROJECTS/SAPPENFIELD_FATHERS/")


## LOAD DF 

df <- read.csv("SAPPENFIELD_STRUCTURAL_LAST_FULL_V3.csv")


vars <- c("S_rsSAL_Lnacc_9", "S_rsSAL_Rnacc_9", "S_rsSAL_Lamy_9", "S_rsSAL_Ramy_9")

df_subset <- df %>%
  filter(!if_all(all_of(vars), is.na))

T9_IDs <- df_subset$subID

#NEW VARIABLE DESIGNATOR
df$LIB_T9_INC <- as.integer(df$subID %in% T9_IDs)

#MAKE SUBSET
T1T9 <- subset(df, LIB_T9_INC == 1)

table(df$LIB_T9_INC)

##SAVE PRIO TO FACTORING 

write.csv(df, "SAPPENFIELD_STRUCTURAL_LAST_FULL_V4.csv", row.names = F)
prepareMplusData(df,"SAPPENFIELD_STRUCTURAL_LAST_FULL_V4.dat", inpfile =T)

write.csv(T1T9, "SAPPENFIELD_STRUCTURAL_LAST_T1T9_V4.csv", row.names = F)
prepareMplusData(T1T9,"SAPPENFIELD_STRUCTURAL_LAST_T1T9_V4.dat", inpfile =T)

################################################################################
###################### ADD IN MORE LATENT VARIABLES 

## LOAD IN SAVED FACTORS
SAL_T1T9 <- read.table("AVG_SAL_SAL_T1T9_LCS_FULLER.txt", header = FALSE, na.strings = "*")
AMY_T1T9 <- read.table("AVG_SAL_AMY_T1T9_LCS_FULLER.txt", header = FALSE, na.strings = "*")
NACC_T1T9 <- read.table("AVG_SAL_NAC_T1T9_LCS_FULLER.txt", header = FALSE, na.strings = "*")

# Subset the data frame to keep only NUMID, and the latent scores you want


SAL_T1T9_FACTORS <- SAL_T1T9[, c("V7", "V9", "V12")]
AMY_T1T9_FACTORS <- AMY_T1T9[, c("V7", "V9", "V12")]
NACC_T1T9_FACTORS <- NACC_T1T9[, c("V7", "V9", "V12")]

# Give names to the variables
colnames(SAL_T1T9_FACTORS) <- c("FUL_fSALSALt1t9_I", "FUL_fSALSALt1t9_C", "NUMID")
colnames(AMY_T1T9_FACTORS) <- c("FUL_fSALAMYt1t9_I", "FUL_fSALAMYt1t9_C", "NUMID")
colnames(NACC_T1T9_FACTORS) <- c("FUL_fSALNACt1t9_I", "FUL_fSALNACt1t9_C", "NUMID")

######################################
## Merge the latent variables together 

ALL_LAT_VARS <- full_join(SAL_T1T9_FACTORS, AMY_T1T9_FACTORS, by = "NUMID") %>%
  full_join(NACC_T1T9_FACTORS, by = "NUMID") 



names(ALL_LAT_VARS)
## VISUALIZE DISTRIBUTIONS PRIOR TO SCALING 
library(dplyr)
library(tidyr)
library(ggplot2)


vars <- c("FUL_fSALSALt1t9_I", "FUL_fSALSALt1t9_C",
          "FUL_fSALAMYt1t9_I", "FUL_fSALAMYt1t9_C",
          "FUL_fSALNACt1t9_I", "FUL_fSALNACt1t9_C")

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

ALL_LAT_VARS$W_FUL_fSALSALt1t9_I <- winsor(ALL_LAT_VARS$FUL_fSALSALt1t9_I, trim = 0.01)
ALL_LAT_VARS$W_FUL_fSALSALt1t9_C <- winsor(ALL_LAT_VARS$FUL_fSALSALt1t9_C, trim = 0.01)
ALL_LAT_VARS$W_FUL_fSALAMYt1t9_I <- winsor(ALL_LAT_VARS$FUL_fSALAMYt1t9_I, trim = 0.01)
ALL_LAT_VARS$W_FUL_fSALAMYt1t9_C <- winsor(ALL_LAT_VARS$FUL_fSALAMYt1t9_C, trim = 0.01)
ALL_LAT_VARS$W_FUL_fSALNACt1t9_I <- winsor(ALL_LAT_VARS$FUL_fSALNACt1t9_I, trim = 0.01)
ALL_LAT_VARS$W_FUL_fSALNACt1t9_C<- winsor(ALL_LAT_VARS$FUL_fSALNACt1t9_C, trim = 0.01)


## CENTER INTERCEPTS 

ALL_LAT_VARS$CW_FUL_fSALSALt1t9_I <- (ALL_LAT_VARS$W_FUL_fSALSALt1t9_I - mean(ALL_LAT_VARS$W_FUL_fSALSALt1t9_I, na.rm = TRUE))
ALL_LAT_VARS$CW_FUL_fSALAMYt1t9_I <- (ALL_LAT_VARS$W_FUL_fSALAMYt1t9_I - mean(ALL_LAT_VARS$W_FUL_fSALAMYt1t9_I, na.rm = TRUE))
ALL_LAT_VARS$CW_FUL_fSALNACt1t9_I<- (ALL_LAT_VARS$W_FUL_fSALNACt1t9_I - mean(ALL_LAT_VARS$W_FUL_fSALNACt1t9_I, na.rm = TRUE))



vars <- c("CW_FUL_fSALSALt1t9_I", "W_FUL_fSALSALt1t9_C",
          "CW_FUL_fSALAMYt1t9_I", "W_FUL_fSALAMYt1t9_C",
          "CW_FUL_fSALNACt1t9_I", "W_FUL_fSALNACt1t9_C")

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
df <- read.csv("SAPPENFIELD_STRUCTURAL_LAST_FULL_V4.csv")

## ADD IN NEW VARS 

df_FULL <- full_join(df, ALL_LAT_VARS, by = "NUMID") 


## SPLIT FOT T1T9 and T1T5
names(df_FULL)

T1T9_OG <- subset(df_FULL, IMG_PT1T9 == 1)

T1T9_LIB <- subset(df_FULL, LIB_T9_INC == 1)


## SAVE DATA 

write.csv(df_FULL, "SAPPENFIELD_STRUCTURAL_LAST_FULL_V5.csv", row.names = F)
prepareMplusData(df_FULL,"SAPPENFIELD_STRUCTURAL_LAST_FULL_V5.dat", inpfile =T)

write.csv(T1T9_OG, "SAPPENFIELD_STRUCTURAL_LAST_T1T9OG_V5.csv", row.names = F)
prepareMplusData(T1T9_OG,"SAPPENFIELD_STRUCTURAL_LAST_T1T9OG_V5.dat", inpfile =T)

write.csv(T1T9_LIB, "SAPPENFIELD_STRUCTURAL_LAST_T1T9LIB_V5.csv", row.names = F)
prepareMplusData(T1T9_LIB,"SAPPENFIELD_STRUCTURAL_LAST_T1T9LIB_V5.dat", inpfile =T)


