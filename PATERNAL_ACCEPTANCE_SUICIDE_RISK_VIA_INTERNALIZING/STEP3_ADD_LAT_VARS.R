#Library

library(dplyr)
library(psych)
library(stringr)
library(purrr)
library(MplusAutomation)

# Set working directory 

setwd("/home/cjh37695/ABCD_PROJECTS/SAPPENFIELD_FATHERS/")

#Library

df <- read.csv("ABCD_SAPPENFIELD_FATHERS_RR1_FINAL_PRE_MEASURE.csv")


################################################################################
#################### CREATE WITHIN-PERSON CENTERED AGE VARIABLES 

age_vars <- grep("^Y_AGE_\\d+$", names(df), value = TRUE)

df <- df %>%
  rowwise() %>%
  mutate(
    AGE_PM = mean(c_across(all_of(age_vars)), na.rm = TRUE)
  ) %>%
  ungroup() %>%
  mutate(
    across(
      all_of(age_vars),
      ~ .x - AGE_PM,
      .names = "PMC_{.col}"
    )
  )

names(df)

################################################################################
#################### ADD IN LATENT VARIABLE SCORES 


## LOAD IN SAVED FACTORS
SAL_T1T9 <- read.table("AVG_SAL_SAL_T1T9_LCS.txt", header = FALSE, na.strings = "*")
AMY_T1T9 <- read.table("AVG_SAL_AMY_T1T9_LCS.txt", header = FALSE, na.strings = "*")
NACC_T1T9 <- read.table("AVG_SAL_NAC_T1T9_LCS.txt", header = FALSE, na.strings = "*")
SAL_T1T5 <- read.table("AVG_SAL_SAL_T1T5_LCS.txt", header = FALSE, na.strings = "*")
AMY_T1T5 <- read.table("AVG_SAL_AMY_T1T5_LCS.txt", header = FALSE, na.strings = "*")
NACC_T1T5 <- read.table("AVG_SAL_NAC_T1T5_LCS.txt", header = FALSE, na.strings = "*")

# Subset the data frame to keep only NUMID, and the latent scores you want

SAL_T1T9_FACTORS <- SAL_T1T9[, c("V10", "V12", "V15")]
AMY_T1T9_FACTORS <- AMY_T1T9[, c("V10", "V12", "V15")]
NACC_T1T9_FACTORS <- NACC_T1T9[, c("V10", "V12", "V15")]
SAL_T1T5_FACTORS <- SAL_T1T5[, c("V10", "V12", "V15")]
AMY_T1T5_FACTORS <- AMY_T1T5[, c("V10", "V12", "V15")]
NACC_T1T5_FACTORS <- NACC_T1T5[, c("V10", "V12", "V15")]

# Give names to the variables
colnames(SAL_T1T9_FACTORS) <- c("SALSALt1t9_I", "SALSALt1t9_C", "NUMID")
colnames(AMY_T1T9_FACTORS) <- c("SALAMYt1t9_I", "SALAMYt1t9_C", "NUMID")
colnames(NACC_T1T9_FACTORS) <- c("SALNACt1t9_I", "SALNACt1t9_C", "NUMID")
colnames(SAL_T1T5_FACTORS) <- c("SALSALt1t5_I", "SALSALt1t5_C", "NUMID")
colnames(AMY_T1T5_FACTORS) <- c("SALAMYt1t5_I", "SALAMYt1t5_C", "NUMID")
colnames(NACC_T1T5_FACTORS) <- c("SALNACt1t5_I", "SALNACt1t5_C", "NUMID")

######################################
## Merge the latent variables together 

ALL_LAT_VARS <- full_join(SAL_T1T5_FACTORS, SAL_T1T9_FACTORS, by = "NUMID") %>%
  full_join(AMY_T1T5_FACTORS, by = "NUMID")  %>%
  full_join(AMY_T1T9_FACTORS, by = "NUMID")  %>%
  full_join(NACC_T1T5_FACTORS, by = "NUMID")  %>%
  full_join(NACC_T1T9_FACTORS, by = "NUMID") 


## VISUALIZE DISTRIBUTIONS PRIOR TO SCALING 
library(dplyr)
library(tidyr)
library(ggplot2)


vars <- c("SALSALt1t5_I", "SALSALt1t5_C", "SALSALt1t9_I", "SALSALt1t9_C",
          "SALAMYt1t5_I", "SALAMYt1t5_C", "SALAMYt1t9_I", "SALAMYt1t9_C",
          "SALNACt1t5_I", "SALNACt1t5_C", "SALNACt1t9_I", "SALNACt1t9_C")

df_long <- ALL_LAT_VARS %>%
  select(all_of(vars)) %>%
  pivot_longer(cols = everything(),
               names_to = "Variable",
               values_to = "Value")

ggplot(df_long, aes(x = Value)) +
  geom_histogram(bins = 30) +
  facet_wrap(~ Variable, scales = "free") +
  theme_minimal()


## LOOKS GOOD 

### MERGE INTO THE FULL DF 

df_FULL <- full_join(df, ALL_LAT_VARS, by = "NUMID")

################################################################################
#################### SAVE THE NEW DF 

write.csv(df_FULL, "ABCD_SAPPENFIELD_FATHERS_RR1_FINAL_STRUCTURAL.csv", row.names = F)
prepareMplusData(df_FULL,"ABCD_SAPPENFIELD_FATHERS_RR1_FINAL_STRUCTURAL.dat", inpfile =T)


