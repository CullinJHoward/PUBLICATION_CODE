#Library

library(dplyr)
library(psych)
library(stringr)
library(purrr)
library(MplusAutomation)

# Set working directory 

setwd("/home/cjh37695/ABCD_PROJECTS/SAPPENFIELD_FATHERS/")

#Library

df <- read.csv("ABCD_SAPPENFIELD_FATHERS_RR_FINAL.csv")


################################################################################
#################### SAMPLE DESCRIPTIVES 
names(df)
## AGE 
describe(df$Y_AGE_1, na.rm = T) # M = 9.95 (.63)
describe(df$Y_AGE_3, na.rm = T) # M = 10.96 (.65)
describe(df$Y_AGE_5, na.rm = T) # M = 12.06 (.67)
describe(df$Y_AGE_7, na.rm = T) # M = 12.95 (.65)
describe(df$Y_AGE_9, na.rm = T) # M = 14.17 (.71)

## SEX 

table(df$Y_SEX)
# 1 (MALE) = 5323 
# 2 (FEMALE) = 4754 (47% female)

## RACE 
table(df$Y_RACE_1)
#1 (HISPANIC) = 2012 (20%)
#2 (WHITE) = 5709 (57%)
#3 (BLACK) = 1155 (11%)
#4 (ASIAN) = 213 (2%)
#5 (OTHER) = 987 (10%)


################################################################################
#################### MEASUREMENT RELIABILITY 

## T1
Facc1_ALPHA <- alpha(df[ , c("D_1acc_1","D_2acc_1","D_3acc_1","D_4acc_1","D_5acc_1")])
print(Facc1_ALPHA) #a = .77

Macc1_ALPHA <- alpha(df[ , c("M_1acc_1","M_2acc_1","M_3acc_1","M_4acc_1","M_5acc_1")])
print(Macc1_ALPHA) # a = .71

## T3
Facc3_ALPHA <- alpha(df[ , c("D_1acc_3","D_2acc_3","D_3acc_3","D_4acc_3","D_5acc_3")])
print(Facc3_ALPHA) #a = .75

Macc3_ALPHA <- alpha(df[ , c("M_1acc_3","M_2acc_3","M_3acc_3","M_4acc_3","M_5acc_3")])
print(Macc3_ALPHA) # a = .72

## T5 NOT AVAILABLE 

## T7
Facc7_ALPHA <- alpha(df[ , c("D_1acc_7","D_2acc_7","D_3acc_7","D_4acc_7","D_5acc_7")])
print(Facc7_ALPHA) #a = .80

Macc7_ALPHA <- alpha(df[ , c("M_1acc_7","M_2acc_7","M_3acc_7","M_4acc_7","M_5acc_7")])
print(Macc7_ALPHA) # a = .82

## T9
Facc9_ALPHA <- alpha(df[ , c("D_1acc_9","D_2acc_9","D_3acc_9","D_4acc_9","D_5acc_9")])
print(Facc9_ALPHA) #a = .83

Macc9_ALPHA <- alpha(df[ , c("M_1acc_9","M_2acc_9","M_3acc_9","M_4acc_9","M_5acc_9")])
print(Macc9_ALPHA) # a = .83

## MEAN ALPHAs

FACC_M <- (.83 + .80 + .75 + .77)/4
FACC_M # .79

MACC_M <- (.83 + .82 + .72 + .71)/4
MACC_M # .77

table(df$IMG_PT1T9)

################################################################################
#################### Rescale variables to improve model convergence

## Here we will increase the scale of the imaging variables (mutiple by 10) to keep
## their scaling, but increase variability and put it on a similar metric as the other
## variables. 

## we will keep parenting on 1 - 3 scale

## We will divide internaliing by 10 to get it on the proper scale 

scale_x_10_prefix_S <- function(df, vars) {
  df %>%
    dplyr::mutate(
      dplyr::across(
        dplyr::all_of(vars),
        ~ ifelse(is.na(.x), NA, .x * 10),
        .names = "S_{.col}"
      )
    )
}

scale_div_10_prefix_S <- function(df, vars) {
  df %>%
    dplyr::mutate(
      dplyr::across(
        dplyr::all_of(vars),
        ~ ifelse(is.na(.x), NA, .x / 10),
        .names = "S_{.col}"
      )
    )
}

## Imaging variables

rsfc_vars <- c(
  "rsMnMOT_1", "rsMnMOT_5", "rsMnMOT_9", 
  "rsSAL_Lamy_1",  "rsSAL_Lamy_5",  "rsSAL_Lamy_9",
  "rsSAL_Lnacc_1", "rsSAL_Lnacc_5", "rsSAL_Lnacc_9",
  "rsSAL_Ltha_1",  "rsSAL_Ltha_5", "rsSAL_Ltha_9",
  "rsSAL_Ramy_1",  "rsSAL_Ramy_5",  "rsSAL_Ramy_9",
  "rsSAL_Rnacc_1", "rsSAL_Rnacc_5", "rsSAL_Rnacc_9",
  "rsSAL_Rtha_1",  "rsSAL_Rtha_5", "rsSAL_Rtha_9",
  "rsSAL_SAL_1",   "rsSAL_SAL_5",   "rsSAL_SAL_9"
)

df_SCALE <- scale_x_10_prefix_S(df, rsfc_vars)


## Psychopathology 

int_vars <- c("CAxDpTP_1", "CAxDpTP_3", "CAxDpTP_5", "CAxDpTP_7", "CAxDpTP_9", "CAxDpTP_11",
              "CIntSTP_1", "CIntSTP_3", "CIntSTP_5", "CIntSTP_7", "CIntSTP_9", "CIntSTP_11",
              "CSomSTP_1", "CSomSTP_3", "CSomSTP_5", "CSomSTP_7", "CSomSTP_9", "CSomSTP_11",
              "CWtDpTP_1", "CWtDpTP_3", "CWtDpTP_5", "CWtDpTP_7", "CWtDpTP_9", "CWtDpTP_11")


df_SCALE <- scale_div_10_prefix_S(df_SCALE, int_vars)

names(df)

### CENTER REDICTORS 

CENTER <- function(x) {
  x - mean(x, na.rm = TRUE)
}

# Applied
names(df_SCALE)

scanner_vars  <- c("rsMnMOT_1", "rsMnMOT_5", "rsMnMOT_9")
age_vars      <- c("Y_AGE_1", "Y_AGE_3", "Y_AGE_5", "Y_AGE_7", "Y_AGE_9", "Y_AGE_11")
par_acc_vars  <- c("D_acc_M_1", "D_acc_M_3", "D_acc_M_7", "D_acc_M_9", "D_acc_M_11",
                   "M_acc_M_1", "M_acc_M_3", "M_acc_M_7", "M_acc_M_9", "M_acc_M_11")
psycho_T1_vars  <- c("S_CIntSTP_1", "S_CSomSTP_1", "S_CAxDpTP_1", "S_CWtDpTP_1")
income_vars  <- c("INCOME6L_1", "INCOME6L_3", "INCOME6L_5", "INCOME6L_7", "INCOME6L_9", "INCOME6L_11")

df_SCALE <- df_SCALE %>%
  mutate(
    across(all_of(scanner_vars), CENTER, .names = "C_{.col}"),
    across(all_of(age_vars), CENTER, .names = "C_{.col}"),
    across(all_of(par_acc_vars), CENTER, .names = "C_{.col}"),
    across(all_of(psycho_T1_vars), CENTER, .names = "C_{.col}"),
    across(all_of(income_vars), CENTER, .names = "C_{.col}")
  )


### AVERAGE CRPSS-HEMISPHERE SCORES AT EACH WAVE 

# SAL - Nacc
df_SCALE$mS_SALNAC1 <- (df_SCALE$S_rsSAL_Lnacc_1  + df_SCALE$S_rsSAL_Rnacc_1)/2
df_SCALE$mS_SALNAC5 <- (df_SCALE$S_rsSAL_Lnacc_5  + df_SCALE$S_rsSAL_Rnacc_5)/2
df_SCALE$mS_SALNAC9 <- (df_SCALE$S_rsSAL_Lnacc_9  + df_SCALE$S_rsSAL_Rnacc_9)/2

#### between hemisphere corelations
cor(df_SCALE$S_rsSAL_Lnacc_1 , df_SCALE$S_rsSAL_Rnacc_1, use = "complete.obs" ) # r = .60
cor(df_SCALE$S_rsSAL_Lnacc_5 , df_SCALE$S_rsSAL_Rnacc_5, use = "complete.obs" ) # r = .61
cor(df_SCALE$S_rsSAL_Lnacc_9 , df_SCALE$S_rsSAL_Rnacc_9, use = "complete.obs" ) # r = .61

# SAL - Amy 
df_SCALE$mS_SALAMY1 <- (df_SCALE$S_rsSAL_Lamy_1  + df_SCALE$S_rsSAL_Ramy_1)/2
df_SCALE$mS_SALAMY5 <- (df_SCALE$S_rsSAL_Lamy_5  + df_SCALE$S_rsSAL_Ramy_5)/2
df_SCALE$mS_SALAMY9 <- (df_SCALE$S_rsSAL_Lamy_9  + df_SCALE$S_rsSAL_Ramy_9)/2

#### between hemisphere corelations
cor(df_SCALE$S_rsSAL_Ramy_1 , df_SCALE$S_rsSAL_Lamy_1, use = "complete.obs" ) # r = .67
cor(df_SCALE$S_rsSAL_Ramy_5 , df_SCALE$S_rsSAL_Lamy_5, use = "complete.obs" ) # r = .66
cor(df_SCALE$S_rsSAL_Ramy_9 , df_SCALE$S_rsSAL_Lamy_9, use = "complete.obs" ) # r = .65


################################################################################
#################### RECOMPUTE SI AND NSSI

# ITEMS 

# SLF_HARM C_7SUIC	You mentioned that in the past, you did some things to hurt yourself, like scratching, cutting, or burning yourself. Were you trying to kill yourself when you did these things?
# SUI ID C_9SUIC	Did you think about how you would do it (even if you had no intention of actually doing it)?
# SUI ID C_10SUIC	Back then, at any point did you have some intention on acting on these thoughts, even if you weren't 100% sure you would do it?
# SUI ID C_11SUIC	Did you think through the details of exactly how you would do it, for instance, decide on a specific place or time?
# SUI ID C_12SUIC	Back then, did you make any preparations for killing yourself?
# SUI ID C_13SUIC	In the past two weeks, did you think through the details of how you would do it, for instance, decide on a specific method, place, or time?
# SLF_HARM C_17SUIC	Sometimes when kids get upset or feel numb, they may do things to hurt themselves, like scratching, cutting, or burning themselves. In the past two weeks, how often have you done any of these things or other things to try to hurt yourself?
# SLF_HARM C_18SUIC	Was there ever a time in the past when you did things to hurt yourself on purpose because you were upset, like cut, scratch or burn yourself?
# C_19SUIC	Was there ever another time in the past when you did things to hurt yourself on purpose, like cut, scratch or burn yourself?
# SUI ID C_20SUIC	In the past two weeks, how often have you wished you were dead or had thoughts that you would be better off dead?
# SUI ID C_21SUIC	Was there ever a time in the past when you often wished you were dead or thought you would be better off dead?
# SUI ID C_22SUIC	Was there ever another time in the past when you often wished you were dead or thought you would be better off dead?
# SUI ID C_25SUIC	Was there ever another time when you thought about wanting to kill yourself?
# SLF_HARM C_29SUIC	You mentioned that in the past 2 weeks you did some things to hurt yourself, like scratching, cutting, or burning yourself. Were you trying to kill yourself when you did these things?
# SUI ID C_31SUIC	You mentioned in the past two weeks you thought about actually wanting to kill yourself. Have you thought about how you would do it (even if you had no intention of actually doing it)?


## UPDATE: We're just going to use the diagnostic categories. No reason to complicate with piecing together symptoms

# df_POMS_LONG_CENT_RED_SUI <- df_POMS_LONG_CENT_RED %>%
#   mutate(SUI_ID_1 = ifelse(rowSums(select(., C_9SUIC1, C_10SUIC1, C_11SUIC1, 
#                                           C_12SUIC1, C_13SUIC1, C_20SUIC1, 
#                                           C_21SUIC1, C_22SUIC1, C_25SUIC1, 
#                                           C_31SUIC1), na.rm = TRUE) > 0, 1, 0))
# df_POMS_LONG_CENT_RED_SUI <- df_POMS_LONG_CENT_RED_SUI %>%
#   mutate(SUI_ID_3 = ifelse(rowSums(select(., C_9SUIC3, C_10SUIC3, C_11SUIC3, 
#                                           C_12SUIC3, C_13SUIC3, C_20SUIC3, 
#                                           C_21SUIC3, C_22SUIC3, C_25SUIC3, 
#                                           C_31SUIC3), na.rm = TRUE) > 0, 1, 0))
# 
# df_POMS_LONG_CENT_RED_SUI <- df_POMS_LONG_CENT_RED_SUI %>%
#   mutate(SUI_ID_5 = ifelse(rowSums(select(., C_9SUIC5, C_10SUIC5, C_11SUIC5, 
#                                           C_12SUIC5, C_13SUIC5, C_20SUIC5, 
#                                           C_21SUIC5, C_22SUIC5, C_25SUIC5, 
#                                           C_31SUIC5), na.rm = TRUE) > 0, 1, 0))
# 
# table(SUI_DF_SUI_ID$SUI_ID_5)
# 
# SUI_DF_SUI_ID$SUI_ID_DUM_TOT <- rowSums(SUI_DF_SUI_ID[, c('SUI_ID_1', 'SUI_ID_3', 'SUI_ID_5')], na.rm = TRUE)
# 
# SUI_DF_SUI_ID$SUI_ID_CHA <- (ifelse(is.na(SUI_DF_SUI_ID$SUI_ID_DUM_TOT), 0, SUI_DF_SUI_ID$SUI_ID_DUM_TOT) - 
#                                ifelse(is.na(SUI_DF_SUI_ID$SUI_ID_1), 0, SUI_DF_SUI_ID$SUI_ID_1))
# 
# SUI_DF_SUI_ID$SUI_ID_CHA <- ifelse(SUI_DF_SUI_ID$SUI_ID_CHA == 2, 1,
#                                    ifelse(SUI_DF_SUI_ID$SUI_ID_CHA == 1,1,0))
# 
# table(SUI_DF_SUI_ID$SUI_ID_CHA)
# 
# ## SELF-HARM
# names(SUI_DF_SUI_ID)
# 
# SUI_DF_SUI_ID <- SUI_DF_SUI_ID %>%
#   mutate(SLF_HRM_1 = ifelse(rowSums(select(., C_7SUIC1, C_17SUIC1, C_18SUIC1, C_19SUIC, C_29SUIC1),
#                                     na.rm = TRUE) > 0, 1, 0))
# SUI_DF_SUI_ID <- SUI_DF_SUI_ID %>%
#   mutate(SLF_HRM_3 = ifelse(rowSums(select(.,C_7SUIC3, C_17SUIC3, C_18SUIC3, C_29SUIC3),
#                                     na.rm = TRUE) > 0, 1, 0))
# 
# SUI_DF_SUI_ID <- SUI_DF_SUI_ID %>%
#   mutate(SLF_HRM_5 = ifelse(rowSums(select(.,C_7SUIC5, C_17SUIC5, C_18SUIC5, C_29SUIC5),
#                                     na.rm = TRUE) > 0, 1, 0))
# 
# 
# 
# SUI_DF_SUI_ID$SLF_HRM_DUM_TOT <- rowSums(SUI_DF_SUI_ID[, c('SLF_HRM_1', 'SLF_HRM_3', 'SLF_HRM_5')], na.rm = TRUE)
# 
# SUI_DF_SUI_ID$SLF_HRM_CHA <- (ifelse(is.na(SUI_DF_SUI_ID$SLF_HRM_DUM_TOT), 0, SUI_DF_SUI_ID$SLF_HRM_DUM_TOT) - 
#                                 ifelse(is.na(SUI_DF_SUI_ID$SLF_HRM_1), 0, SUI_DF_SUI_ID$SLF_HRM_1))
# 
# SUI_DF_SUI_ID$SLF_HRM_CHA <- ifelse(SUI_DF_SUI_ID$SLF_HRM_CHA == 2, 1,
#                                     ifelse(SUI_DF_SUI_ID$SLF_HRM_CHA == 1,1,0))
# 
# names(SUI_DF_SUI_ID)

################################################################################
######## NEW SUICIDAL RISK VARIABLES - CANONICAL CODING FOR SURVIVAL ANALYSIS 
## People wjp experienced the event must "drop out"

#### ACTIVE SUICIDAL IDEATIONS 

# Active Suicidal Ideations (SuIdActT_*)
table(df_SCALE$SuIdActT_1) #0 = 8006, 1 = 48
table(df_SCALE$SuIdActT_3) #0 = 7660, 1 = 54
table(df_SCALE$SuIdActT_5) #0 = 7435, 1 = 83
table(df_SCALE$SuIdActT_7) #0 = 7074, 1 = 120
table(df_SCALE$SuIdActT_9) #0 = 6188, 1 = 165
table(df_SCALE$SuIdActT_11) #0 = 5903, 1 = 26

## Set base true variables
df_SCALE$SV_SuIdActT_1  <- df_SCALE$SuIdActT_1
df_SCALE$SV_SuIdActT_3  <- df_SCALE$SuIdActT_3
df_SCALE$SV_SuIdActT_5  <- df_SCALE$SuIdActT_5
df_SCALE$SV_SuIdActT_7  <- df_SCALE$SuIdActT_7
df_SCALE$SV_SuIdActT_9  <- df_SCALE$SuIdActT_9
df_SCALE$SV_SuIdActT_11 <- df_SCALE$SuIdActT_11

# Code so scores go NA after event is observed 
df_SCALE$SV_SuIdActT_3[df_SCALE$SV_SuIdActT_1 == 1] <- NA

# T5 becomes NA if onset at T1 or T3
df_SCALE$SV_SuIdActT_5[
  df_SCALE$SV_SuIdActT_1 == 1 |
    df_SCALE$SV_SuIdActT_3 == 1
] <- NA

# T7 becomes NA if onset earlier
df_SCALE$SV_SuIdActT_7[
  df_SCALE$SV_SuIdActT_1 == 1 |
    df_SCALE$SV_SuIdActT_3 == 1 |
    df_SCALE$SV_SuIdActT_5 == 1
] <- NA

# T9
df_SCALE$SV_SuIdActT_9[
  df_SCALE$SV_SuIdActT_1 == 1 |
    df_SCALE$SV_SuIdActT_3 == 1 |
    df_SCALE$SV_SuIdActT_5 == 1 |
    df_SCALE$SV_SuIdActT_7 == 1
] <- NA

# T11
df_SCALE$SV_SuIdActT_11[
  df_SCALE$SV_SuIdActT_1 == 1 |
    df_SCALE$SV_SuIdActT_3 == 1 |
    df_SCALE$SV_SuIdActT_5 == 1 |
    df_SCALE$SV_SuIdActT_7 == 1 |
    df_SCALE$SV_SuIdActT_9 == 1
] <- NA

# Active Suicidal Ideation EVENTS [post-transfmormation] (SuIdActT_*)
table(df_SCALE$SV_SuIdActT_1) #0 = 8006, 1 = 48
table(df_SCALE$SV_SuIdActT_3) #0 = 7618, 1 = 51
table(df_SCALE$SV_SuIdActT_5) #0 = 7353, 1 = 77
table(df_SCALE$SV_SuIdActT_7) #0 = 6940, 1 = 98
table(df_SCALE$SV_SuIdActT_9) #0 = 6002, 1 = 128
table(df_SCALE$SV_SuIdActT_11) #0 = 5616, 1 = 14


#### PASSIVE SUICIDAL IDEATIONS 

# PASSIVE SUICIDAL IDEATION (SuIdPasT_*)
table(df_SCALE$SuIdPasT_1) #0 = 7575, 1 = 479
table(df_SCALE$SuIdPasT_3) #0 = 7272, 1 = 442
table(df_SCALE$SuIdPasT_5) #0 = 7069, 1 = 449
table(df_SCALE$SuIdPasT_7) #0 = 6568, 1 = 626
table(df_SCALE$SuIdPasT_9) #0 = 5640, 1 = 713
table(df_SCALE$SuIdPasT_11) #0 = 5820, 1 = 109

## Set base true variables
df_SCALE$SV_SuIdPasT_1  <- df_SCALE$SuIdPasT_1
df_SCALE$SV_SuIdPasT_3  <- df_SCALE$SuIdPasT_3
df_SCALE$SV_SuIdPasT_5  <- df_SCALE$SuIdPasT_5
df_SCALE$SV_SuIdPasT_7  <- df_SCALE$SuIdPasT_7
df_SCALE$SV_SuIdPasT_9  <- df_SCALE$SuIdPasT_9
df_SCALE$SV_SuIdPasT_11 <- df_SCALE$SuIdPasT_11

# Code so scores go NA after event is observed
df_SCALE$SV_SuIdPasT_3[df_SCALE$SV_SuIdPasT_1 == 1] <- NA

# T5 becomes NA if onset at T1 or T3
df_SCALE$SV_SuIdPasT_5[
  df_SCALE$SV_SuIdPasT_1 == 1 |
    df_SCALE$SV_SuIdPasT_3 == 1
] <- NA

# T7 becomes NA if onset earlier
df_SCALE$SV_SuIdPasT_7[
  df_SCALE$SV_SuIdPasT_1 == 1 |
    df_SCALE$SV_SuIdPasT_3 == 1 |
    df_SCALE$SV_SuIdPasT_5 == 1
] <- NA

# T9
df_SCALE$SV_SuIdPasT_9[
  df_SCALE$SV_SuIdPasT_1 == 1 |
    df_SCALE$SV_SuIdPasT_3 == 1 |
    df_SCALE$SV_SuIdPasT_5 == 1 |
    df_SCALE$SV_SuIdPasT_7 == 1
] <- NA

# T11
df_SCALE$SV_SuIdPasT_11[
  df_SCALE$SV_SuIdPasT_1 == 1 |
    df_SCALE$SV_SuIdPasT_3 == 1 |
    df_SCALE$SV_SuIdPasT_5 == 1 |
    df_SCALE$SV_SuIdPasT_7 == 1 |
    df_SCALE$SV_SuIdPasT_9 == 1
] <- NA

# PASSIVE SUICIDAL IDEATION EVENTS (SuIdPasT_*)
table(df_SCALE$SV_SuIdPasT_1) #0 = 7575, 1 = 479
table(df_SCALE$SV_SuIdPasT_3) #0 = 6938, 1 = 316
table(df_SCALE$SV_SuIdPasT_5) #0 = 6513, 1 = 264
table(df_SCALE$SV_SuIdPasT_7) #0 = 5902, 1 = 343
table(df_SCALE$SV_SuIdPasT_9) #0 = 4891, 1 = 331
table(df_SCALE$SV_SuIdPasT_11) #0 = 4591, 1 = 37
331/(4891+331)
#### NONSUICIDAL SELF-INJURY 

# NONSUICIDAL SELF-INJURY (NSSIT_*)
table(df_SCALE$NSSIT_1) #0 = 7628, 1 = 426
table(df_SCALE$NSSIT_3) #0 = 7371, 1 = 343
table(df_SCALE$NSSIT_5) #0 = 7236, 1 = 282
table(df_SCALE$NSSIT_7) #0 = 6685, 1 = 509
table(df_SCALE$NSSIT_9) #0 = 5754, 1 = 599
table(df_SCALE$NSSIT_11) #0 = 5759, 1 = 170

## Set base true variables
df_SCALE$SV_NSSIT_1  <- df_SCALE$NSSIT_1
df_SCALE$SV_NSSIT_3  <- df_SCALE$NSSIT_3
df_SCALE$SV_NSSIT_5  <- df_SCALE$NSSIT_5
df_SCALE$SV_NSSIT_7  <- df_SCALE$NSSIT_7
df_SCALE$SV_NSSIT_9  <- df_SCALE$NSSIT_9
df_SCALE$SV_NSSIT_11 <- df_SCALE$NSSIT_11

# Code so scores go NA after event is observed
df_SCALE$SV_NSSIT_3[df_SCALE$SV_NSSIT_1 == 1] <- NA

# T5 becomes NA if onset at T1 or T3
df_SCALE$SV_NSSIT_5[
  df_SCALE$SV_NSSIT_1 == 1 |
    df_SCALE$SV_NSSIT_3 == 1
] <- NA

# T7 becomes NA if onset earlier
df_SCALE$SV_NSSIT_7[
  df_SCALE$SV_NSSIT_1 == 1 |
    df_SCALE$SV_NSSIT_3 == 1 |
    df_SCALE$SV_NSSIT_5 == 1
] <- NA

# T9
df_SCALE$SV_NSSIT_9[
  df_SCALE$SV_NSSIT_1 == 1 |
    df_SCALE$SV_NSSIT_3 == 1 |
    df_SCALE$SV_NSSIT_5 == 1 |
    df_SCALE$SV_NSSIT_7 == 1
] <- NA

# T11
df_SCALE$SV_NSSIT_11[
  df_SCALE$SV_NSSIT_1 == 1 |
    df_SCALE$SV_NSSIT_3 == 1 |
    df_SCALE$SV_NSSIT_5 == 1 |
    df_SCALE$SV_NSSIT_7 == 1 |
    df_SCALE$SV_NSSIT_9 == 1
] <- NA


# NONSUICIDAL SELF-INJURY EVENTS (NSSIT_*)
table(df_SCALE$SV_NSSIT_1) #0 = 7628, 1 = 426
table(df_SCALE$SV_NSSIT_3) #0 = 7063, 1 = 250
table(df_SCALE$SV_NSSIT_5) #0 = 6727, 1 = 167
table(df_SCALE$SV_NSSIT_7) #0 = 6139, 1 = 305
table(df_SCALE$SV_NSSIT_9) #0 = 5104, 1 = 313
table(df_SCALE$SV_NSSIT_11) #0 = 4735, 1 = 83


#### REMOVE NEEDLESS VARIABLES

names(df_SCALE)

df_SCALE_RED <- df_SCALE %>%
  select(-c("NSSIT_1",         "NSSIT_3",         "NSSIT_5",         "NSSIT_7",        
            "NSSIT_9",         "NSSIT_11",        "NSSIpr_1",        "NSSIpr_3",        "NSSIpr_5",        "NSSIpr_7",        "NSSIpr_9",       
            "NSSIpr_11",       "NSSIpt_1",        "NSSIpt_3",        "NSSIpt_5",        "NSSIpt_7",        "NSSIpt_9",        "NSSIpt_11",      
            "SuIdActT_1",      "SuIdActT_3",      "SuIdActT_5",      "SuIdActT_7",      "SuIdActT_9",      "SuIdActT_11",     "SuIdActpr_1",    
            "SuIdActpr_3",     "SuIdActpr_5",     "SuIdActpr_7",     "SuIdActpr_9",     "SuIdActpr_11",    "SuIdActpt_1",     "SuIdActpt_3",    
            "SuIdActpt_5",     "SuIdActpt_7",     "SuIdActpt_9",     "SuIdActpt_11",    "SuIdPasT_1",      "SuIdPasT_3",      "SuIdPasT_5",     
            "SuIdPasT_7",      "SuIdPasT_9",      "SuIdPasT_11",     "SuIdPaspr_1",     "SuIdPaspr_3",     "SuIdPaspr_5",     "SuIdPaspr_7",    
            "SuIdPaspr_9",     "SuIdPaspr_11",    "SuIdPaspt_1",     "SuIdPaspt_3",     "SuIdPaspt_5",     "SuIdPaspt_7",     "SuIdPaspt_9",    
            "SuIdPaspt_11",    "rsSAL_Lamy_1",    "rsSAL_Lamy_5",    "rsSAL_Lamy_9",    "rsSAL_Lnacc_1",   "rsSAL_Lnacc_5",  
            "rsSAL_Lnacc_9",   "rsSAL_Ltha_1",    "rsSAL_Ltha_5",    "rsSAL_Ltha_9",    "rsSAL_Ramy_1",    "rsSAL_Ramy_5",    "rsSAL_Ramy_9",   
            "rsSAL_Rnacc_1",   "rsSAL_Rnacc_5",   "rsSAL_Rnacc_9",   "rsSAL_Rtha_1",    "rsSAL_Rtha_5",    "rsSAL_Rtha_9",    "rsSAL_SAL_1",    
            "rsSAL_SAL_5",     "rsSAL_SAL_9",     "CAxDpTP_1",       "CAxDpTP_3",       "CAxDpTP_5",      
            "CAxDpTP_7",       "CAxDpTP_9",       "CAxDpTP_11",      "CIntSTP_1",       "CIntSTP_3",       "CIntSTP_5",       "CIntSTP_7",      
            "CIntSTP_9",       "CIntSTP_11",      "CSomSTP_1",       "CSomSTP_3",       "CSomSTP_5",       "CSomSTP_7",       "CSomSTP_9",      
            "CSomSTP_11",      "CWtDpTP_1",       "CWtDpTP_3",       "CWtDpTP_5",       "CWtDpTP_7",       "CWtDpTP_9",       "CWtDpTP_11"  ))

names(df_SCALE_RED)



################################################################################
#################### SAVE DATA FOR SEM 


## SAVE BOTH DFs
write.csv(df_SCALE, "ABCD_SAPPENFIELD_FATHERS_RR1_FINAL_PRE_MEASURE.csv", row.names = F)


# MPLUS PREP 
prepareMplusData(df_SCALE_RED,"ABCD_SAPPENFIELD_FATHERS_RR1_FINAL_PRE_MEASURE.dat", inpfile =T)

################################################################################
#################### WAVE MISSINGNESS 

library(dplyr)

waves <- c("_1", "_3", "_5", "_7", "_9", "_11")

wave_qc <- lapply(waves, function(w) {
  
  # Identify variables ending in the wave suffix
  vars_w <- grep(paste0(w, "$"), names(df_SCALE_RED), value = TRUE)
  
  cat("\n=============================\n")
  cat("Wave suffix:", w, "\n")
  cat("Variables identified:\n")
  print(vars_w)
  cat("Number of variables:", length(vars_w), "\n")
  
  if (length(vars_w) == 0) {
    return(data.frame(
      wave = w,
      n_vars = 0,
      n_all_na = NA,
      n_non_na = NA,
      total_n = nrow(df_SCALE_RED)
    ))
  }
  
  # Rows that are ALL NA for this wave
  all_na <- df_SCALE_RED %>%
    select(all_of(vars_w)) %>%
    apply(1, function(x) all(is.na(x)))
  
  data.frame(
    wave = w,
    n_vars = length(vars_w),
    n_all_na = sum(all_na),
    n_non_na = sum(!all_na),
    total_n = nrow(df_SCALE_RED)
  )
})

wave_qc <- bind_rows(wave_qc)
wave_qc

N <- 6232
NO <- 4735
YES <- 83

N - (NO+YES)

