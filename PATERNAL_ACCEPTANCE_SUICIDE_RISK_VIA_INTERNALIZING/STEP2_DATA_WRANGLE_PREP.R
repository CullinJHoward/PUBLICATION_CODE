# Set working directory 


POUNDTOWN <- 1

# set working directory
if (POUNDTOWN == 1) {
  work_dir <- 'C:\\Users\\cjh37695\\Dropbox\\ABCD_FATHERS\\ABCD_NeuroXInternalizing_Suicide\\ANALYSIS\\'
} else {
  work_dir <- 'C:\\Users\\0910h\\Dropbox\\ABCD_FATHERS\\ABCD_NeuroXInternalizing_Suicide\\ANALYSIS\\'
}

setwd(work_dir)

#Library

library(dplyr)
library(MplusAutomation)
library(ggplot2)
# Load data frame 

## BAUER DATA
df <- read.csv("BAUER_COLLABORATION//ABCD_FATHERSxNEURO_SUICIDE_1.9.24.csv")
## OG BROWN DATA 

OG_DF <- read.csv("C:\\Users\\cjh37695\\Dropbox\\ABCD_FATHERS\\ABCD_FATHE.PROSOC.SUIC\\ABCD_FATHER_PROSOC_SUI_RESCHA_2.1.25.csv")

################################################################################
############ BRING IN THE NEEDED DATA FROM BROWN - REDUCE ######################
################################################################################

names(df)

DF_SUB <- df %>%
  select("subID",        "SCANNER", 
           "LamLLvN1",    "RamLLvN1",         
           "LamLvN1",      "RamLvN1",      "SalCin1",      "SalSal1",    
           "CinCin1",      "CinLamy1",     "CinRamy1",         
           "SalLamy1",     "SalRamy1",     "DMNCON1",      "DMNCPN1",
           "DMNDMN1",      "DMNFPN1",      "DMN_NA1",      "DMNSAL1",     
           "DMNVAN1",      "CONFPN1",      "CPNFPN1",         
           "FPNFPN1",      "FPNLnac1",     "FPNRnac1",     "FPNLamy1",    
           "FPNRamy1",     "SALLnac1",     "SALRnac1",     
           "DMNLnac1",     "DMNRnac1",     "DMNLamy1",    "DMNRamy1",    
           "CONLnac1",     "CONRnac1",     "CPNLnac1",        
           "CPNRnac1",     "RS_MOT1",      "MID_MOT1",     
           "INTPRT1",      "INTPRT3",      "INTPRT5",     
           "EXTPRT1",      "EXTPRT3",      "EXTPRT5",
           "INT_M_Y2",     "INT_M_Y3",     "INT_M_Y4",  "INT_M_Y5",     "INT_M_Y6",     "INT_M_Y7",
           "C_FAMTHR", "C_NEIGHTHT", "C_HRSHSES" )

names(OG_DF)

OG_SUB <- OG_DF %>%
  select("subid",  "siteid", "ysex", "familyid", "Par.Edu", "ppensity",
         "yrace",  "marital", "income", "NUMID", "yage1", "yage3", "yage5",
         "SUI_ID_1", "SUI_ID_3", "SUI_ID_5", "SUI_ID_DUM_TOT",  "SUI_ID_CHA", 
         "SLF_HRM_1", "SLF_HRM_3", "SLF_HRM_5", "SLF_HRM_DUM_TOT", "SLF_HRM_CHA",
         "Facc_RS", "Macc_RS", "D_acc_M1", "M_acc_M1")

## RENAME THE ID VARIABLE TO MATCH 

OG_SUB <- OG_SUB %>%
  rename(subID = subid)


## MERGE DFs together 

JOINED_DF <- OG_SUB %>%
  left_join(DF_SUB, by = "subID")


################################################################################
################## RESIDUALIZED CHANGE SCORES ##################################
################################################################################

## INTERNALIZING/EXTERNALIZING PROBLEMS - PARENT REPORT 
names(JOINED_DF)

REG_DF <- JOINED_DF %>%
  filter(complete.cases(INTPRT3, INTPRT1, EXTPRT1, EXTPRT3))

# REGRESS OUT PRIOR TIME POINTS  

EXT_P_REG <- lm(EXTPRT3 ~  EXTPRT1, data = REG_DF)

summary(EXT_P_REG) 

INT_P_REG <- lm(INTPRT3 ~  INTPRT1, data = REG_DF)
summary(INT_P_REG) 

# Extract the residuals

REG_DF$Pext_RS <- resid(EXT_P_REG)
REG_DF$Pint_RS <- resid(INT_P_REG)

# Subset just the unique residual variables I just made to simplify the final merge

P_PSY_RES <- REG_DF %>%
  select(subID, Pint_RS, Pext_RS)


# Merge the residuals dataframe back into the original dataframe
JOINED_DF_PSYres <- JOINED_DF %>%
  left_join(P_PSY_RES, by = "subID")


names(JOINED_DF_PSYres)


## NEUROIMAING SCORES 

names(JOINED_DF_PSYres)

IMG_REG_DF <- JOINED_DF_PSYres %>%
  filter(complete.cases(SalLamy1, SalRamy1, SALLnac1,       
                        SALRnac1, RS_MOT1, ysex))

# REGRESS OUT HEMISPHERES 

SN_RAMY_REG <- lm(SalRamy1 ~  SalLamy1 + RS_MOT1 + factor(ysex), data = IMG_REG_DF)

summary(SN_RAMY_REG) 

SN_LAMY_REG <- lm(SalLamy1 ~  SalRamy1 + RS_MOT1 + factor(ysex), data = IMG_REG_DF)

summary(SN_LAMY_REG)  

SN_RNACC_REG <- lm(SALRnac1 ~  SALLnac1 + RS_MOT1 + factor(ysex), data = IMG_REG_DF)

summary(SN_RNACC_REG) 

SN_LNACC_REG <- lm(SALLnac1 ~  SALRnac1 + RS_MOT1 + factor(ysex), data = IMG_REG_DF)

summary(SN_LNACC_REG)  

# Extract the residuals

IMG_REG_DF$SNLamyRS <- resid(SN_LAMY_REG)
IMG_REG_DF$SNRamyRS <- resid(SN_RAMY_REG)
IMG_REG_DF$SNLnaccRS <- resid(SN_LNACC_REG)
IMG_REG_DF$SNRnaccRS <- resid(SN_RNACC_REG)

# Subset just the unique residual variables I just made to simplify the final merge

IMG_RES <- IMG_REG_DF %>%
  select(subID, SNLamyRS, SNRamyRS, SNLnaccRS, SNRnaccRS)


# Merge the residuals dataframe back into the original dataframe
ALL_RES_DF <- JOINED_DF_PSYres %>%
  left_join(IMG_RES, by = "subID")

names(ALL_RES_DF)
hist(ALL_RES_DF$SALRnac1)
cor(ALL_RES_DF$SALLnac1, ALL_RES_DF$SNRnaccRS, use = "complete.obs")

hist(ALL_RES_DF$SalLamy1)


################################################################################
###################### POMS TRANSFORMATION 
### Time-Invariant POMS

##### POMS over single timepoints #####

SINGLE_POMS <- function(data, var_name) {
  # Check if the variable exists in the data
  if (!(var_name %in% names(data))) {
    stop("Variable not found in the dataset.")
  }
  
  # Convert the variable to numeric if it's not already
  data[[var_name]] <- as.numeric(data[[var_name]])
  
  # Calculate the maximum and minimum values for the variable
  max_val <- max(data[[var_name]], na.rm = TRUE)
  min_val <- min(data[[var_name]], na.rm = TRUE)
  
  # Calculate the POMS score for each observation and multiply by 100
  poms_score <- ifelse(
    is.na(data[[var_name]]), 
    NA, 
    ((data[[var_name]] - min_val) / (max_val - min_val)) * 100
  )
  
  # Add the POMS variable to the dataset with a prefix "p"
  new_var_name <- paste0("p", var_name)
  data[[new_var_name]] <- poms_score
  
  # Return the modified dataframe
  return(data)
}

## Applied TIME-INVARIANT POMS
names(ALL_RES_DF)



#DEMOGRAPHICS
JOINED_DF_PSYres_POMS <- SINGLE_POMS(ALL_RES_DF, var_name = "Par.Edu")
JOINED_DF_PSYres_POMS <- SINGLE_POMS(JOINED_DF_PSYres_POMS, var_name = "income")
#PARENTING MEANS
JOINED_DF_PSYres_POMS <- SINGLE_POMS(JOINED_DF_PSYres_POMS, var_name = "D_acc_M1")
JOINED_DF_PSYres_POMS <- SINGLE_POMS(JOINED_DF_PSYres_POMS, var_name = "M_acc_M1")
#IMAGING - ONLY FOR VARIABLES WITH NO ZERO/NEGATIVE CORRELATIONS
JOINED_DF_PSYres_POMS <- SINGLE_POMS(JOINED_DF_PSYres_POMS, var_name = "SalSal1")
JOINED_DF_PSYres_POMS <- SINGLE_POMS(JOINED_DF_PSYres_POMS, var_name = "CinCin1")
JOINED_DF_PSYres_POMS <- SINGLE_POMS(JOINED_DF_PSYres_POMS, var_name = "DMNDMN1")
JOINED_DF_PSYres_POMS <- SINGLE_POMS(JOINED_DF_PSYres_POMS, var_name = "FPNFPN1")
JOINED_DF_PSYres_POMS <- SINGLE_POMS(JOINED_DF_PSYres_POMS, var_name = "RS_MOT1")
JOINED_DF_PSYres_POMS <- SINGLE_POMS(JOINED_DF_PSYres_POMS, var_name = "MID_MOT1")



##### TIME-VARYING POMS ######


LONG_POMS <- function(data, var_range) {
  # Extracting the variable names based on the given range
  var_names <- names(data)[var_range]
  
  # Extracting the common part of variable names up to the last character (change if different!)
  common_part <- substr(var_names[1], 1, nchar(var_names[1]) - 1)
  
  # Identifying variables that belong to the same measure group
  var_groups <- split(var_names, sapply(var_names, function(x) substr(x, 1, nchar(x) - 1)))
  
  for (var_group in var_groups) {
    # Converting variables to numeric if they are not already
    data[, var_group] <- lapply(data[, var_group], as.numeric)
    
    # Calculating maximum and minimum values across all variables in the group for all participants
    group_max_vals <- apply(data[, var_group], 1, function(x) if (all(is.na(x))) NA else max(x, na.rm = TRUE))
    group_min_vals <- apply(data[, var_group], 1, function(x) if (all(is.na(x))) NA else min(x, na.rm = TRUE))
    
    # Creating new variables with normalized values multiplied by 100
    new_var_names <- paste0("p", var_group)
    
    for (i in 1:length(var_group)) {
      data[[new_var_names[i]]] <- ifelse(
        is.na(data[[var_group[i]]]), 
        NA, 
        ((data[[var_group[i]]] - min(group_min_vals, na.rm = TRUE)) / 
           (max(group_max_vals, na.rm = TRUE) - min(group_min_vals, na.rm = TRUE))) * 100
      )
    }
  }
  
  return(data)
}

## APPLIED LONG-POMS
names(JOINED_DF_PSYres_POMS)

JOINED_DF_PSYres_POMS <- LONG_POMS(JOINED_DF_PSYres_POMS, 
                                  var_range = c(11:13, 66:77))

##### POMPS ON RESIDUALIZED SCORES ######
# This creates a POMPS score ( 0 -100) but retains 0's spot (by centering around original 0 values)
#  because it is meaningful (i.e., no change). 

RESIDUALIZED_POMS <- function(data, var_name) {
  # Check if the variable exists in the data
  if (!(var_name %in% names(data))) {
    stop("Variable not found in the dataset.")
  }
  
  # Convert the variable to numeric if it's not already
  data[[var_name]] <- as.numeric(data[[var_name]])
  
  # Calculate the maximum and minimum values for the variable
  max_val <- max(data[[var_name]], na.rm = TRUE)
  min_val <- min(data[[var_name]], na.rm = TRUE)
  
  # Calculate the POMS score for each observation and multiply by 100
  poms_score <- ifelse(
    is.na(data[[var_name]]), 
    NA, 
    ((data[[var_name]] - min_val) / (max_val - min_val)) * 100
  )
  
  # Compute the POMS value at residual 0 (i.e., the "baseline" value)
  poms_at_zero <- ((0 - min_val) / (max_val - min_val)) * 100
  
  # Adjust the POMS scores by subtracting the value for 0 (no change)
  adjusted_poms_score <- poms_score - poms_at_zero
  
  # Add the adjusted POMS variable to the dataset with a prefix "p"
  new_var_name <- paste0("p", var_name)
  data[[new_var_name]] <- adjusted_poms_score
  
  # Return the modified dataset
  return(data)
}

## APPLIED RESIDUALIZED POMS 
names(JOINED_DF_PSYres_POMS)
# PARENTING 
JOINED_DF_PSYres_POMS <- RESIDUALIZED_POMS(JOINED_DF_PSYres_POMS, var_name = "Macc_RS")
JOINED_DF_PSYres_POMS <- RESIDUALIZED_POMS(JOINED_DF_PSYres_POMS, var_name = "Facc_RS")

#INTERNALIZING/EXTERNALIZING
JOINED_DF_PSYres_POMS <- RESIDUALIZED_POMS(JOINED_DF_PSYres_POMS, var_name = "Pint_RS")
JOINED_DF_PSYres_POMS <- RESIDUALIZED_POMS(JOINED_DF_PSYres_POMS, var_name = "Pext_RS")

#IMAGING VARIABLES
JOINED_DF_PSYres_POMS <- RESIDUALIZED_POMS(JOINED_DF_PSYres_POMS, var_name = "SNLamyRS")
JOINED_DF_PSYres_POMS <- RESIDUALIZED_POMS(JOINED_DF_PSYres_POMS, var_name = "SNRamyRS")
JOINED_DF_PSYres_POMS <- RESIDUALIZED_POMS(JOINED_DF_PSYres_POMS, var_name = "SNLnaccRS")
JOINED_DF_PSYres_POMS <- RESIDUALIZED_POMS(JOINED_DF_PSYres_POMS, var_name = "SNRnaccRS")

JOINED_DF_PSYres_POMS <- RESIDUALIZED_POMS(JOINED_DF_PSYres_POMS, var_name = "LamLLvN1")
JOINED_DF_PSYres_POMS <- RESIDUALIZED_POMS(JOINED_DF_PSYres_POMS, var_name = "RamLLvN1")
JOINED_DF_PSYres_POMS <- RESIDUALIZED_POMS(JOINED_DF_PSYres_POMS, var_name = "LamLvN1")
JOINED_DF_PSYres_POMS <- RESIDUALIZED_POMS(JOINED_DF_PSYres_POMS, var_name = "RamLvN1")

JOINED_DF_PSYres_POMS <- RESIDUALIZED_POMS(JOINED_DF_PSYres_POMS, var_name = "SalCin1")
JOINED_DF_PSYres_POMS <- RESIDUALIZED_POMS(JOINED_DF_PSYres_POMS, var_name = "CinLamy1")
JOINED_DF_PSYres_POMS <- RESIDUALIZED_POMS(JOINED_DF_PSYres_POMS, var_name = "CinRamy1")
JOINED_DF_PSYres_POMS <- RESIDUALIZED_POMS(JOINED_DF_PSYres_POMS, var_name = "SalLamy1")
JOINED_DF_PSYres_POMS <- RESIDUALIZED_POMS(JOINED_DF_PSYres_POMS, var_name = "SalRamy1")
JOINED_DF_PSYres_POMS <- RESIDUALIZED_POMS(JOINED_DF_PSYres_POMS, var_name = "DMNCON1")
JOINED_DF_PSYres_POMS <- RESIDUALIZED_POMS(JOINED_DF_PSYres_POMS, var_name = "DMNCPN1")

JOINED_DF_PSYres_POMS <- RESIDUALIZED_POMS(JOINED_DF_PSYres_POMS, var_name = "DMNFPN1")
JOINED_DF_PSYres_POMS <- RESIDUALIZED_POMS(JOINED_DF_PSYres_POMS, var_name = "DMN_NA1")
JOINED_DF_PSYres_POMS <- RESIDUALIZED_POMS(JOINED_DF_PSYres_POMS, var_name = "DMNSAL1")
JOINED_DF_PSYres_POMS <- RESIDUALIZED_POMS(JOINED_DF_PSYres_POMS, var_name = "DMNVAN1")
JOINED_DF_PSYres_POMS <- RESIDUALIZED_POMS(JOINED_DF_PSYres_POMS, var_name = "CONFPN1")
# ADVERSITY VARIABLES 
JOINED_DF_PSYres_POMS <- RESIDUALIZED_POMS(JOINED_DF_PSYres_POMS, var_name = "C_FAMTHR")
JOINED_DF_PSYres_POMS <- RESIDUALIZED_POMS(JOINED_DF_PSYres_POMS, var_name = "C_HRSHSES")
JOINED_DF_PSYres_POMS <- RESIDUALIZED_POMS(JOINED_DF_PSYres_POMS, var_name = "C_HRSHSES")


################################################################################
####################### SAVE DATA - MPLUS 



names(JOINED_DF_PSYres_POMS)

###### FULL SAMPLE 

## CSV 

write.csv(JOINED_DF_PSYres_POMS, "ABCD_BAUER_COLLAB_5.9.25.csv", row.names=FALSE, na="")

## MPLUS

prepareMplusData(JOINED_DF_PSYres_POMS,"ABCD_BAUER_COLLAB_5.9.25.dat")

