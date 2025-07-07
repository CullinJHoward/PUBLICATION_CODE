
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
####################### LOAD IN ALL SURVEY AND LATENT VARIABLES 

df <- read.csv ("ARCHIVE\\DISSERTATION_REDUCED_SURVEY_AND_LAT_3.7.25.csv")

names(df)


################################################################################
####################### POMPS CONVERT OUTCOMES AND COVARIATES 

##### TIME-VARYING POMPS ######

# USED-FOR MEAN LEVEL SCORES 

LONG_POMPS <- function(data, var_range) {
  # Extracting the variable names based on the given range
  var_names <- names(data)[var_range]
  
  # Identifying variables that belong to the same measure group
  var_groups <- split(var_names, sapply(var_names, function(x) substr(x, 1, nchar(x) - 1)))
  
  # Print identified groups
  cat("Identified groups for longitudinal POMPS values:\n")
  print(var_groups)
  
  for (var_group in var_groups) {
    # Converting variables to numeric if they are not already
    data[, var_group] <- lapply(data[, var_group], as.numeric)
    
    # Calculating maximum and minimum values across all variables in the group for all participants
    group_max_vals <- apply(data[, var_group], 1, function(x) if (all(is.na(x))) NA else max(x, na.rm = TRUE))
    group_min_vals <- apply(data[, var_group], 1, function(x) if (all(is.na(x))) NA else min(x, na.rm = TRUE))
    
    # Creating new variables with normalized values multiplied by 100
    new_var_names <- paste0("p", var_group)
    
    for (i in seq_along(var_group)) {
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

## APPLIED LONG-POMPS
names(df)

df_LONGPOMPS <- LONG_POMPS(df, var_range = c(49:52,   #P-FACTOR 
                                             57:58,   # SOCIAL ACCEPTANCE
                                             61:64    #EMPATHIC CONCERN AND PERSPECTIVE TAKING 
                                             ))

##### POMPS ON RESIDUALIZED CHANGE SCORE ######

## CREATES POMPS SCORES, BUT IDENTIFIES 0 (I.E., NO CHANGE) AND RETAINS ITS LOCATION

RECENTERED_POMPS <- function(data, var_name, center_value) {
  # Check if the variable exists in the data
  if (!(var_name %in% names(data))) {
    stop("Variable not found in the dataset.")
  }
  
  # Convert the variable to numeric if it's not already
  data[[var_name]] <- as.numeric(data[[var_name]])
  
  # Calculate the maximum and minimum values for the variable
  max_val <- max(data[[var_name]], na.rm = TRUE)
  min_val <- min(data[[var_name]], na.rm = TRUE)
  
  # Calculate the POMPS score for each observation and multiply by 100
  poms_score <- ifelse(
    is.na(data[[var_name]]), 
    NA, 
    ((data[[var_name]] - min_val) / (max_val - min_val)) * 100
  )
  
  # Compute the POMPS value at the specified center_value
  poms_at_center <- ((center_value - min_val) / (max_val - min_val)) * 100
  
  # Adjust the POMPS scores by centering on the specified value
  adjusted_poms_score <- poms_score - poms_at_center
  
  # Add the adjusted POMPS variable to the dataset with a prefix "p"
  new_var_name <- paste0("C_p", var_name)
  data[[new_var_name]] <- adjusted_poms_score
  
  # Return the modified dataset
  return(data)
}

## APPLIED RECENTERED POMPS ON RESIDUALIZED CHANGE SCORES 

names(df_LONGPOMPS)

df_LONGPOMPS_RESPOMP <- RECENTERED_POMPS(df_LONGPOMPS, var_name = "Pfac_PLC", center_value = 0) # PARENT P-FACTOR
df_LONGPOMPS_RESPOMP <- RECENTERED_POMPS(df_LONGPOMPS_RESPOMP, var_name = "Pfac_YLC", center_value = 0) # CHILD P-FACTOR
df_LONGPOMPS_RESPOMP <- RECENTERED_POMPS(df_LONGPOMPS_RESPOMP, var_name = "SOACC_LC", center_value = 0) # SOCIAL ACCEPTANCE
df_LONGPOMPS_RESPOMP <- RECENTERED_POMPS(df_LONGPOMPS_RESPOMP, var_name = "PERTK_LC", center_value = 0) # PERSPECTIVE TAKING 
df_LONGPOMPS_RESPOMP <- RECENTERED_POMPS(df_LONGPOMPS_RESPOMP, var_name = "EMPCN_LC", center_value = 0) # EMPATHIC CONCERN 

## APPLIED RECENTERED POMPS ON INTERCPET VALUES - CENTERED AT THE AVERAGE INTERCEPT VALUE

df_LONGPOMPS_RESPOMP_INT <- RECENTERED_POMPS(df_LONGPOMPS_RESPOMP, var_name = "Pfac_PLI", center_value = 1.46) # PARENT P-FACTOR
df_LONGPOMPS_RESPOMP_INT <- RECENTERED_POMPS(df_LONGPOMPS_RESPOMP_INT, var_name = "Pfac_YLI", center_value = 3.59) # CHILD P-FACTOR
df_LONGPOMPS_RESPOMP_INT <- RECENTERED_POMPS(df_LONGPOMPS_RESPOMP_INT, var_name = "SOACC_LI", center_value = 3.46) # SOCIAL ACCEPTANCE
df_LONGPOMPS_RESPOMP_INT <- RECENTERED_POMPS(df_LONGPOMPS_RESPOMP_INT, var_name = "PERTK_LI", center_value = 3.49) # PERSPECTIVE TAKING 
df_LONGPOMPS_RESPOMP_INT <- RECENTERED_POMPS(df_LONGPOMPS_RESPOMP_INT, var_name = "EMPCN_LI", center_value = 3.99) # EMPATHIC CONCERN 


## APPLIED RECENTERED POMPS ON COVARIATES - ALL CENTERED AT THEIR MEANS
names(df_LONGPOMPS_RESPOMP_INT)

df_LONGPOMPS_RESPOMP_INT_COV <- RECENTERED_POMPS(df_LONGPOMPS_RESPOMP_INT, var_name = "P_EDU", center_value = 4.76) # Parent education 
df_LONGPOMPS_RESPOMP_INT_COV <- RECENTERED_POMPS(df_LONGPOMPS_RESPOMP_INT_COV, var_name = "C_AGE", center_value = 12.89) #Child age, W1
df_LONGPOMPS_RESPOMP_INT_COV <- RECENTERED_POMPS(df_LONGPOMPS_RESPOMP_INT_COV, var_name = "INCOME_W1", center_value = 64.44) #Family income 


## REORDER VARIABLES
names(df_LONGPOMPS_RESPOMP_INT_COV)

FINAL_DF <- df_LONGPOMPS_RESPOMP_INT_COV %>%
  select(c(ID, P_SEX, C_SEX , P_EDU, C_RACE, C_AGE, MAR_STAT, INCOME_W1,
           C_pP_EDU, C_pC_AGE, C_pINCOME_W1,
           CTQCW1_T, P_CTQ_T, FOOD_INS, SUBSEScm, SUBSESus, 
           pTHR_Y, pDEP_Y, pUNP_Y, pCUM_Y, pTHR_P, pDEP_P, pCUM_P,     
           C_pTHR_Y, Q_C_pTHR_Y, C_pDEP_Y, Q_C_pDEP_Y, C_pUNP_Y, Q_C_pUNP_Y, C_pCUM_Y, Q_C_pCUM_Y,
           C_pTHR_P, Q_C_pTHR_P, C_pDEP_P, Q_C_pDEP_P, C_pCUM_P, Q_C_pCUM_P, 
           Pfac_P1, Pfac_P2, Pfac_PLI, Pfac_PLC, pPfac_P1, pPfac_P2, C_pPfac_PLI, C_pPfac_PLC,
           Pfac_Y1, Pfac_Y2, Pfac_YLI, Pfac_YLC, pPfac_Y1, pPfac_Y2, C_pPfac_YLC, C_pPfac_YLI,
           SOACC_1, SOACC_2, SOACC_LI, SOACC_LC, pSOACC_1, pSOACC_2, C_pSOACC_LC, C_pSOACC_LI,
           PERTK_1, PERTK_2, PERTK_LI, PERTK_LC, pPERTK_1, pPERTK_2, C_pPERTK_LC, C_pPERTK_LI,
           EMPCN_1, EMPCN_2, EMPCN_LI, EMPCN_LC, pEMPCN_1, pEMPCN_2, C_pEMPCN_LC, C_pEMPCN_LI))

## SAVE DF 

write.csv(FINAL_DF, "DISS_REDUCED_SURVEY_POMPS_3.8.25.csv", row.names=FALSE, na="")


