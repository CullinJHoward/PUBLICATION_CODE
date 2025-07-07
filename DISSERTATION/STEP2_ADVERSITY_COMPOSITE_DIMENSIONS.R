
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
####################### COMPUTE Parnet and Child ELA DIMESNISONS

df <- read.csv("SURVEY_DATA\\DORRY_CHILD_DIM_ADV.csv")

# REMOVE DUPLICATE IDs (but fill in any missing values)

df <- df %>%
  group_by(ID) %>%
  summarise(across(everything(), ~ na.omit(.)[1])) %>%
  ungroup()

## CREATE A USER IDENTIFIED RANGE POMPS FUNCTION 
# THIS IS WHEN THERE IS A SET (AND MEANINGFUL) MINIMUM AND MAXIMUM
# BECAUSE I AM SUMMING ACROSS SCALES, I WANT THEM TO HOLD THEIR WEIGHTING, SO I
# AM ENTERING THE SCALE MINIMUM AND MAXIMUM AND THEN COMPUTING ITS PERCENTILE BASED 
# ON THESE VALUES. THIS IS DIFFERENT THAN THE RELATIVE POMPS SCORE I HAVE PREVIOUSLY USED
# WHERE THE BOUNDS ARE THE OBSERVED RANGE OF THE DATA 


## USAGE: DF_POMS <- SINGLE_POMP_SET(df, "CRIM_ACT", min_val = 1, max_val = 20)
# MIN IS THE LOWEST POSSIBLE SCORE, MAX IS THE MAX POSSIBLE (NOT OBSERVED) SCORE

SINGLE_POMP_SET <- function(data, var_name, min_val, max_val) {
  # Check if the variable exists in the data
  if (!(var_name %in% names(data))) {
    stop("Variable not found in the dataset.")
  }
  
  # Convert the variable to numeric if it's not already
  data[[var_name]] <- as.numeric(data[[var_name]])
  
  # Check if min_val and max_val are valid
  if (min_val >= max_val) {
    stop("min_val must be less than max_val.")
  }
  
  # Calculate the POMPS score for each observation using the provided min and max
  poms_score <- ifelse(
    is.na(data[[var_name]]), 
    NA, 
    (((data[[var_name]] - min_val) / (max_val - min_val))) * 100
  )
  
  # Add the POMPS variable to the dataset with a prefix "p"
  new_var_name <- paste0("p", var_name)
  data[[new_var_name]] <- poms_score
  
  # Return the modified dataframe
  return(data)
}

## SET UP RECODING FUNCTION 

RECODE_CTS <- function(data, var_name) {
  # Check if the variable exists in the data
  if (!(var_name %in% names(data))) {
    stop("Variable not found in the dataset.")
  }
  
  # Apply recoding using nested ifelse
  data[[var_name]] <- ifelse(is.na(data[[var_name]]), NA,
                             ifelse(data[[var_name]] == 0, 0,
                                    ifelse(data[[var_name]] == 7, 1,
                                           ifelse(data[[var_name]] == 1, 2,
                                                  ifelse(data[[var_name]] == 2, 3,
                                                         ifelse(data[[var_name]] == 3, 4,
                                                                ifelse(data[[var_name]] == 4, 5,
                                                                       ifelse(data[[var_name]] == 5, 6,
                                                                              ifelse(data[[var_name]] == 6, 7, NA)))))))))

# Return modified data
return(data)
}
  

################ CHILD ELA

# Items were categorized into their most relevant dimension
# Items getting at the same underlying construct (but having novel elements) were parceled using summing
# Each dimensional subconstruct was POMP-scored, keeping their rank order relative to objective scoring, but
# placing each on a comparable scale 
# subdimensional subconstruct scores were summed. 
# the final dimensional sum score was POMPs scored, putting the highest observed value as the max. 


################################ THREAT ########################################


######## PARCEL LIKE ITEMS TOGETHER THROUGH SUMMING 

########## CRIMINAL ACTIVITY 

df$CRIM_ACT <- apply(df[, c("NER5PW1", "NER9PW1")], 1, function(x) {
  # If any value is missing, set total to NA
  if (any(is.na(x))) {
    return(NA)
  } else {
    # Otherwise, sum the non-missing values
    return(sum(x, na.rm = TRUE))
  }
})


min(df$CRIM_ACT, na.rm = TRUE) #2


########## WRECKLESS ACTIVITY 

df$WRECK_ACT <- apply(df[, c("NER13PW1", "NER14PW1", "NER6PW1")], 1, function(x) {
  # Count the number of missing values
  missing_count <- sum(is.na(x))
  
  # If more than 1 value is missing, set total to NA
  if (missing_count > 1) {
    return(NA)
  } else {
    # Otherwise, sum the non-missing values
    return(sum(x, na.rm = TRUE))
  }
})


min(df$WRECK_ACT, na.rm = TRUE)#2

########## NEIGHBORHOOD SAFETY

df$UNSAFE_NEI <- apply(df[, c("NBC3PW1", "NBC5PW1", "NBC9PW1")], 1, function(x) {
  # Count the number of missing values
  missing_count <- sum(is.na(x))
  
  # If more than 1 value is missing, set total to NA
  if (missing_count > 1) {
    return(NA)
  } else {
    # Otherwise, sum the non-missing values
    return(sum(x, na.rm = TRUE))
  }
})


min(df$UNSAFE_NEI, na.rm = TRUE)#3

########## PHYSICAL DISCIPLINE 

## RECODE VARIABLES 
df <- RECODE_CTS(df, "CTS3PW1")
df <- RECODE_CTS(df, "CTS4PW1")
df <- RECODE_CTS(df, "CTS17PW1")

df$PHY_DISP <- apply(df[, c("CTS3PW1", "CTS4PW1", "CTS17PW1")], 1, function(x) {
  # Count the number of missing values
  missing_count <- sum(is.na(x))
  
  # If more than 1 value is missing, set total to NA
  if (missing_count > 1) {
    return(NA)
  } else {
    # Otherwise, sum the non-missing values
    return(sum(x, na.rm = TRUE))
  }
})


min(df$PHY_DISP, na.rm = TRUE)#0

########## PSYCHOLOGICAL AGGRESSION 

## RECODE VARIABLES 
df <- RECODE_CTS(df, "CTS6PW1")
df <- RECODE_CTS(df, "CTS10PW1")

df$PSY_AGG <- apply(df[, c("CTS6PW1", "CTS10PW1")], 1, function(x) {
  # If there is any missing value, set the total to NA
  if (any(is.na(x))) {
    return(NA)
  } else {
    # Otherwise, sum the values
    return(sum(x))
  }
})


min(df$PSY_AGG, na.rm = TRUE)#0

######## CREATE POMPS SCORES 

DF_POMPS <- SINGLE_POMP_SET(df, "CRIM_ACT", min_val = 2, max_val = 20) # CRIMAL ACTIVITY (two items, scale of 1-10)
DF_POMPS <- SINGLE_POMP_SET(DF_POMPS, "WRECK_ACT", min_val = 2, max_val = 30) #WRECKLESS ACTIVITY (three items, scale of 1-10)
DF_POMPS <- SINGLE_POMP_SET(DF_POMPS, "UNSAFE_NEI", min_val = 3, max_val = 30) #UNSAFE NEIGHBORHOOD (three items, scale of 1-10)

DF_POMPS <- SINGLE_POMP_SET(DF_POMPS, "PHY_DISP", min_val = 0, max_val = 21) #PARENT PSYCHOLOGICAL AGGRESSION (three itemes, scale of 0 - 7)
DF_POMPS <- SINGLE_POMP_SET(DF_POMPS, "PSY_AGG", min_val = 0, max_val = 14) #PHYSICAL DISCIPLINE (two itemes, scale of 0 - 7)

DF_POMPS <- SINGLE_POMP_SET(DF_POMPS, "CTQCW1_PA", min_val = 5, max_val = 25) #PHYSICAL ABUSE
DF_POMPS <- SINGLE_POMP_SET(DF_POMPS, "CTQCW1_SA", min_val = 5, max_val = 25) #SEXUAL ABUSE
DF_POMPS <- SINGLE_POMP_SET(DF_POMPS, "CTQCW1_EA", min_val = 5, max_val = 25) #EMOTIONAL ABUSE

########CREATE A TOTAL THREAT SCORE 

DF_POMPS$THR_Y <- apply(DF_POMPS[, c("pCRIM_ACT", "pWRECK_ACT", "pUNSAFE_NEI", "pPHY_DISP", "pPSY_AGG",
                                        "pCTQCW1_PA", "pCTQCW1_SA", "pCTQCW1_EA")], 1, function(x) {
                                          # Count the number of missing values
                                          missing_count <- sum(is.na(x))
                                          
                                          # If more than 4 values are missing, set total to NA
                                          if (missing_count > 4) {
                                            return(NA)
                                          } else {
                                            # Otherwise, sum the non-missing values
                                            return(sum(x, na.rm = TRUE))
                                          }
                                        })

## MIN IS ZERO (SINCE THEY ARE POMPS SCORES, A ZERO ON ALL IS POSSIBLE, JUST UNLIKELY)
## MAX IS 800 (SINCE POMPS SCORES, RANGE IS 0 - 100. 8 X 100 = 800)

DF_POMPS <- SINGLE_POMP_SET(DF_POMPS, "THR_Y", min_val = 0, max_val = 800) #TOTAL THREAT ON POMPS SCALE, max obs is the new max

hist(DF_POMPS$pTHR_Y)

############################# DEPRIVATION #####################################


########## INCOME TO POVERTY RATIO 

# I use the 2018 guidelines, with the poverty line set for 200% the poverty line. 

# convert the annual income in 1000s to the full value 
DF_POMPS$INCOME_FULL <- (DF_POMPS$INCOME_W1*1000)

# FUNCTION TO INDEX INCOME TO NEEDS BASED ON FEDERAL GUIDELINES 

# Function to compute income-to-needs ratio
INCOME_TO_NEEDS <- function(data, income_var, household_var) {
  # Define the poverty threshold index
  poverty_thresholds <- c(24280, 32920, 41560, 50200, 58840, 67480, 76120, 84760) #HOUSE HOLD NUMBER 1 - 8 PERSONS
  
  # Ensure household size does not exceed the index range
  data$poverty_threshold <- ifelse(is.na(data[[household_var]]), NA, 
                                   ifelse(data[[household_var]] >= 8, 
                                          poverty_thresholds[8], 
                                          poverty_thresholds[pmin(data[[household_var]], length(poverty_thresholds))]))
  
  # Compute the income-to-needs ratio, keeping NA values
  data$INC_TO_NEED <- ifelse(is.na(data[[income_var]]) | is.na(data$poverty_threshold), 
                                 NA, 
                                 data[[income_var]] / data$poverty_threshold)
  
  # Return updated data frame
  return(data)
}


DF_POMPS <- INCOME_TO_NEEDS(DF_POMPS, "INCOME_FULL", "HOUSE_NUMW1")

# REVERSE SCORE SO HIGH VALUES = POVERTY

DF_POMPS$INC_TO_NEED_R <- DF_POMPS$INC_TO_NEED*(-1) + max(DF_POMPS$INC_TO_NEED, na.rm = TRUE)

max(DF_POMPS$INC_TO_NEED_R, na.rm = TRUE) # MIN = 0, MAX = 2.34

########## Parent Education 

#Reverse score to make it negative
DF_POMPS$PAR_EDU_R <- DF_POMPS$PAR_EDU*(-1) + max(DF_POMPS$PAR_EDU, na.rm = TRUE)

hist(DF_POMPS$PAR_EDU_R) # MIN = 0, MAX = 5 (original scale is up to 1 - 6, now it is 0-5!)

########## MARITAL STATUS 

# RECODE TO DUMMY, MARRIED OR NOT 

DF_POMPS$NOT_MAR <- ifelse(DF_POMPS$MAR_STAT == 1, 0, 1)


########## LACKING QUALITY ATTACHMENT


DF_POMPS$POOR_ATT <- apply(df[, c("ECR3CW1", "ECR4CW1", "ECR6CW1", "ECR7CW1")], 1, function(x) {
  # Count the number of missing values
  missing_count <- sum(is.na(x))
  
  # If there are 2 or more missing values, set total to NA
  if (missing_count >= 2) {
    return(NA)
  } else {
    # Otherwise, sum the non-missing values
    return(sum(x, na.rm = TRUE))
  }
})

min(DF_POMPS$POOR_ATT, na.rm = TRUE) #3

########## PHYSICAL NEGLECT

## DOES NOT REQUIRE CHANGE AT THIS LEVEL 
names(DF_POMPS)
min(DF_POMPS$CTQCW1_PN, na.rm = TRUE)

########## EMOTIONAL NEGLECT

## DOES NOT REQUIRE CHANGE AT THIS LEVEL 

min(DF_POMPS$CTQCW1_EN, na.rm = TRUE)

########## FOOD INSECURITY

## DOES NOT REQUIRE CHANGE AT THIS LEVEL 
names(DF_POMPS)
max(DF_POMPS$FOOD_INS, na.rm = TRUE) #MIN 0 MAX 6

########## INNATENTIVE PARENTING
DF_POMPS <- RECODE_CTS(DF_POMPS, "CTS27PW1")
DF_POMPS <- RECODE_CTS(DF_POMPS, "CTS28PW1")

DF_POMPS$INN_PAR <- apply(df[, c("CTS27PW1", "CTS28PW1")], 1, function(x) {
  # Count the number of missing values
  missing_count <- sum(is.na(x))
  
  # If there are 1 or more missing values, set total to NA
  if (missing_count >= 1) {
    return(NA)
  } else {
    # Otherwise, sum the non-missing values
    return(sum(x, na.rm = TRUE))
  }
})

min(DF_POMPS$INN_PAR, na.rm = TRUE) #0


######## CREATE POMPS SCORES 

DF_POMPS_DEP <- SINGLE_POMP_SET(DF_POMPS, "INC_TO_NEED_R", min_val = 0, max_val = 2.34) # REVERSE SCORED INCOME TO NEEDS NO MAX, USED OBSEVED MAX. 
DF_POMPS_DEP <- SINGLE_POMP_SET(DF_POMPS_DEP, "PAR_EDU_R", min_val = 0, max_val = 5) #REVERSE PARENT EDUCATION (OG range is 0 -6. After recoding [and no 1 - 2] the possible range is 0 -5). 
DF_POMPS_DEP <- SINGLE_POMP_SET(DF_POMPS_DEP, "NOT_MAR", min_val = 0, max_val = 1) #MARITAL STATUS (0 = married, 1 = mot married for any reason) 
DF_POMPS_DEP <- SINGLE_POMP_SET(DF_POMPS_DEP, "POOR_ATT", min_val = 3, max_val = 28) #POOR ATTACHMENT (ECR, 4 items, 1 - 28 range)
DF_POMPS_DEP <- SINGLE_POMP_SET(DF_POMPS_DEP, "CTQCW1_PN", min_val = 5, max_val = 25) #PHYSICAL NEGLECT - CTQ (5 items, each 1 - 5)
DF_POMPS_DEP <- SINGLE_POMP_SET(DF_POMPS_DEP, "CTQCW1_EN", min_val = 5, max_val = 25) #EMOTONAL NEGLECT - CTQ (5 items, each 1 - 5)
DF_POMPS_DEP <- SINGLE_POMP_SET(DF_POMPS_DEP, "FOOD_INS", min_val = 0, max_val = 6) #FOOD INSECURITY ( 6 items, re-coded to 0 - 1)
DF_POMPS_DEP <- SINGLE_POMP_SET(DF_POMPS_DEP, "INN_PAR", min_val = 0, max_val = 14) #INNATENTIVE PARENTING - CTS ( 2 items 0 - 7)

########CREATE A TOTAL DEPRIVATION SCORE 

DF_POMPS_DEP$DEP_Y <- apply(DF_POMPS_DEP[, c("pINC_TO_NEED_R", "pPAR_EDU_R", "pNOT_MAR", "pPOOR_ATT", 
                                             "pCTQCW1_PN", "pCTQCW1_EN", "pFOOD_INS", "pINN_PAR")], 1, function(x) {
                                          # Count the number of missing values
                                          missing_count <- sum(is.na(x))
                                          
                                          # If more than 4 values are missing, set total to NA
                                          if (missing_count > 4) {
                                            return(NA)
                                          } else {
                                            # Otherwise, sum the non-missing values
                                            return(sum(x, na.rm = TRUE))
                                          }
                                        })


min(DF_POMPS_DEP$DEP_Y, na.rm = TRUE) #0
max(DF_POMPS_DEP$DEP_Y, na.rm = TRUE) #387.45

## MIN IS ZERO (SINCE THEY ARE POMPS SCORES, A ZERO ON ALL IS POSSIBLE, JUST UNLIKELY)
## MAX IS 800 (SINCE POMPS SCORES, RANGE IS 0 - 100. 8 X 100 = 800)

DF_POMPS_DEP <- SINGLE_POMP_SET(DF_POMPS_DEP, "DEP_Y", min_val = 0, max_val = 800) 


hist(DF_POMPS_DEP$pDEP_Y)

#############################  UNPREDICTABILITY ############################# 

########## NEIGHBORHOOD DISORDER

DF_POMPS_DEP$NEI_DIS <- apply(DF_POMPS_DEP[, c("NER1PW1", "NER2PW1", "NER3PW1", "NER4PW1", "NER10PW1", "NER11PW1")], 1, function(x) {
                           # Count the number of missing values
                           missing_count <- sum(is.na(x))
                           
                           # If more than 3 values are missing, set total to NA
                           if (missing_count > 3) {
                             return(NA)
                           } else {
                             # Otherwise, sum the non-missing values
                             return(sum(x, na.rm = TRUE))
                           }
                         })


min(DF_POMPS_DEP$NEI_DIS, na.rm = TRUE) #5
max(DF_POMPS_DEP$NEI_DIS, na.rm = TRUE) #58


########## FAMILY CHAOS 

## NO CHANGES AT THIS LEVEL 

min(DF_POMPS_DEP$FAC1_CH, na.rm = TRUE) #7
max(DF_POMPS_DEP$FAC1_CH, na.rm = TRUE) #27


########## PARENT ALCOHOL ABUCE 

## NO CHANGES AT THIS LEVEL 

min(DF_POMPS_DEP$ALCPW1_PB, na.rm = TRUE) #0
max(DF_POMPS_DEP$ALCPW1_PB, na.rm = TRUE) #14

########## PARENT DRUG USE  

## NO CHANGES AT THIS LEVEL 

min(DF_POMPS_DEP$DRGPW1_T, na.rm = TRUE) #0
max(DF_POMPS_DEP$DRGPW1_T, na.rm = TRUE) #12


######## CREATE POMPS SCORES 

DF_POMPS_UNP <- SINGLE_POMP_SET(DF_POMPS_DEP, "NEI_DIS", min_val = 5, max_val = 60) # ITEMS FROM NEIGHBORHOOD ENVIONMENT PARCEL (six items, scale of 1-10 each)
DF_POMPS_UNP <- SINGLE_POMP_SET(DF_POMPS_UNP, "FAC1_CH", min_val = 7, max_val = 35) #FAMILY CHAOS RAW SCALE - FACES (7 items, 1 - 5 scale)
DF_POMPS_UNP <- SINGLE_POMP_SET(DF_POMPS_UNP, "ALCPW1_PB", min_val = 0, max_val = 28) #ALCOHOL PROBLEMS SCALE (TOTAL 0 - 28 point scale) 
DF_POMPS_UNP <- SINGLE_POMP_SET(DF_POMPS_UNP, "DRGPW1_T", min_val = 0, max_val = 44) #DRUG USE FREQUENCY & PROBLEMS (DUDIT, 0 - 44 total scale)

########CREATE A TOTAL UNPREDICTABILITY SCORE 

DF_POMPS_UNP$UNP_Y <- apply(DF_POMPS_UNP[, c("pNEI_DIS", "pFAC1_CH", "pALCPW1_PB", "pDRGPW1_T")], 1, function(x) {
                                               # Count the number of missing values
                                               missing_count <- sum(is.na(x))
                                               
                                               # If more than 2 values are missing, set total to NA
                                               if (missing_count > 2) {
                                                 return(NA)
                                               } else {
                                                 # Otherwise, sum the non-missing values
                                                 return(sum(x, na.rm = TRUE))
                                               }
                                             })


min(DF_POMPS_UNP$UNP_Y, na.rm = TRUE) #0
max(DF_POMPS_UNP$UNP_Y, na.rm = TRUE) #148.44
## MIN IS ZERO (SINCE THEY ARE POMPS SCORES, A ZERO ON ALL IS POSSIBLE, JUST UNLIKELY)
## MAX IS 800 (SINCE POMPS SCORES, RANGE IS 0 - 100. 4 X 100 = 400)

DF_POMPS_UNP <- SINGLE_POMP_SET(DF_POMPS_UNP, "UNP_Y", min_val = 0, max_val = 400) #TOTAL DEP ON POMPS SCALE, max obs is the new max


hist(DF_POMPS_UNP$pUNP_Y)

#############################  CUMULATIVE ADV #############################


DF_POMPS_UNP$CUM_Y <- apply(DF_POMPS_UNP[, c("pCRIM_ACT", "pWRECK_ACT", "pUNSAFE_NEI", "pPHY_DISP", 
                                             "pPSY_AGG", "pCTQCW1_PA", "pCTQCW1_SA", "pCTQCW1_EA",
                                             "pINC_TO_NEED_R", "pPAR_EDU_R", "pNOT_MAR", "pPOOR_ATT", 
                                             "pCTQCW1_PN", "pCTQCW1_EN", "pFOOD_INS", "pINN_PAR",
                                             "pNEI_DIS", "pFAC1_CH", "pALCPW1_PB", "pDRGPW1_T")], 1, function(x) {
                                               # Count the number of missing values
                                               missing_count <- sum(is.na(x))
                                               
                                               # If more than 10 values are missing, set total to NA
                                               if (missing_count > 10) {
                                                 return(NA)
                                               } else {
                                                 # Otherwise, sum the non-missing values
                                                 return(sum(x, na.rm = TRUE))
                                               }
                                             })

min(DF_POMPS_UNP$CUM_Y, na.rm = TRUE) #0
max(DF_POMPS_UNP$CUM_Y, na.rm = TRUE) #581.59

## MIN IS ZERO (SINCE THEY ARE POMPS SCORES, A ZERO ON ALL IS POSSIBLE, JUST UNLIKELY)
## MAX IS 2000 (SINCE POMPS SCORES, RANGE IS 0 - 100. 20 X 100 = 2000)

DF_POMPS_UNP <- SINGLE_POMP_SET(DF_POMPS_UNP, "CUM_Y", min_val = 0, max_val = 2000) #TOTAL DEP ON POMPS SCALE, max obs is the new max


hist(DF_POMPS_UNP$pCUM_Y)


################ PARENT ELA

##############################  THREAT ##################################

DF_POMPS_UNP$THR_P <- apply(DF_POMPS_UNP[, c("P_CTQ_EA", "P_CTQ_PA", "P_CTQ_SA")], 1, function(x) {
  # Count the number of missing values
  missing_count <- sum(is.na(x))
  
  # If more than 1 values are missing, set total to NA
  if (missing_count > 1) {
    return(NA)
  } else {
    # Otherwise, sum the non-missing values
    return(sum(x, na.rm = TRUE))
  }
})

DF_POMPS_UNP <- SINGLE_POMP_SET(DF_POMPS_UNP, "THR_P", min_val = 15, max_val = 75) # TOTAL PARENT THREAT - CTQ, 3 scale of 5 - 25 each (min 15, max 75)


#############################  DEPRIVATION #############################

DF_POMPS_UNP$DEP_P <- apply(DF_POMPS_UNP[, c("P_CTQ_PN", "P_CTQ_EN", "P_CTQ_MD")], 1, function(x) {
  # Count the number of missing values
  missing_count <- sum(is.na(x))
  
  # If more than 1 values are missing, set total to NA
  if (missing_count > 1) {
    return(NA)
  } else {
    # Otherwise, sum the non-missing values
    return(sum(x, na.rm = TRUE))
  }
})

DF_POMPS_UNP <- SINGLE_POMP_SET(DF_POMPS_UNP, "DEP_P", min_val = 10, max_val = 53) # TOTAL PARENT THREAT - CTQ, 3 scale of 2 are 5 - 25, 1 is 0 - 3 


#############################  CUMULATIVE  #############################

DF_POMPS_UNP$CUM_P <- apply(DF_POMPS_UNP[, c("P_CTQ_PN", "P_CTQ_EN", "P_CTQ_MD",
                                             "P_CTQ_EA", "P_CTQ_PA", "P_CTQ_SA")], 1, function(x) {
  # Count the number of missing values
  missing_count <- sum(is.na(x))
  
  # If more than 1 values are missing, set total to NA
  if (missing_count > 3) {
    return(NA)
  } else {
    # Otherwise, sum the non-missing values
    return(sum(x, na.rm = TRUE))
  }
})

DF_POMPS_UNP <- SINGLE_POMP_SET(DF_POMPS_UNP, "CUM_P", min_val = 25, max_val = 128) # TOTAL PARENT THREAT - CTQ, 3 scale of 2 are 5 - 25, 1 is 0 - 3 

### VARIABLE PREP ON P-ELA and Y-ELA

## CENTER VARIABLES 

###### Data centering 

CENTERFUCKS <- function(x) {
  mean_x <- mean(x, na.rm = TRUE)
  x <- x - mean_x
  return(x)
}

#Applied centering 

# INDEX VARIABLES TO CENTER 
ADV_CEN_VARS <- c("pTHR_Y", "pUNP_Y", "pDEP_Y", "pCUM_Y", "pTHR_P", "pDEP_P", "pCUM_P")  

# CENTER DEM SHITS 
DF_POMPS_UNP_CENT <- DF_POMPS_UNP %>%
  mutate(across(all_of(ADV_CEN_VARS), CENTERFUCKS, .names = "C_{.col}"))

## QUADRTIFY VARIBLES 

QUADFUCKS <- function(x) {
  return(x^2)
}

# INDEX VARIABLES TO QUAD 
ADV_QUAD_VARS <- c("C_pTHR_Y", "C_pUNP_Y", "C_pDEP_Y", "C_pCUM_Y", "C_pTHR_P", "C_pDEP_P", "C_pCUM_P")  

# QUAD DEM CENTERED SHITS 
DF_POMPS_UNP_CENT_QUAD <- DF_POMPS_UNP_CENT %>%
  mutate(across(all_of(ADV_QUAD_VARS), QUADFUCKS, .names = "Q_{.col}"))

######################################################################
################## SAVE THE ADVERSITY INCLUSIVE DF 


write.csv(DF_POMPS_UNP_CENT_QUAD, "DISSERTATION_ADVERSITY_DORRY_7.6.25.csv", row.names=FALSE, na="")




