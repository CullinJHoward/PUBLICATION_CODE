#Library
library(dplyr)
library(MplusAutomation)

# Set working directory 


POUNDTOWN <- 1

# set working directory
if (POUNDTOWN == 1) {
  work_dir <- 'C:\\Users\\cjh37695\\OneDrive - University of Georgia\\Suicide sleep paper\\ANALYSIS\\ANALYSIS_R&R\\'
} else {
  work_dir <- 'C:\\Users\\0910h\\OneDrive - University of Georgia\\Suicide sleep paper\\ANALYSIS\\ANALYSIS_R&R\\MEASUREMENT_MODEL\\'
}

setwd(work_dir)



###############################################################################
############# MERGE IN SAVED LATENT VARIABLES #################################
###############################################################################

## ADVERSITY VARIABLES 

#### Load exported factor TXT file (from Mplus)


HSES_DF <- read.table("HARSH_SES_SAVE_LAT.txt", header = FALSE)

# Subset the data frame to keep only NUMID, and the latent scores you want
names(HSES_DF)

HSES_FACTOR <- HSES_DF[, c("V7", "V8", "V10")]

# Give names to the variables
colnames(HSES_FACTOR) <- c("HSES", "HSES_SE", "NUMID")

names(HSES_FACTOR)

###### Merge back into the MasterDF

# Load the Master data frame 

df<- read.csv("C:\\Users\\0910h\\OneDrive - University of Georgia\\Suicide sleep paper\\ANALYSIS\\ANALYSIS_R&R\\DATA_WRANGLE\\ABCD_SLP_DMN_SUI_RR_6.11.25.csv")

# DROP OLD HSES THAT DIDN"T HAVE UNEMPLOYMENT
df_RED <- df %>%
  select(-c("H_SES", "C_H_SES", "Q_C_H_SES"))

FULL_MASTER <- df_RED %>%
  left_join(HSES_FACTOR, by = "NUMID")

################################################################################
############################# CENTER SHIT  

setwd("C:\\Users\\0910h\\OneDrive - University of Georgia\\Suicide sleep paper\\ANALYSIS\\ANALYSIS_R&R\\")


names(FULL_MASTER)
FULL_MASTER$C_HSES <- FULL_MASTER$HSES - mean(FULL_MASTER$HSES, na.rm = TRUE)
FULL_MASTER$C_CAxDpRP_1 <- FULL_MASTER$CAxDpRP_1 - mean(FULL_MASTER$CAxDpRP_1, na.rm = TRUE)
FULL_MASTER$C_CAxDpTP_1 <- FULL_MASTER$CAxDpTP_1 - mean(FULL_MASTER$CAxDpTP_1, na.rm = TRUE)
FULL_MASTER$C_CIntSRP_1 <- FULL_MASTER$CIntSRP_1 - mean(FULL_MASTER$CIntSRP_1, na.rm = TRUE)
FULL_MASTER$C_CIntSTP_1 <- FULL_MASTER$CIntSTP_1 - mean(FULL_MASTER$CIntSTP_1, na.rm = TRUE)
FULL_MASTER$C_CWtDpRP_1<- FULL_MASTER$CWtDpRP_1 - mean(FULL_MASTER$CWtDpRP_1, na.rm = TRUE)
FULL_MASTER$C_CWtDpTP_1 <- FULL_MASTER$CWtDpTP_1 - mean(FULL_MASTER$CWtDpTP_1, na.rm = TRUE)
FULL_MASTER$C_CTotPTP_1 <- FULL_MASTER$CTotPTP_1 - mean(FULL_MASTER$CTotPTP_1, na.rm = TRUE)
FULL_MASTER$C_CTotPRP_1 <- FULL_MASTER$CTotPRP_1 - mean(FULL_MASTER$CTotPRP_1, na.rm = TRUE)

################################################################################
############################# ID CLINICAL CUTOFFS

FULL_MASTER$CLI_INT <- ifelse(FULL_MASTER$CIntSTP_1 >= 63, 1, 0)

################################################################################
############################# POMPS SCORES 


## POMPS - BUT RETAIN THE MEANING OF ZERO

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
    ((data[[var_name]] - min_val) / (max_val - min_val))
  )
  
  # Compute the POMPS value at the specified center_value
  poms_at_center <- ((center_value - min_val) / (max_val - min_val))
  
  # Adjust the POMPS scores by centering on the specified value
  adjusted_poms_score <- poms_score - poms_at_center
  
  # Add the adjusted POMPS variable to the dataset with a prefix "p"
  new_var_name <- paste0("C_p", var_name)
  data[[new_var_name]] <- adjusted_poms_score
  
  # Return the modified dataset
  return(data)
}


## POMPS BUT NOT CENTERED 

# Modified RECENTERED_POMPS function to calculate relative POMP scores without centering
INVARIANT_POMPS <- function(data, var_name) {
  # Check if the variable exists in the data
  if (!(var_name %in% names(data))) {
    stop("Variable not found in the dataset.")
  }
  
  # Convert the variable to numeric if it's not already
  data[[var_name]] <- as.numeric(data[[var_name]])
  
  # Calculate the maximum and minimum values for the variable
  max_val <- max(data[[var_name]], na.rm = TRUE)
  min_val <- min(data[[var_name]], na.rm = TRUE)
  
  # Calculate the relative POMP score for each observation (0 to 100)
  poms_score <- ifelse(
    is.na(data[[var_name]]), 
    NA, 
    ((data[[var_name]] - min_val) / (max_val - min_val)) * 100
  )
  
  # Add the relative POMP score variable to the dataset with a prefix "p"
  new_var_name <- paste0("p", var_name)
  data[[new_var_name]] <- poms_score
  
  # Return the modified dataset
  return(data)
}


## APPLIED RECENTERED POMPS 

names(FULL_MASTER)

FULL_MASTER_POMP <- RECENTERED_POMPS(FULL_MASTER, var_name = "C_HSES", center_value = 0) # PARENT P-FACTOR
FULL_MASTER_POMP <- RECENTERED_POMPS(FULL_MASTER_POMP, var_name = "C_dmn_dmn_W1", center_value = 0) # PARENT P-FACTOR
FULL_MASTER_POMP <- RECENTERED_POMPS(FULL_MASTER_POMP, var_name = "DMNrs_L", center_value = 0) # PARENT P-FACTOR
FULL_MASTER_POMP <- RECENTERED_POMPS(FULL_MASTER_POMP, var_name = "C_sleep_1_p_W1", center_value = 0) # PARENT P-FACTOR
FULL_MASTER_POMP <- RECENTERED_POMPS(FULL_MASTER_POMP, var_name = "C_sds_p_ss_dims_W1", center_value = 0) # PARENT P-FACTOR
FULL_MASTER_POMP <- RECENTERED_POMPS(FULL_MASTER_POMP, var_name = "PoSLpDrs_L", center_value = 0) # PARENT P-FACTOR
FULL_MASTER_POMP <- RECENTERED_POMPS(FULL_MASTER_POMP, var_name = "PoSLpDrs_S", center_value = 0) # PARENT P-FACTOR
FULL_MASTER_POMP <- RECENTERED_POMPS(FULL_MASTER_POMP, var_name = "SLPROrs_L", center_value = 0) # PARENT P-FACTOR
FULL_MASTER_POMP <- RECENTERED_POMPS(FULL_MASTER_POMP, var_name = "SLPROrs_S", center_value = 0) # PARENT P-FACTOR
FULL_MASTER_POMP <- RECENTERED_POMPS(FULL_MASTER_POMP, var_name = "C_Y_AGE1", center_value = 0) # PARENT P-FACTOR
FULL_MASTER_POMP <- RECENTERED_POMPS(FULL_MASTER_POMP, var_name = "C_Y_AGE3", center_value = 0) # PARENT P-FACTOR
FULL_MASTER_POMP <- RECENTERED_POMPS(FULL_MASTER_POMP, var_name = "C_Y_AGE5", center_value = 0) # PARENT P-FACTOR
FULL_MASTER_POMP <- RECENTERED_POMPS(FULL_MASTER_POMP, var_name = "C_CAxDpRP_1", center_value = 0) # PARENT P-FACTOR
FULL_MASTER_POMP <- RECENTERED_POMPS(FULL_MASTER_POMP, var_name = "C_CAxDpTP_1", center_value = 0) # PARENT P-FACTOR
FULL_MASTER_POMP <- RECENTERED_POMPS(FULL_MASTER_POMP, var_name = "C_CIntSTP_1", center_value = 0) # PARENT P-FACTOR
FULL_MASTER_POMP <- RECENTERED_POMPS(FULL_MASTER_POMP, var_name = "C_CIntSRP_1", center_value = 0) # PARENT P-FACTOR
FULL_MASTER_POMP <- RECENTERED_POMPS(FULL_MASTER_POMP, var_name = "C_CWtDpRP_1", center_value = 0) # PARENT P-FACTOR
FULL_MASTER_POMP <- RECENTERED_POMPS(FULL_MASTER_POMP, var_name = "C_CWtDpTP_1", center_value = 0) #
FULL_MASTER_POMP <- RECENTERED_POMPS(FULL_MASTER_POMP, var_name = "C_LAG_1_3", center_value = 0) # PARENT P-FACTOR
FULL_MASTER_POMP <- RECENTERED_POMPS(FULL_MASTER_POMP, var_name = "C_LAG_1_5", center_value = 0) # PARENT P-FACTOR
FULL_MASTER_POMP <- RECENTERED_POMPS(FULL_MASTER_POMP, var_name = "C_LAG_3_5", center_value = 0) #
FULL_MASTER_POMP <- RECENTERED_POMPS(FULL_MASTER_POMP, var_name = "C_mean_motion_W1", center_value = 0) #
FULL_MASTER_POMP <- RECENTERED_POMPS(FULL_MASTER_POMP, var_name = "C_CTotPTP_1", center_value = 0) #
FULL_MASTER_POMP <- RECENTERED_POMPS(FULL_MASTER_POMP, var_name = "C_CTotPRP_1", center_value = 0) #


################################################################################
############################# SAVE IT! 

## ADD SCANNER 

write.csv(FULL_MASTER_POMP, "ABCD_SLP_DMN_SUI_RR_6.25.25.csv", row.names=FALSE, na="")

prepareMplusData(FULL_MASTER_POMP,"ABCD_SLP_DMN_SUI_RR_6.25.25.dat")


################################################################################
############################# ADD IN OTHER SLEEP VARIABLES 

df <- read.csv("ABCD_SLP_DMN_SUI_RR_6.25.25.csv")

SLP_DF <- read.csv("DATA_WRANGLE\\sleepdata.csv")

SLP_DF_RED <- SLP_DF %>%
  select(c(src_subject_id, 
            sleepdisturb1_p_t1, sleepdisturb1_p_t3,
            sleepdisturb2_p_t1, sleepdisturb2_p_t3,
            sleepdisturb3_p_t1, sleepdisturb3_p_t3,
            sleepdisturb4_p_t1, sleepdisturb4_p_t3,
            sleepdisturb5_p_t1, sleepdisturb5_p_t3,
            sleepdisturb10_p_t1, sleepdisturb10_p_t3,
            sleepdisturb11_p_t1, sleepdisturb11_p_t3))

SLP_DF_RED <- SLP_DF_RED %>%
  rename(subID = src_subject_id,
         SLDUR1 = sleepdisturb1_p_t1, # Sleep duration 
         SLDUR3 = sleepdisturb1_p_t3,
         SLLAT1 = sleepdisturb2_p_t1, #Sleep Latency 
         SLLAT3 = sleepdisturb2_p_t3, 
         BEDREL1 = sleepdisturb3_p_t1, # Reluctance To go to bed 
         BEDREL3 = sleepdisturb3_p_t3,
         DIFFSLP1 = sleepdisturb4_p_t1, # Difficutly falling asleep 
         DIFFSLP3 = sleepdisturb4_p_t3,
         ANXSLP1 = sleepdisturb5_p_t1, #Amnxiety when falling asleep 
         ANXSLP3 = sleepdisturb5_p_t3,
         NGTWAK1 = sleepdisturb10_p_t1, # Night wakings
         NGTWAK3 = sleepdisturb10_p_t3,
         DFSLaWA1 = sleepdisturb11_p_t1, # Difficulty falling asleep after waking 
         DFSLaWA3 = sleepdisturb11_p_t3)

names(SLP_DF_RED)

################################################################################
################ Create residualized change scores 

SLEEP_RED <- SLP_DF_RED %>%
  filter(complete.cases(SLDUR1,   SLDUR3,   SLLAT1,   SLLAT3,   BEDREL1,
                        BEDREL3,  DIFFSLP1, DIFFSLP3, ANXSLP1, 
                        ANXSLP3,  NGTWAK1,  NGTWAK3,  DFSLaWA1, DFSLaWA3))


# Regression models 

SLPDRmodel <- lm(SLDUR3 ~  SLDUR1, data = SLEEP_RED)
summary(SLPDRmodel) 

SLLATmodel <- lm(SLLAT3 ~  SLLAT1, data = SLEEP_RED)
summary(SLLATmodel) 

BEDRLmodel <- lm(BEDREL3 ~  BEDREL1, data = SLEEP_RED)
summary(BEDRLmodel) 

DIFFSLPmodel <- lm(DIFFSLP3 ~  DIFFSLP1, data = SLEEP_RED)
summary(DIFFSLPmodel) 

ANXSLPmodel <- lm(ANXSLP3 ~  ANXSLP1, data = SLEEP_RED)
summary(ANXSLPmodel) 

NGTWAKmodel <- lm(NGTWAK3 ~  NGTWAK1, data = SLEEP_RED)
summary(NGTWAKmodel) 

DFSLmodel <- lm(DFSLaWA3 ~  DFSLaWA1, data = SLEEP_RED)
summary(DFSLmodel) 

# Extract the residuals
SLEEP_RED$SLDURrs <- resid(SLPDRmodel)
SLEEP_RED$SLLATrs <- resid(SLLATmodel)
SLEEP_RED$BEDRELrs <- resid(BEDRLmodel)
SLEEP_RED$DIFFSLrs <- resid(DIFFSLPmodel)
SLEEP_RED$ANXSLPrs <- resid(ANXSLPmodel)
SLEEP_RED$NGTWAKrs <- resid(NGTWAKmodel)
SLEEP_RED$DFSLaWArs <- resid(DFSLmodel)

# Subset just the unique residual variables I just made to simplify the final merge

SLEEP_RED_RES <- SLEEP_RED %>%
  select(subID, SLDURrs, SLLATrs, BEDRELrs, DIFFSLrs, ANXSLPrs, NGTWAKrs, DFSLaWArs)

# Merge the residuals dataframe back into the original dataframe
FULL_SLP_VARS <- SLP_DF_RED %>%
  left_join(SLEEP_RED_RES, by = "subID")

names(FULL_SLP_VARS)


## Mean Center Sleep Variables 

FULL_SLP_VARS$C_SLDUR1 <- (FULL_SLP_VARS$SLDUR1 - mean(FULL_SLP_VARS$SLDUR1, na.rm = TRUE))
FULL_SLP_VARS$C_SLLAT1 <- (FULL_SLP_VARS$SLLAT1 - mean(FULL_SLP_VARS$SLLAT1, na.rm = TRUE))
FULL_SLP_VARS$C_BEDREL1 <- (FULL_SLP_VARS$BEDREL1 - mean(FULL_SLP_VARS$BEDREL1, na.rm = TRUE))
FULL_SLP_VARS$C_DIFFSLP1 <- (FULL_SLP_VARS$DIFFSLP1 - mean(FULL_SLP_VARS$DIFFSLP1, na.rm = TRUE))
FULL_SLP_VARS$C_ANXSLP1 <- (FULL_SLP_VARS$ANXSLP1 - mean(FULL_SLP_VARS$ANXSLP1, na.rm = TRUE))
FULL_SLP_VARS$C_NGTWAK1 <- (FULL_SLP_VARS$NGTWAK1 - mean(FULL_SLP_VARS$NGTWAK1, na.rm = TRUE))
FULL_SLP_VARS$C_DFSLaWA1 <- (FULL_SLP_VARS$DFSLaWA1 - mean(FULL_SLP_VARS$DFSLaWA1, na.rm = TRUE))


## POMPS TRANSFORM NEW VARIABLES 
FULL_SLP_VARS_POMP <- RECENTERED_POMPS(FULL_SLP_VARS, var_name = "C_SLDUR1", center_value = 0)
FULL_SLP_VARS_POMP <- RECENTERED_POMPS(FULL_SLP_VARS_POMP, var_name = "C_SLLAT1", center_value = 0) 
FULL_SLP_VARS_POMP <- RECENTERED_POMPS(FULL_SLP_VARS_POMP, var_name = "C_BEDREL1", center_value = 0)
FULL_SLP_VARS_POMP <- RECENTERED_POMPS(FULL_SLP_VARS_POMP, var_name = "C_DIFFSLP1", center_value = 0)
FULL_SLP_VARS_POMP <- RECENTERED_POMPS(FULL_SLP_VARS_POMP, var_name = "C_ANXSLP1", center_value = 0) 
FULL_SLP_VARS_POMP <- RECENTERED_POMPS(FULL_SLP_VARS_POMP, var_name = "C_NGTWAK1", center_value = 0)
FULL_SLP_VARS_POMP <- RECENTERED_POMPS(FULL_SLP_VARS_POMP, var_name = "C_DFSLaWA1", center_value = 0)

FULL_SLP_VARS_POMP <- INVARIANT_POMPS(FULL_SLP_VARS_POMP, var_name = "SLDUR3") 
FULL_SLP_VARS_POMP <- INVARIANT_POMPS(FULL_SLP_VARS_POMP, var_name = "SLLAT3") 
FULL_SLP_VARS_POMP <- INVARIANT_POMPS(FULL_SLP_VARS_POMP, var_name = "BEDREL3") 
FULL_SLP_VARS_POMP <- INVARIANT_POMPS(FULL_SLP_VARS_POMP, var_name = "DIFFSLP3") 
FULL_SLP_VARS_POMP <- INVARIANT_POMPS(FULL_SLP_VARS_POMP, var_name = "ANXSLP3") 
FULL_SLP_VARS_POMP <- INVARIANT_POMPS(FULL_SLP_VARS_POMP, var_name = "NGTWAK3") 
FULL_SLP_VARS_POMP <- INVARIANT_POMPS(FULL_SLP_VARS_POMP, var_name = "DFSLaWA3") 

FULL_SLP_VARS_POMP <- RECENTERED_POMPS(FULL_SLP_VARS_POMP, var_name = "SLDURrs", center_value = 0) 
FULL_SLP_VARS_POMP <- RECENTERED_POMPS(FULL_SLP_VARS_POMP, var_name = "SLLATrs", center_value = 0) 
FULL_SLP_VARS_POMP <- RECENTERED_POMPS(FULL_SLP_VARS_POMP, var_name = "BEDRELrs", center_value = 0) 
FULL_SLP_VARS_POMP <- RECENTERED_POMPS(FULL_SLP_VARS_POMP, var_name = "DIFFSLrs", center_value = 0) 
FULL_SLP_VARS_POMP <- RECENTERED_POMPS(FULL_SLP_VARS_POMP, var_name = "ANXSLPrs", center_value = 0) 
FULL_SLP_VARS_POMP <- RECENTERED_POMPS(FULL_SLP_VARS_POMP, var_name = "NGTWAKrs", center_value = 0) 
FULL_SLP_VARS_POMP <- RECENTERED_POMPS(FULL_SLP_VARS_POMP, var_name = "DFSLaWArs", center_value = 0) 

names(FULL_SLP_VARS_POMP)

################################################################################
#################### MERGE WITH THE EXISTING DATA FRAME 

FULL_DF <- df %>%
  left_join(FULL_SLP_VARS_POMP, by = "subID")


FULL_DF_RED <- FULL_DF %>%
  select(-c("sleep_1_p_W1", "sleep_1_p_W3", "sds_p_ss_dims_W1", "sds_p_ss_dims_W3",
           "C_sleep_1_p_W1", "C_sds_p_ss_dims_W1",  "DMNrs_L", "DMNrs_S",
           "PoSLpDrs_S", "SLPROrs_S", "Gen_LATF", "Soc_LATF", "Pnat_LATF",
           "CAxDpRP_5", "CAxDpRP_7", "CAxDpRP_9", "CAxDpTP_5", "CAxDpTP_7", 
           "CAxDpTP_9", "CIntSRP_5", "CIntSRY_2", "CIntSRY_3", "CIntSRY_4",
           "CIntSRY_5", "CIntSTP_5", "CIntSTP_7", "CIntSTP_9", "CORE_VD_1",
           "CORE_VD_2", "CORE_VD_3", "CORE_VD_4", "CORE_VD_5", "CSomSRP_1",
           "CSomSTP_1", "CSomSTP_3", "CSomSTP_5", "CSomSTP_7", "CSomSTP_9", 
           "CTotPRP_5", "CTotPRP_7", "CTotPRP_9", "CTotPTP_5", "CTotPTP_7",
           "CTotPTP_9", "CWtDpRP_5", "CWtDpRP_7", "CWtDpRP_9", "CWtDpTP_5",
           "CWtDpTP_7", "CWtDpTP_9", "INC_FND_5", "INC_FND_7", "INC_FND_9",
           "C_pC_sleep_1_p_W1", "C_pC_sds_p_ss_dims_W1", "C_pPoSLpDrs_S" ,
           "C_pSLPROrs_S"))

names(FULL_DF_RED)

###############################################################################
#################### SAVE IT 

## ADD SCANNER 

write.csv(FULL_DF_RED, "ABCD_SLP_DMN_SUI_RR_7.1.25.csv", row.names=FALSE, na="")

prepareMplusData(FULL_DF_RED,"ABCD_SLP_DMN_SUI_RR_7.1.25.dat")




