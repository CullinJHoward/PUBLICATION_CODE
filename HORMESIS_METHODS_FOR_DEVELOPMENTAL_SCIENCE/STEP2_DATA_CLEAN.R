#Library

library(dplyr)
library(MplusAutomation)

# Set working directory 

setwd("/home/cjh37695/ABCD_PROJECTS/HORMESIS_METHODS/")

#Library

df <- read.csv("ABCD_HORM_METH_3.30.26.csv")


################################################################################
##################### PREP ADVERSITY VARIABLES  
names(df)

### EARLY LIFE COGNITIVE RESOURCES (HOUSEHOLD SHORT FORM)
# THESE ARE SCALED WEIRD, SOME ARE 0 - 3 AND SOME 0 - 5. 
# Coded so that bad are high

df$ErlHs1PR <- 
  dplyr::recode(df$ErlHs1P_11,
                `0` = 4,
                `1` = 3,
                `2` = 2,
                `3` = 1)
df$ErlHs2PR <- 
  dplyr::recode(df$ErlHs2P_11,
                `0` = 4,
                `1` = 3,
                `2` = 2,
                `3` = 1)
df$ErlHs3PR <- 
  dplyr::recode(df$ErlHs3P_11,
                `0` = 4,
                `1` = 3,
                `2` = 2,
                `3` = 1)
df$ErlHs4PR <- 
  dplyr::recode(df$ErlHs4P_11,
                `0` = 6,
                `1` = 5,
                `2` = 4,
                `3` = 3,
                `4` = 2,
                `5` = 1)
df$ErlHs5PR <- 
  dplyr::recode(df$ErlHs5P_11,
                `0` = 4,
                `1` = 3,
                `2` = 2,
                `3` = 1)
# ITEM 6
df$ErlHs6PR <- 
  dplyr::recode(df$ErlHs6P_11,
                `0` = 4,
                `1` = 3,
                `2` = 2,
                `3` = 1)

#ITEM 7 
df$ErlHs7PR <- (df$ErlHs7P_11*(-1)) + max(df$ErlHs7P_11, na.rm = TRUE)

#ITEM 8 
df$ErlHs8PR <- (df$ErlHs8P_11*(-1)) + max(df$ErlHs8P_11, na.rm = TRUE)

#ITEM 9
df$ErlHs9PR <- (df$ErlHs9P_11*(-1)) + max(df$ErlHs9P_11, na.rm = TRUE)

### NEIGHBORHOOD SAFETY
#W1
df$NBHsaf1P_1R <- 
  max(df$NBHsaf1P_1, na.rm = TRUE) +
  min(df$NBHsaf1P_1, na.rm = TRUE) -
  df$NBHsaf1P_1
df$NBHsaf2P_1R <- 
  max(df$NBHsaf2P_1, na.rm = TRUE) +
  min(df$NBHsaf2P_1, na.rm = TRUE) -
  df$NBHsaf2P_1
df$NBHsaf3P_1R <- 
  max(df$NBHsaf3P_1, na.rm = TRUE) +
  min(df$NBHsaf3P_1, na.rm = TRUE) -
  df$NBHsaf3P_1
#W3
df$NBHsaf1P_3R <- 
  max(df$NBHsaf1P_3, na.rm = TRUE) +
  min(df$NBHsaf1P_3, na.rm = TRUE) -
  df$NBHsaf1P_3
df$NBHsaf2P_3R <- 
  max(df$NBHsaf2P_3, na.rm = TRUE) +
  min(df$NBHsaf2P_3, na.rm = TRUE) -
  df$NBHsaf2P_3
df$NBHsaf3P_3R <- 
  max(df$NBHsaf3P_3, na.rm = TRUE) +
  min(df$NBHsaf3P_3, na.rm = TRUE) -
  df$NBHsaf3P_3
#W5
df$NBHsaf1P_5R <- 
  max(df$NBHsaf1P_5, na.rm = TRUE) +
  min(df$NBHsaf1P_5, na.rm = TRUE) -
  df$NBHsaf1P_5
df$NBHsaf2P_5R <- 
  max(df$NBHsaf2P_5, na.rm = TRUE) +
  min(df$NBHsaf2P_5, na.rm = TRUE) -
  df$NBHsaf2P_5
df$NBHsaf3P_5R <- 
  max(df$NBHsaf3P_5, na.rm = TRUE) +
  min(df$NBHsaf3P_5, na.rm = TRUE) -
  df$NBHsaf3P_5
#W7
df$NBHsaf1P_7R <- 
  max(df$NBHsaf1P_7, na.rm = TRUE) +
  min(df$NBHsaf1P_7, na.rm = TRUE) -
  df$NBHsaf1P_7
df$NBHsaf2P_7R <- 
  max(df$NBHsaf2P_7, na.rm = TRUE) +
  min(df$NBHsaf2P_7, na.rm = TRUE) -
  df$NBHsaf2P_7
df$NBHsaf3P_7R <- 
  max(df$NBHsaf3P_7, na.rm = TRUE) +
  min(df$NBHsaf3P_7, na.rm = TRUE) -
  df$NBHsaf3P_7
#W9
df$NBHsaf1P_9R <- 
  max(df$NBHsaf1P_9, na.rm = TRUE) +
  min(df$NBHsaf1P_9, na.rm = TRUE) -
  df$NBHsaf1P_9
df$NBHsaf2P_9R <- 
  max(df$NBHsaf2P_9, na.rm = TRUE) +
  min(df$NBHsaf2P_9, na.rm = TRUE) -
  df$NBHsaf2P_9
df$NBHsaf3P_9R <- 
  max(df$NBHsaf3P_9, na.rm = TRUE) +
  min(df$NBHsaf3P_9, na.rm = TRUE) -
  df$NBHsaf3P_9
#W11
df$NBHsaf1P_11R <- 
  max(df$NBHsaf1P_11, na.rm = TRUE) +
  min(df$NBHsaf1P_11, na.rm = TRUE) -
  df$NBHsaf1P_11
df$NBHsaf2P_11R <- 
  max(df$NBHsaf2P_11, na.rm = TRUE) +
  min(df$NBHsaf2P_11, na.rm = TRUE) -
  df$NBHsaf2P_11
df$NBHsaf3P_11R <- 
  max(df$NBHsaf3P_11, na.rm = TRUE) +
  min(df$NBHsaf3P_11, na.rm = TRUE) -
  df$NBHsaf3P_11

### SCHOOL SAFETY
#W1
df$SCHenv1Y_1R <- 
  max(df$SCHenv1Y_1, na.rm = TRUE) +
  min(df$SCHenv1Y_1, na.rm = TRUE) -
  df$SCHenv1Y_1
df$SCHenv2Y_1R <- 
  max(df$SCHenv2Y_1, na.rm = TRUE) +
  min(df$SCHenv2Y_1, na.rm = TRUE) -
  df$SCHenv2Y_1
df$SCHenv3Y_1R <- 
  max(df$SCHenv3Y_1, na.rm = TRUE) +
  min(df$SCHenv3Y_1, na.rm = TRUE) -
  df$SCHenv3Y_1
df$SCHenv4Y_1R <- 
  max(df$SCHenv4Y_1, na.rm = TRUE) +
  min(df$SCHenv4Y_1, na.rm = TRUE) -
  df$SCHenv4Y_1
df$SCHenv5Y_1R <- 
  max(df$SCHenv5Y_1, na.rm = TRUE) +
  min(df$SCHenv5Y_1, na.rm = TRUE) -
  df$SCHenv5Y_1
df$SCHenv6Y_1R <- 
  max(df$SCHenv6Y_1, na.rm = TRUE) +
  min(df$SCHenv6Y_1, na.rm = TRUE) -
  df$SCHenv6Y_1
#W3
df$SCHenv1Y_3R <- 
  max(df$SCHenv1Y_3, na.rm = TRUE) +
  min(df$SCHenv1Y_3, na.rm = TRUE) -
  df$SCHenv1Y_3
df$SCHenv2Y_3R <- 
  max(df$SCHenv2Y_3, na.rm = TRUE) +
  min(df$SCHenv2Y_3, na.rm = TRUE) -
  df$SCHenv2Y_3
df$SCHenv3Y_3R <- 
  max(df$SCHenv3Y_3, na.rm = TRUE) +
  min(df$SCHenv3Y_3, na.rm = TRUE) -
  df$SCHenv3Y_3
df$SCHenv4Y_3R <- 
  max(df$SCHenv4Y_3, na.rm = TRUE) +
  min(df$SCHenv4Y_3, na.rm = TRUE) -
  df$SCHenv4Y_3
df$SCHenv5Y_3R <- 
  max(df$SCHenv5Y_3, na.rm = TRUE) +
  min(df$SCHenv5Y_3, na.rm = TRUE) -
  df$SCHenv5Y_3
df$SCHenv6Y_3R <- 
  max(df$SCHenv6Y_3, na.rm = TRUE) +
  min(df$SCHenv6Y_3, na.rm = TRUE) -
  df$SCHenv6Y_3
#W5
df$SCHenv1Y_5R <- 
  max(df$SCHenv1Y_5, na.rm = TRUE) +
  min(df$SCHenv1Y_5, na.rm = TRUE) -
  df$SCHenv1Y_5
df$SCHenv2Y_5R <- 
  max(df$SCHenv2Y_5, na.rm = TRUE) +
  min(df$SCHenv2Y_5, na.rm = TRUE) -
  df$SCHenv2Y_5
df$SCHenv3Y_5R <- 
  max(df$SCHenv3Y_5, na.rm = TRUE) +
  min(df$SCHenv3Y_5, na.rm = TRUE) -
  df$SCHenv3Y_5
df$SCHenv4Y_5R <- 
  max(df$SCHenv4Y_5, na.rm = TRUE) +
  min(df$SCHenv4Y_5, na.rm = TRUE) -
  df$SCHenv4Y_5
df$SCHenv5Y_5R <- 
  max(df$SCHenv5Y_5, na.rm = TRUE) +
  min(df$SCHenv5Y_5, na.rm = TRUE) -
  df$SCHenv5Y_5
df$SCHenv6Y_5R <- 
  max(df$SCHenv6Y_5, na.rm = TRUE) +
  min(df$SCHenv6Y_5, na.rm = TRUE) -
  df$SCHenv6Y_5
#W7
df$SCHenv1Y_7R <- 
  max(df$SCHenv1Y_7, na.rm = TRUE) +
  min(df$SCHenv1Y_7, na.rm = TRUE) -
  df$SCHenv1Y_7
df$SCHenv2Y_7R <- 
  max(df$SCHenv2Y_7, na.rm = TRUE) +
  min(df$SCHenv2Y_7, na.rm = TRUE) -
  df$SCHenv2Y_7
df$SCHenv3Y_7R <- 
  max(df$SCHenv3Y_7, na.rm = TRUE) +
  min(df$SCHenv3Y_7, na.rm = TRUE) -
  df$SCHenv3Y_7
df$SCHenv4Y_7R <- 
  max(df$SCHenv4Y_7, na.rm = TRUE) +
  min(df$SCHenv4Y_7, na.rm = TRUE) -
  df$SCHenv4Y_7
df$SCHenv5Y_7R <- 
  max(df$SCHenv5Y_7, na.rm = TRUE) +
  min(df$SCHenv5Y_7, na.rm = TRUE) -
  df$SCHenv5Y_7
df$SCHenv6Y_7R <- 
  max(df$SCHenv6Y_7, na.rm = TRUE) +
  min(df$SCHenv6Y_7, na.rm = TRUE) -
  df$SCHenv6Y_7
#W9
df$SCHenv1Y_9R <- 
  max(df$SCHenv1Y_9, na.rm = TRUE) +
  min(df$SCHenv1Y_9, na.rm = TRUE) -
  df$SCHenv1Y_9
df$SCHenv2Y_9R <- 
  max(df$SCHenv2Y_9, na.rm = TRUE) +
  min(df$SCHenv2Y_9, na.rm = TRUE) -
  df$SCHenv2Y_9
df$SCHenv3Y_9R <- 
  max(df$SCHenv3Y_9, na.rm = TRUE) +
  min(df$SCHenv3Y_9, na.rm = TRUE) -
  df$SCHenv3Y_9
df$SCHenv4Y_9R <- 
  max(df$SCHenv4Y_9, na.rm = TRUE) +
  min(df$SCHenv4Y_9, na.rm = TRUE) -
  df$SCHenv4Y_9
df$SCHenv5Y_9R <- 
  max(df$SCHenv5Y_9, na.rm = TRUE) +
  min(df$SCHenv5Y_9, na.rm = TRUE) -
  df$SCHenv5Y_9
df$SCHenv6Y_9R <- 
  max(df$SCHenv6Y_9, na.rm = TRUE) +
  min(df$SCHenv6Y_9, na.rm = TRUE) -
  df$SCHenv6Y_9

##### HARSH SES 

#Reverse Score SOCIAL MOBILITY 
df$SocMob_1R <- 
  max(df$SocMob_1, na.rm = TRUE) +
  min(df$SocMob_1, na.rm = TRUE) -
  df$SocMob_1



# Reverse score parents highest education
df$HiParEdu_1R <- 
  dplyr::recode(df$HiParEdu_1,
                `1` = 5,
                `2` = 4,
                `3` = 3,
                `4` = 2,
                `5` = 1)



##### COMPUTE FAMILY CONFLICT 
# HELPER FUNCTION FOR CONFLICT SUM 

compute_scale <- function(df, items, max_na = 1) {
  na_count <- rowSums(is.na(df[items]))
  score    <- rowSums(df[items], na.rm = TRUE)
  ifelse(na_count <= max_na, score, NA_real_)
}

# DEFINE WAVES WE WANT 
waves <- c(1, 3, 5, 7, 9, 11)


# YOUTH REPORT 
for (w in waves) {
  items <- c(
    paste0("Fcon1Y_", w),
    paste0("Fcon2rY_", w),
    paste0("Fcon3Y_", w),
    paste0("Fcon4rY_", w),
    paste0("Fcon5Y_", w),
    paste0("Fcon6Y_", w),
    paste0("Fcon7rY_", w),
    paste0("Fcon8Y_", w),
    paste0("Fcon9Y_", w)
  )
  
  df[[paste0("FAMCONY_", w)]] <- compute_scale(df, items)
}


################################################################################
##################### CLEAN UP

names(df)
df_RED <- df %>%
  select(-c(
    #DROP NEIGHBORHOOD VALUES NOT AT EARLY TIMEPOINTS 
            "NBHcoh1P_5",    "NBHcoh1P_9",   
            "NBHcoh2P_5",    "NBHcoh2P_9",    "NBHcoh3P_5",    "NBHcoh3P_9",    "NBHcoh4P_5",    "NBHcoh4P_9",   
            "NBHcoh5P_5",    "NBHcoh5P_9",    "NBHisc1P_5",    "NBHisc1P_9",    "NBHisc2P_5",    "NBHisc2P_9",   
            "NBHisc3P_5",    "NBHisc3P_9",    "NBHisc4P_5",    "NBHisc4P_9",    "NBHisc5P_5",    "NBHisc5P_9",
    #DROP ORIGINAL-SCALED HOME VARIABLES 
            "ErlHs1P_11",      "ErlHs2P_11",      "ErlHs3P_11",      "ErlHs4P_11",      "ErlHs5P_11",     
            "ErlHs6P_11",      "ErlHs7P_11",      "ErlHs8P_11",      "ErlHs9P_11"))

################################################################################
##################### SAVE DATA 

# FIX FAMILY ID 0 TO BE NOT ZERO 
df_RED$FamilyID <- ifelse(df_RED$FamilyID == 0, 11884, df_RED$FamilyID)

FILENAME <- "ABCD_HORM_METH_PREP_4.23.26"

write.csv(df_RED, paste0(FILENAME, ".csv"), row.names=FALSE, na="")



