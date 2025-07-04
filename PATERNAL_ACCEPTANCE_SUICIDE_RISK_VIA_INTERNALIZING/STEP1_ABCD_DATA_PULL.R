#title: "ABCD Data Pull"



################################################################################
###########Load packages and create directory/variable objects 

library(dplyr)
library(tidyverse)
library(fauxnaif)
library(summarytools)
library(psych)
library(DescTools)
library(MplusAutomation)
library(purrr)



################################################################################
########### SET UP THE REQUIRED DATA LOCATIONS FOR THE PULL ####################
################################################################################


# REQUESTORS DIRECTORY 

setwd("C:\\Users\\0910h\\OneDrive - University of Georgia\\YDI-Shared\\DATA_MANAGEMENT\\DATA_REQUESTS\\CULLIN\\ABCD_FATHER_ACCxNEURO_SUIC\\")

# LOAD DATA REQUEST FORM

REQUEST <- read.csv("VAR_REQUEST.csv")

# IDENTIFY THE LOCATION OF THE ABCD DATA FILES 

ABCD_DIRECTORY <- "F:\\ABCD 5.0\\abcd 5.0\\core\\"

# LIST THE TABLES CONTAINING THE REQUESTED DATA

csv_files <- paste0(unique(REQUEST$Tables), ".csv") 

# ADD TABLES THAT HAVE DEMOGRAPHICS AND OTHER DEFAULT DATA IN THEM

csv_files <- unique(c(csv_files, "abcd_y_lt.csv", "abcd_p_demo.csv", "abcd_y_lf.csv"))

# IMAGING CONSIDERATIONS
### Check if any csv file starting with "mri_" is present
if (any(grepl("^mri_", csv_files))) {
  # Add mri_y_qc_incl.csv and mri_y_qc_motion.csv to csv_files if present
  csv_files <- c(csv_files, "mri_y_qc_incl.csv", "mri_y_qc_motion.csv", "mri_y_adm_info.csv")
}

# LIST THE REQUESTED VARIABLES 

ABCD_VARS <- unique(REQUEST$Variable.Name)

## LOAD IN G-CDS DATA DICTIONARY

NAME_DIC <- read.csv("C:\\Users\\0910h\\OneDrive - University of Georgia\\YDI-Shared\\DATA_MANAGEMENT\\DATA_REQUESTS\\G-CDS_ABCD.VARIABLE.NAME_DICTIONARY-HDFS-1TR13Q3.csv")

#### FIND MISSING VARIABLES FROM REQUEST TO DICTIONARY ####

# List variable names in df1$Variable.Name that are not in df2$OG_NAME
not_found <- setdiff(REQUEST$Variable.Name, NAME_DIC$OG_NAME)

# Print the list of variable names not found
if (length(not_found) > 0) {
  cat("The following variable names are not found in df2$OG_NAME:\n")
  print(not_found)
} else {
  cat("All variable names in df1$Variable.Name are found in df2$OG_NAME.\n")
}


################################################################################
########### Conditional additions: Check included imaging, add QC variables

#Resting state QC
if (any(grepl("^rsfmri_", ABCD_VARS))) {
  # Add imgincl_rsfmri_include and rsfmri_meanmotion to ABCD_VARS if present
  ABCD_VARS <- c(ABCD_VARS, "imgincl_rsfmri_include", "rsfmri_meanmotion", "mri_info_deviceserialnumber")
}

#dMRI (DTI or RSI)
if (any(grepl("^mri_y_dti_|^mri_y_rsi_", csv_files))) {
  # Add imgincl_dmri_include and dmri_meanmotion to ABCD_VARS if present
  ABCD_VARS <- c(ABCD_VARS, "imgincl_dmri_include", "dmri_meanmotion", "mri_info_deviceserialnumber")
}

#Functional Task - MID 
if (any(grepl("^mri_y_tfmr_mid_", csv_files))) {
  # Add imgincl_mid_include and tfmri_mid_all_meanmotion to ABCD_VARS if present
  ABCD_VARS <- c(ABCD_VARS, "imgincl_mid_include", "tfmri_mid_all_meanmotion", "mri_info_deviceserialnumber")
}

#Functional Task - NBACK 
if (any(grepl("^mri_y_tfmr_nback_", csv_files))) {
  # Add imgincl_nback_include and tfmri_sst_all_meanmotion to ABCD_VARS if present
  ABCD_VARS <- c(ABCD_VARS, "imgincl_nback_include", "tfmri_nback_all_meanmotion", "mri_info_deviceserialnumber")
}

#Functional Task - SST 
if (any(grepl("^mri_y_tfmr_sst_", csv_files))) {
  # Add imgincl_sst_include and tfmri_nback_all_meanmotion to ABCD_VARS if present
  ABCD_VARS <- c(ABCD_VARS, "imgincl_sst_include", "tfmri_sst_all_meanmotion", "mri_info_deviceserialnumber")
}

#Structural - T1 weighted

if (any(grepl("^mri_y_smr_t1_", csv_files))) {
  # Add imgincl_t1w_include 
  ABCD_VARS <- c(ABCD_VARS, "imgincl_t1w_include", "mri_info_deviceserialnumber")
}

#Structural - T2 weighted

if (any(grepl("^mri_y_smr_t2_", csv_files))) {
  # Add imgincl_t2w_include 
  ABCD_VARS <- c(ABCD_VARS, "imgincl_t2w_include", "mri_info_deviceserialnumber")
}



################################################################################
######### RECURSIVELY SEARCH THE FILE STRUCTURE FOR NEEDED TABLES ##############
################################################################################


## CREATE SEARCHING FUNCTION 

find_csv <- function(directory, file_list) {
  # Get list of files and directories in current directory
  files <- list.files(path = directory, full.names = TRUE)
  
  # Filter out directories
  files <- files[!file.info(files)$isdir]
  
  # Iterate over files
  for (file in files) {
    # Check if the file is in the list of CSV files
    if (basename(file) %in% file_list) {
      # Load CSV file into a data frame
      df_name <- tools::file_path_sans_ext(basename(file))
      assign(df_name, read.csv(file), envir = .GlobalEnv)
      cat("Loaded", df_name, "from", file, "\n")
      # Remove the found file from the list
      file_list <- file_list[file_list != basename(file)]
    }
  }
  
  # Get list of directories in current directory
  directories <- list.dirs(path = directory, full.names = TRUE, recursive = FALSE)
  
  # Recursively search subdirectories
  for (subdir in directories) {
    find_csv(subdir, file_list)
  }
}


# APPLY THE SEARCHING FUNCTION 

find_csv(ABCD_DIRECTORY, csv_files)



## MERGE DATA CSVs TOGETHER 


# List all data frames in the environment
DFs <- Filter(is.data.frame, mget(ls()))

# Exclude the "REQUEST" and "NAME_DIC" data frames if they exist
DFs <- DFs[!names(DFs) %in% c("REQUEST", "NAME_DIC")]

# MERGE THEM TO A FULL DF 
MERGED_ALLTIME_ALLVAR <- Reduce(function(x, y) 
  merge(x, y, by = c("src_subject_id", "eventname"), 
        all = TRUE, suffixes = c("", "")), DFs)


# ENSURE WE FOUND ALL REQUESTED VARIABLES

non_existing <- setdiff(ABCD_VARS, colnames(MERGED_ALLTIME_ALLVAR))

print(non_existing)


# REDUCE DF TO JUST THE REQUESTED VARIABLES AND DEFAULT DEMOS 

MERGED_ALLTIME_SPECVAR <- MERGED_ALLTIME_ALLVAR[, c("src_subject_id", "eventname", "site_id_l",
                                                    "interview_age", "demo_sex_v2", "latent_factor_ss_general_ses",
                                                    "latent_factor_ss_social", "latent_factor_ss_perinatal", "acs_raked_propensity_score",
                                                    "race_ethnicity", "demo_prnt_marital_v2", "demo_prnt_ed_v2",
                                                    "demo_prtnr_ed_v2", "demo_comb_income_v2",
                                                    "rel_family_id",
                                                    ABCD_VARS)]


# ENSURE MISSINGNESS IS THE SAME ACROSS ALL VARIABLES

MERGED_ALLTIME_SPECVAR <- MERGED_ALLTIME_SPECVAR %>%
  mutate(
    across(everything(), ~ifelse(. %in% c(999, 777, ""), NA, .)))


################################################################################
###########  Neuroimaging Quality Control 

# Resting state 

if (any(grepl("^rsfmri_", colnames(MERGED_ALLTIME_SPECVAR)))) {
  CLEAN <- MERGED_ALLTIME_SPECVAR %>%
    mutate(across(
      starts_with("rsfmri_"),  # Select columns that start with "rsfmri_"
      ~case_when(
        imgincl_rsfmri_include != 1  ~ NA_real_,  # Apply NA to the specified conditions
        TRUE ~ .)
    ))
}


# NBACK 
if (any(grepl("^tfmri_", colnames(MERGED_ALLTIME_SPECVAR)))) {
  CLEAN <- MERGED_ALLTIME_SPECVAR %>%
    mutate(across(
      starts_with("tfmri_"),  # Select columns that start with "rsfmri_"
      ~case_when(
        imgincl_nback_include != 1  ~ NA_real_,  # Apply NA to the specified conditions
        TRUE ~ .)
    ))
}

#dMRI

if (any(grepl("^rsfmri_", colnames(MERGED_ALLTIME_SPECVAR)))) {
  CLEAN <- MERGED_ALLTIME_SPECVAR %>%
    mutate(across(
      starts_with("rsfmri_"),  # Select columns that start with "rsfmri_"
      ~case_when(
        imgincl_dmri_include != 1  ~ NA_real_,  # Apply NA to the specified conditions
        TRUE ~ .)
    ))
}




################################################################################
########### RENAME VARIABLES 

###### Rename base demographic variables ######

#CLEAN_NAMED <- CLEAN %>% 
#  rename(subid=src_subject_id) %>%
#  rename(siteid=site_id_l) %>%
#  rename(yage=interview_age) %>%
#  rename(ysex=demo_sex_v2) %>% #1 = male, 2 = female
#  rename(general_lf=latent_factor_ss_general_ses) %>%
#  rename(social_lf=latent_factor_ss_social) %>%
#  rename(perinatal_lf=latent_factor_ss_perinatal) %>% 
#  rename(ppensity=acs_raked_propensity_score) %>%
#  rename(yrace=race_ethnicity) %>% #1=white, 2=black, 3=hispanic, 4=asian, 5=other
#  rename(marital=demo_prnt_marital_v2) %>% #1=married, 2=widowed, 3=divorced, 4=separated, 5=nevermarried, 6=living with partner, 777=refuse
#  rename(pedu=demo_prnt_ed_v2) %>%
#  rename(pedu2=demo_prtnr_ed_v2) %>%
#  rename(income=demo_comb_income_v2) %>%
#  rename(incomel=demo_comb_income_v2_l) %>%
#  rename(familyid=rel_family_id) %>%
#  rename(MOTION=rsfmri_meanmotion) %>%
#  rename(SCANNER = mri_info_deviceserialnumber)

#income
#1 = Less than $5,000 ; 2 = $5,000 through $11,999 ; 3 = $12,000 through $15,999 ; 
#4 = $16,000 through $24,999 ; 5 = $25,000 through $34,999 ; 6 = $35,000 through $49,999 ; 
#7 = $50,000 through $74,999 ; 8 = $75,000 through $99,999 ; 9 = $100,000 through $199,999 ;
#10 = $200,000 and greater ; 999, Don't know ; 777, Refuse to answer


## REDUCE THE SITE TO A NUMERIC CHARACTER 

MERGED_ALLTIME_SPECVAR$site_id_l <- sub("site", "", MERGED_ALLTIME_SPECVAR$site_id_l)

## CHANGE THE NAMES TO BE REDUCED AND MATCH CENTER CONVENTIONS

#CREATE DICTIONARY 
name_mapping <- setNames(NAME_DIC$G.CDS_NAME, NAME_DIC$OG_NAME)

#RENAME BASED ON MATCHING IN THE DICTIONARY 
names(MERGED_ALLTIME_SPECVAR) <- ifelse(names(MERGED_ALLTIME_SPECVAR) %in% names(name_mapping), 
                                        name_mapping[names(MERGED_ALLTIME_SPECVAR)], names(MERGED_ALLTIME_SPECVAR))


# Drop QC Variable 

#CLEAN_NAMED <- CLEAN_NAMED %>% select(-imgincl_rsfmri_include)

#names(CLEAN_NAMED)

################################################################################
########### SUBSET INTO WAVES 


#dat1 <- subset(MERGED_ALLTIME_SPECVAR, eventname == "baseline_year_1_arm_1")
dat2 <- subset(MERGED_ALLTIME_SPECVAR, eventname == "6_month_follow_up_arm_1")
dat3 <- subset(MERGED_ALLTIME_SPECVAR, eventname == "1_year_follow_up_y_arm_1")
dat4 <- subset(MERGED_ALLTIME_SPECVAR, eventname == "18_month_follow_up_arm_1")
dat5 <- subset(MERGED_ALLTIME_SPECVAR, eventname == "2_year_follow_up_y_arm_1")
dat6 <- subset(MERGED_ALLTIME_SPECVAR, eventname == "30_month_follow_up_arm_1")
dat7 <- subset(MERGED_ALLTIME_SPECVAR, eventname == "3_year_follow_up_y_arm_1")
dat8 <- subset(MERGED_ALLTIME_SPECVAR, eventname == "42_month_follow_up_arm_1")
dat9 <- subset(MERGED_ALLTIME_SPECVAR, eventname == "4_year_follow_up_y_arm_1")

# remove empty columns from each wave-------------------------------------------

#dat1 <- janitor::remove_empty(dat1, which = "cols")
dat2 <- janitor::remove_empty(dat2, which = "cols")
dat3 <- janitor::remove_empty(dat3, which = "cols")
dat4 <- janitor::remove_empty(dat4, which = "cols")
dat5 <- janitor::remove_empty(dat5, which = "cols")
dat6 <- janitor::remove_empty(dat6, which = "cols")
dat7 <- janitor::remove_empty(dat7, which = "cols")
dat8 <- janitor::remove_empty(dat8, which = "cols")
dat9 <- janitor::remove_empty(dat9, which = "cols")


## remove variables from each wave ----------------------------------------------
# Lose eventname and any other variables that are redundant or not necessary 

names(dat9)
#dat1 <- dat1 %>% 
#  select(-"eventname")
dat2 <- dat2 %>% 
  select(-c("eventname", "SiteID", "Y_AGE"))
dat3 <- dat3 %>% 
  select(-c("eventname", "SiteID",    "Y_AGE",     "ppensity",  "Y_RACE"))
dat4 <- dat4 %>% 
select(-c("eventname", "SiteID", "Y_AGE"))
dat5 <- dat5 %>% 
  select(-c("eventname", "SiteID", "Y_AGE"))
dat6 <- dat6 %>% 
  select(-c("eventname", "SiteID", "Y_AGE"))
dat7 <- dat7 %>% 
  select(-c("eventname", "SiteID", "Y_AGE"))
dat8 <- dat8 %>% 
  select(-c("eventname", "SiteID", "Y_AGE"))
dat9 <- dat9 %>% 
  select(-c("eventname", "SiteID", "Y_AGE"))

# add _wave # to end of vars per wave  -----------------------------------------

# include all except subid column, which is column #1
## colnames(dat9)[2:92] <- paste(colnames(dat9)[2:92], "1", sep = "_")


## Move variables to put the time-varying variables together 

dat1 <- dat1 %>%
  relocate(yage, .after = 15)

names(dat1)
#colnames(dat1)[2:15] <- paste(colnames(dat1)[2:15], "1", sep = "")
colnames(dat2)[2:8] <- paste(colnames(dat2)[2:8], "2", sep = "")
colnames(dat3)[2:8] <- paste(colnames(dat3)[2:8], "3", sep = "")
colnames(dat4)[2:8] <- paste(colnames(dat4)[2:8], "4", sep = "")
colnames(dat5)[2:8] <- paste(colnames(dat5)[2:8], "5", sep = "")
colnames(dat6)[2:8] <- paste(colnames(dat6)[2:8], "6", sep = "")
colnames(dat7)[2:8] <- paste(colnames(dat7)[2:8], "7", sep = "")
colnames(dat8)[2:7] <- paste(colnames(dat8)[2:7], "8", sep = "")
colnames(dat9)[2:8] <- paste(colnames(dat9)[2:8], "9", sep = "")

# merge waves ------------------------------------------------------------------


### Remove empty participants

dat2_ALL <- dat2[!apply(dat2[, -1], 1, function(x) all(is.na(x))), ]
dat3_ALL <- dat3[!apply(dat3[, -1], 1, function(x) all(is.na(x))), ]
dat4_ALL <- dat4[!apply(dat4[, -1], 1, function(x) all(is.na(x))), ]
dat5_ALL <- dat5[!apply(dat5[, -1], 1, function(x) all(is.na(x))), ]
dat6_ALL <- dat6[!apply(dat6[, -1], 1, function(x) all(is.na(x))), ]
dat7_ALL <- dat7[!apply(dat7[, -1], 1, function(x) all(is.na(x))), ]

### Remove dupliates
dat2_CLEAN <- dat2_ALL %>% distinct(subID, .keep_all = TRUE)
dat3_CLEAN <- dat3_ALL %>% distinct(subID, .keep_all = TRUE)
dat4_CLEAN <- dat4_ALL %>% distinct(subID, .keep_all = TRUE)
dat5_CLEAN <- dat5_ALL %>% distinct(subID, .keep_all = TRUE)
dat6_CLEAN <- dat6_ALL %>% distinct(subID, .keep_all = TRUE)
dat7_CLEAN <- dat7_ALL %>% distinct(subID, .keep_all = TRUE)

dat_all <- full_join(dat2_CLEAN, dat3_CLEAN, by = "subID") %>%
  full_join(dat4_CLEAN, by = "subID") %>%
  full_join(dat5_CLEAN, by = "subID")%>%
  full_join(dat6_CLEAN, by = "subID")%>%
  full_join(dat7_CLEAN, by = "subID")

## Remove duplicate rows 

# Check make a new DF of unique ID observations (elimiate duplicates)

dat_all_CLEAN <- dat_all %>% distinct(subID, .keep_all = TRUE)

##### Reorder the variables to be with like variables 

# Extract characters except the last one of variable names
var_key <- sapply(names(dat_all), function(name) gsub("\\d+$", "", name))

# Get unique keys
unique_key <- unique(var_key)

# Sort the variable names based on the extracted keys
sorted_names <- unlist(sapply(unique_key, function(key) {
  grep(paste0("^", key), names(dat_all), value = TRUE)
}))

# Reorder the data frame columns
dat_all_REORDER <- dat_all[, sorted_names]
names(dat_all_REORDER)

### MERGE INTO EXISTING DATASET 

OG_FILE <- read.csv("C:\\Users\\0910h\\Dropbox\\ABCD_FATHERS\\ABCD_NeuroXInternalizing_Suicide\\ANALYSIS\\ABCD_SUICIDE_FULL_W.LAT_7.19.24.csv")

OG_FILE <- OG_FILE %>%
  rename(subID = subid)

FINAL_DF <- left_join(OG_FILE, dat_all_REORDER, by = "subID")

################# SAVE DATA 

write.csv(FINAL_DF, "ABCD_FATHERSxNEURO_SUICIDE_1.9.24.csv", row.names=FALSE, na="")

prepareMplusData(FINAL_DF,"ABCD_FATHERSxNEURO_SUICIDE_1.9.24.dat")



