#Library

library(dplyr)
library(psych)
library(stringr)
library(purrr)
library(MplusAutomation)

# Set working directory 

setwd("/home/cjh37695/ABCD_PROJECTS/SAPPENFIELD_FATHERS/")

#Library

df <- read.csv("DATA_WRANGLE/dataset.csv", na.strings = c("999", "555", "888", "777", "-999"))
NAME_DIC <- read.csv("/home/cjh37695/ABCD_PROJECTS/G-CDS_ABCD.VARIABLE.NAME_DICTIONARY.csv")

################################################################################
#####################  RENAME THE VARIABLES 

#### FIND ANY MISSING VARIABLES FROM REQUEST TO DICTIONARY ####

not_found <- setdiff(names(df), NAME_DIC$OG_NAME_6.1)

# PRINT ANY VARIABLES NOT FOUND IN THE DATA DICTIONARY 

if (length(not_found) > 0) {
  cat("The following variable names are not found in the Center Data Dictionary\n")
  print(not_found)
} else {
  cat("All variable names in the data frame are found in the Center Data Dictionary\n")
}

## PROCEED TO RE-NAMING VARIABLES WHEN ALL ARE FOUND 

# CREATE NAMING DICTIONARY 
name_mapping <- setNames(NAME_DIC$G.CDS_NAME, NAME_DIC$OG_NAME_6.1)

#RENAME BASED ON MATCHING IN THE DICTIONARY 
names(df) <- ifelse(names(df) %in% names(name_mapping), 
                                        name_mapping[names(df)], names(df))

## CHECK HOW IT DID 
names(df)

### RECODE VISIT DATE TO BE USEFUL 

#df_NAMED <- dplyr::mutate(df_NAMED, YvDate = {
#  if ("YvDate" %in% names(df_NAMED)) {
#    date <- suppressWarnings(lubridate::ymd_hms(YvDate))
#    year <- lubridate::year(date)
#    month <- lubridate::month(date)
#    year + (month - 1) / 12
#  } else {
#    YvDate  # Return unchanged if YvDate not present
#  }
#})


################################################################################
##################### NEUROIMAING INCLUSION REQUIREMENT 

# INITIALIZE LISTS FOR TRACKING EXCLUSION
qc_na_subids <- list(
  resting_state = character(),
  nback = character(),
  mid = character(),
  dmri = character()
)

qc_na_counts <- list(
  resting_state = 0L,
  nback = 0L,
  mid = 0L,
  dmri = 0L
)

# NEW: counts per wave (eventname) for each modality
qc_na_counts_by_wave <- list(
  resting_state = integer(),
  nback = integer(),
  mid = integer(),
  dmri = integer()
)

##### Resting State functional connectivity #####
if (any(grepl("^rs", colnames(df)))) {
  
  # Rows explicitly failing QC (not included, and not missing)
  affected_rows <- df$rs_INC != 1 & !is.na(df$rs_INC)
  
  # Log excluded IDs + total count
  qc_na_subids$resting_state <- df$subID[affected_rows]
  qc_na_counts$resting_state <- sum(affected_rows)
  
  # NEW: excluded count per wave
  qc_na_counts_by_wave$resting_state <- table(df$eventname[affected_rows])
  
  # Add exclusion flag + total count column(s)
  df <- df %>%
    mutate(
      rs_excluded   = affected_rows,
      rs_excluded_n = qc_na_counts$resting_state
    ) %>%
    mutate(across(
      starts_with("rs"),
      ~case_when(
        rs_INC != 1 ~ NA_real_,
        TRUE ~ .
      )
    ))
}

# TODO: replicate the same pattern for NBACK / MID / dMRI once those blocks are implemented


################################################################################

## EXCLUSION SUMMARY 
cat("\nQC Exclusions Summary (rows made NA)\n-----------------------------------\n")

summarize_qc_wave <- function(label, total_n, by_wave_tbl) {
  cat(sprintf("%s: %d excluded\n", label, total_n))
  if (!is.null(by_wave_tbl) && length(by_wave_tbl) > 0) {
    cat("  By wave:\n")
    # print as "wave: n" lines, sorted by wave name
    by_wave_tbl <- by_wave_tbl[order(names(by_wave_tbl))]
    for (w in names(by_wave_tbl)) {
      cat(sprintf("   - %s: %d\n", w, as.integer(by_wave_tbl[[w]])))
    }
  } else {
    cat("  By wave: none\n")
  }
}

summarize_qc_wave(
  "fMRI Resting State",
  qc_na_counts$resting_state,
  qc_na_counts_by_wave$resting_state
)

summarize_qc_wave(
  "fMRI NBACK",
  qc_na_counts$nback,
  qc_na_counts_by_wave$nback
)

summarize_qc_wave(
  "fMRI MID",
  qc_na_counts$mid,
  qc_na_counts_by_wave$mid
)

summarize_qc_wave(
  "dMRI (DTI/RSI)",
  qc_na_counts$dmri,
  qc_na_counts_by_wave$dmri
)


### DROP IMAGING QC VARIABLES

qc_drop_patterns <- c(
  "excluded",   # drop anything with "excluded" anywhere in the name
  "_INC$"       # drop inclusion flags like rs_INC, nb_INC, etc. (ends with _INC)
  # "^imgincl_", # optionally re-add if you still want these gone too
  # "^[tr]fmri_.*meanmotion$"  # optionally re-add motion vars
)

# REMOVE WORKING VARIABLES 

df <- df %>%
  dplyr::select(-dplyr::matches(paste(qc_drop_patterns, collapse = "|"), ignore.case = TRUE))



################################################################################
##################### ID MOTHER AND FATHER ACCEPTANCE 

# Fill in the missing wave 1 CRPBI parent variable with the visit parent

df <- df %>%
  mutate(
    CRPBI_PARENT = if_else(
      eventname == "ses-00A",
      PRI_REL,
      CRPBI_PARENT
    )
  )

## Make father acceptance variable 

df$D_1acc <- case_when(
  is.na(df$CRPBI_PARENT) & is.na(df$SecCare) ~ NA_integer_,
  df$CRPBI_PARENT == 2 ~ df$C1acc1y,
  df$SecCare == 2 ~ df$C2acc1y,
  TRUE ~ NA_integer_
)

df$D_2acc <- case_when(
  is.na(df$CRPBI_PARENT) & is.na(df$SecCare) ~ NA_integer_,
  df$CRPBI_PARENT == 2 ~ df$C1acc2y,
  df$SecCare == 2 ~ df$C2acc2y,
  TRUE ~ NA_integer_
)

df$D_3acc <- case_when(
  is.na(df$CRPBI_PARENT) & is.na(df$SecCare) ~ NA_integer_,
  df$CRPBI_PARENT == 2 ~ df$C1acc3y,
  df$SecCare == 2 ~ df$C2acc3y,
  TRUE ~ NA_integer_
)


df$D_4acc <- case_when(
  is.na(df$CRPBI_PARENT) & is.na(df$SecCare) ~ NA_integer_,
  df$CRPBI_PARENT == 2 ~ df$C1acc4y,
  df$SecCare == 2 ~ df$C2acc4y,
  TRUE ~ NA_integer_
)

df$D_5acc <- case_when(
  is.na(df$CRPBI_PARENT) & is.na(df$SecCare) ~ NA_integer_,
  df$CRPBI_PARENT == 2 ~ df$C1acc5y,
  df$SecCare == 2 ~ df$C2acc5y,
  TRUE ~ NA_integer_
)

## Make mother acceptance variable 

df$M_1acc <- case_when(
  is.na(df$CRPBI_PARENT) & is.na(df$SecCare) ~ NA_integer_,
  df$CRPBI_PARENT == 1 ~ df$C1acc1y,
  df$SecCare == 1 ~ df$C2acc1y,
  TRUE ~ NA_integer_
)

df$M_2acc <- case_when(
  is.na(df$CRPBI_PARENT) & is.na(df$SecCare) ~ NA_integer_,
  df$CRPBI_PARENT == 1 ~ df$C1acc2y,
  df$SecCare == 1 ~ df$C2acc2y,
  TRUE ~ NA_integer_
)

df$M_3acc <- case_when(
  is.na(df$CRPBI_PARENT) & is.na(df$SecCare) ~ NA_integer_,
  df$CRPBI_PARENT == 1 ~ df$C1acc3y,
  df$SecCare == 1 ~ df$C2acc3y,
  TRUE ~ NA_integer_
)


df$M_4acc <- case_when(
  is.na(df$CRPBI_PARENT) & is.na(df$SecCare) ~ NA_integer_,
  df$CRPBI_PARENT == 1 ~ df$C1acc4y,
  df$SecCare == 1 ~ df$C2acc4y,
  TRUE ~ NA_integer_
)

df$M_5acc <- case_when(
  is.na(df$CRPBI_PARENT) & is.na(df$SecCare) ~ NA_integer_,
  df$CRPBI_PARENT == 1 ~ df$C1acc5y,
  df$SecCare == 1 ~ df$C2acc5y,
  TRUE ~ NA_integer_
)


## Summary means of mom and dad 

df$D_acc_M <- rowMeans(df[, c("D_1acc", "D_2acc", "D_3acc", "D_4acc", "D_5acc")])

df$M_acc_M <- rowMeans(df[, c("M_1acc", "M_2acc", "M_3acc", "M_4acc", "M_5acc")])

## Make a father primary caregiver variable

df$DADPRIM <- ifelse(df$PRI_REL == 2, 1, 0)

## Remove the general parental acceptance variable 

df <- df %>%
  select(-c("PRI_REL", "C2acc1y", "C2acc2y", "C2acc3y",     
            "C2acc4y", "C2acc5y", "CRPBI_PARENT",
            "C1acc1y", "C1acc2y",  "C1acc3y",
            "C1acc4y", "C1acc5y", "SecCare"))

################################################################################
##################### Compute SI & NSSI broad dummy scores 

# If youth reports past or present, then the value is 1, otherwise 0 (or NA)
df$SuIdActT <- ifelse(
  df$SuIdActpt == 1 | df$SuIdActpr == 1, 
  1,
  ifelse(is.na(df$SuIdActpt) & is.na(df$SuIdActpr), NA, 0)
)

df$SuIdPasT <- ifelse(
  df$SuIdPaspt == 1 | df$SuIdPaspr == 1, 
  1,
  ifelse(is.na(df$SuIdPaspt) & is.na(df$SuIdPaspr), NA, 0)
)

df$NSSIT <- ifelse(
  df$NSSIpt == 1 | df$NSSIpr == 1, 
  1,
  ifelse(is.na(df$NSSIpt) & is.na(df$NSSIpr), NA, 0)
)


################################################################################
#####################  SUBSET INTO WAVES 

table(df$eventname)

SRNR <- subset(df, eventname == "ses-00S") 
W1 <- subset(df, eventname == "ses-00A") 
W1_5 <- subset(df, eventname == "ses-00M") 
W2 <- subset(df, eventname == "ses-01A") 
W2_5 <- subset(df, eventname == "ses-01M") 
W3 <- subset(df, eventname == "ses-02A") 
W3_5 <- subset(df, eventname == "ses-02M") 
W4 <- subset(df, eventname == "ses-03A") 
W4_5 <- subset(df, eventname == "ses-03M") 
W5 <- subset(df, eventname == "ses-04A") 
W5_5 <- subset(df, eventname == "ses-04M") 
W6 <- subset(df, eventname == "ses-05A") 
W6_5 <- subset(df, eventname == "ses-05M") 
W7 <- subset(df, eventname == "ses-06A") 

# REMOVE ANY WAVE df THAT HAS ZERO OBSERVATIONS 

rm(list = keep(ls(), ~ {
  obj <- get(.x)
  is.data.frame(obj) && nrow(obj) == 0
}))


##################### REMOVE EMPTY COLUMNS FROM WAVE dfs

# LIST ALL DFs
WAVE_DFS <- ls(envir = .GlobalEnv)[
  sapply(ls(envir = .GlobalEnv), function(x) {
    is.data.frame(get(x, envir = .GlobalEnv)) &&
      !x %in% c("df", "NAME_DIC") # ADD NAMES HERE IF OTHER NON-WAVE DFS are PRESENT! 
  })
]

# REMOVE EMPTY COLUMNS FROM THE LIST 

for (nm in WAVE_DFS) {
  assign(
    nm,
    janitor::remove_empty(get(nm, envir = .GlobalEnv), which = "cols"),
    envir = .GlobalEnv
  )
}

##################### REMOVE TIME-INVARIANT & IRRELEVANT VARIABLES

# PULL DFs INTO A LIST 
all_dats <- lapply(WAVE_DFS, get, envir = .GlobalEnv)

# REMOVE ANY VARIABLES THAT SHOULD BE EXCLUDED ENTIRELY (FROM ALL WAVES)
all_dats <- lapply(all_dats, function(df) {
  dplyr::select(df, -any_of("eventname")) # ADD MORE HERE AS NECESSARY
})

# REMOVE INVARIANT VARIABLES FROM LATTER WAVES

## LIST THEM 
vars_to_remove <- c("SiteID", "Y_SEX", "FamilyID", "Y_HISP", "ppensity") # ADD VARIABLES AS NECESSARY

## REMOVE THEM 
if (length(all_dats) > 1) {
  all_dats[-1] <- lapply(all_dats[-1], function(df) {
    dplyr::select(df, -any_of(vars_to_remove))
  })
}

## RETURN LIST TO ENVIRONMENT 
for (i in seq_along(WAVE_DFS)) {
  assign(WAVE_DFS[i], all_dats[[i]], envir = .GlobalEnv)
}


##################### ADD WAVE NUMBER TO VARIABLES 


# LIST OF VARIABLES THAT DO NOT NEED WAVE SPECIFIERS
no_wave_vars <- c("subID", "SiteID", "Y_SEX", "FamilyID", "ppensity", 
                  "Y_HISP") # ADD AS NEEDED 

# MAP WAVES TO NUMBERS 
wave_map <- c(
  W1=1, W1_5=2, W2=3, W2_5=4, W3=5, W3_5=6, W4=7, W4_5=8,
  W5=9, W5_5=10, W6=11, W6_5=12, W7=13, SRNR=0
)

# MAKE A LIST OF WAVES TO CONSIDER 
df_names <- names(wave_map)

# RENAME DF VARIABLES DYNAMICALLY 
for (df_name in df_names) {
  # Check if data frame exists and has rows
  if (exists(df_name, envir = .GlobalEnv)) {
    df <- get(df_name, envir = .GlobalEnv)
    if (nrow(df) > 0) {
      wave_num <- wave_map[df_name]
      
      df <- df %>%
        rename_with(
          .cols = -any_of(no_wave_vars),  # Exclude variables that should not change
          .fn = ~ paste0(., "_", wave_num)  # Append wave number
        )
      
      assign(df_name, df, envir = .GlobalEnv)
    }
  }
}


#################### IDENTIFY DUPLICATE IDs WITHIN WAVES


# INITIALIZE A LIST 
duplicates_report <- list()

for (df_name in WAVE_DFS) {
  # Only process if data frame exists and has rows
  if (exists(df_name, envir = .GlobalEnv)) {
    df <- get(df_name, envir = .GlobalEnv)
    if (nrow(df) > 0 && "subID" %in% names(df)) {
      # Find duplicated subIDs
      dup_ids <- df$subID[duplicated(df$subID)]
      
      # Store results if any duplicates found
      if (length(dup_ids) > 0) {
        duplicates_report[[df_name]] <- unique(dup_ids)
      }
    }
  }
}

# Report results
if (length(duplicates_report) == 0) {
  message("No duplicate participant IDs were found within any wave.")
} else {
  message("⚠️ The following wave data frames have duplicate IDs that must be resolved prior to proceeding:")
  print(duplicates_report)
}

#################### LIST ONLY THE DFS THAT HAVE DATA OF INTEREST IN THEM 


# GET COLUMN COUNTS 
col_counts <- sapply(WAVE_DFS, function(df_name) ncol(get(df_name, envir = .GlobalEnv)))

# PRINT COLUMN COUNTS 
cat("\nColumn Counts for All Waves:\n")
cat("--------------------------------\n")
for (i in seq_along(col_counts)) {
  cat(sprintf("%s: %d columns\n", WAVE_DFS[i], col_counts[i]))
}

# KEEP DFs WITH MORE THAN A SPECIFIED NUMBER (CHANGE AS NEEDED)
VALID_WAVES <- WAVE_DFS[col_counts >= 5]

# Print which dataframes were kept for merging
if (length(VALID_WAVES) > 0) {
  cat("\nDataframes with 5 or more columns (ready for merging):\n")
  cat("------------------------------------------------------\n")
  cat(paste(VALID_WAVES, collapse = ", "), "\n")
} else {
  cat("\n⚠️ No wave data frames have 5 or more columns for merging.\n")
}


################################################################################
##################### MERGE DATA TOGETHER 


# Pull actual data frames from VALID_WAVES
ALL_WAVES_LIST <- lapply(VALID_WAVES, function(df_name) get(df_name, envir = .GlobalEnv))
names(ALL_WAVES_LIST) <- VALID_WAVES

# REMOVE DUPLICATE IDs (THIS IS A SAFETY STEP, IT SHOULD ALREADY BE RESOLVED!)
ALL_WAVES_LIST <- lapply(ALL_WAVES_LIST, function(df) {
  if ("subID" %in% names(df)) {
    df %>% distinct(subID, .keep_all = TRUE)
  } else {
    df
  }
})

# JOIN WAVES TOGETHER BY SubID
if (length(ALL_WAVES_LIST) > 0) {
  MERGED_WAVES <- ALL_WAVES_LIST[[1]]
  
  if (length(ALL_WAVES_LIST) > 1) {
    for (i in 2:length(ALL_WAVES_LIST)) {
      MERGED_WAVES <- full_join(MERGED_WAVES, ALL_WAVES_LIST[[i]], by = "subID")
    }
  }
  
  # SUMMARY OF THE MERGE
  if (exists("MERGED_WAVES")) {
    
    cat("Summary of Merged Waves (wide-format by subID):\n")
    cat("------------------------------------------------\n")
    
    # Rows and columns
    n_rows <- nrow(MERGED_WAVES)
    n_cols <- ncol(MERGED_WAVES)
    cat("Number of participants (rows):", n_rows, "\n")
    cat("Number of columns:", n_cols, "\n")
    
    # Check for duplicate subIDs
    n_dup_ids <- sum(duplicated(MERGED_WAVES$subID))
    if (n_dup_ids > 0) {
      cat("⚠️ Duplicate participant IDs detected:", n_dup_ids, "\n")
      dup_ids <- MERGED_WAVES$subID[duplicated(MERGED_WAVES$subID)]
      cat("Duplicated IDs:", paste(unique(dup_ids), collapse = ", "), "\n")
    } else {
      cat("No duplicate participant IDs detected.\n")
    }
    
    # Check for duplicated column names
    dup_cols <- names(MERGED_WAVES)[duplicated(names(MERGED_WAVES))]
    if (length(dup_cols) > 0) {
      cat("⚠️ Duplicate column names detected:", paste(dup_cols, collapse = ", "), "\n")
    } else {
      cat("No duplicate column names detected.\n")
    }
    
    # Optionally, report missing IDs from each wave
    if (exists("ALL_WAVES_LIST")) {
      for (wave_name in names(ALL_WAVES_LIST)) {
        missing_ids <- setdiff(MERGED_WAVES$subID, ALL_WAVES_LIST[[wave_name]]$subID)
        if (length(missing_ids) > 0) {
          cat(sprintf("⚠️ Wave '%s' missing %d participant(s)\n", wave_name, length(missing_ids)))
        }
      }
    }
    
  } else {
    cat("⚠️ MERGED_WAVES does not exist. Please perform the merge first.\n")
  }
}



###################### REMOVE WAVE 7 (_13, IT IS A PARTIAL WAVE - BE INTENTIAL WHEN KEPT

MERGED_WAVES_RED <- dplyr::select(MERGED_WAVES, -dplyr::ends_with("_13"))


################################################################################
##################### REORDER VARIABLES TO BE WITH LIKE VARIABLES


# DETECT INVARIANT VARIABLES (no _# suffix)
all_names <- names(MERGED_WAVES_RED)
id_names <- all_names[!grepl("_[0-9]+$", all_names)]  # invariant columns first


# EXTRACT VARIABLE KEYS AND WAVE NUMBERS
var_names <- setdiff(all_names, id_names)
var_key <- gsub("_[0-9]+$", "", var_names)
wave_num <- as.numeric(sub(".*_([0-9]+)$", "\\1", var_names))


# CREATE DATA FRAME FOR SORTING
var_df <- data.frame(var = var_names, key = var_key, wave = wave_num, stringsAsFactors = FALSE)


# SORT: first by key alphabetically, then by wave numerically
var_df <- var_df %>%
  arrange(key, wave)


# COMBINE IDENTIFIERS FIRST, THEN CLUSTERED VARIABLES
sorted_names <- c(id_names, var_df$var)


# REORDER THE DATAFRAME
ALL_WAVES_REORDER <- MERGED_WAVES_RED[, sorted_names, drop = FALSE]


# PRINT SUMMARY
cat("Reordering Summary:\n")
cat("------------------\n")
cat("Columns in MERGED_WAVES_RED (original):", ncol(MERGED_WAVES_RED), "\n")
cat("Columns in ALL_WAVES_REORDER (reordered):", ncol(ALL_WAVES_REORDER), "\n")
if (ncol(MERGED_WAVES_RED) != ncol(ALL_WAVES_REORDER)) {
  cat("⚠️ Warning: Column counts differ! Check for dropped or duplicated columns.\n")
}
cat("\nReordered Column Names:\n")
cat("----------------------\n")
print(names(ALL_WAVES_REORDER))


hist(ALL_WAVES_REORDER$rsSAL_Ramy_1)


################################################################################
##################### REDUCE THE SAMPLE TO NUCLEAR FAMILIES WITH IMAGING 

# vector of variables to check
DAD_VARS <- c("D_1acc_1", "D_1acc_3",  "D_1acc_7", "D_1acc_9",
                   "D_2acc_1", "D_2acc_3",  "D_2acc_7", "D_2acc_9",     
                   "D_3acc_1", "D_3acc_3", "D_3acc_7",  "D_3acc_9",      
                   "D_4acc_1",  "D_4acc_3", "D_4acc_7", "D_4acc_9",  "D_5acc_1",     
                   "D_5acc_3",  "D_5acc_7", "D_5acc_9", "D_acc_M_1", "D_acc_M_3",    
                   "D_acc_M_7", "D_acc_M_9")

MOM_VARS <- c("M_1acc_1", "M_1acc_3",  "M_1acc_7", "M_1acc_9",
               "M_2acc_1", "M_2acc_3",  "M_2acc_7", "M_2acc_9",     
               "M_3acc_1", "M_3acc_3", "M_3acc_7",  "M_3acc_9",      
               "M_4acc_1",  "M_4acc_3", "M_4acc_7", "M_4acc_9",  "M_5acc_1",     
               "M_5acc_3",  "M_5acc_7", "M_5acc_9", "M_acc_M_1", "M_acc_M_3",    
               "M_acc_M_7", "M_acc_M_9")


NUC_FAMS <- ALL_WAVES_REORDER %>%
  filter(
    if_any(all_of(DAD_VARS), ~ !is.na(.x)) &
      if_any(all_of(MOM_VARS), ~ !is.na(.x))
  )




T1_IMG <- c("rsSAL_Lamy_1",  "rsSAL_Lnacc_1", 
            "rsSAL_Ltha_1",  "rsSAL_Ramy_1",  
            "rsSAL_Rnacc_1", "rsSAL_Rtha_1",  
            "rsSAL_SAL_1")
T5_IMG <- c("rsSAL_Lamy_5",  "rsSAL_Lnacc_5", 
            "rsSAL_Ltha_5",  "rsSAL_Ramy_5",  
            "rsSAL_Rnacc_5", "rsSAL_Rtha_5",  
            "rsSAL_SAL_5")
T9_IMG <- c("rsSAL_Lamy_9",  "rsSAL_Lnacc_9", 
            "rsSAL_Ltha_9",  "rsSAL_Ramy_9",  
            "rsSAL_Rnacc_9", "rsSAL_Rtha_9",  
            "rsSAL_SAL_9")


NUC_FAMS <- NUC_FAMS %>%
  mutate(
    IMG_PT1T5 = as.integer(
      if_any(all_of(T1_IMG), ~ !is.na(.x)) &
        if_any(all_of(T5_IMG), ~ !is.na(.x))
    ),
    
    IMG_PT1T9 = as.integer(
      if_any(all_of(T1_IMG), ~ !is.na(.x)) &
        if_any(all_of(T9_IMG), ~ !is.na(.x))
    )
  )

## BREAKDOWN 

#T1 - T5 imagers 

table(NUC_FAMS$IMG_PT1T5)
# 0 = 4937 (49%)
# 1 = 5140 (51%)

#T1 - T9 imagers 
table(NUC_FAMS$IMG_PT1T9)
# 0 = 5831 (58%)
# 1 = 4246 (42%)


##ADD NUMERIC NUMBER ID

NUC_FAMS$NUMID <- row.names(NUC_FAMS)
names(NUC_FAMS)


################################################################################
##################### SAVE DATA 
names(NUC_FAMS)


# IMPUTE MISSING PPESNITY VALUE WITH THE MEAN 
NUC_FAMS$ppensity[NUC_FAMS$subID == "sub-P2G0PXCM"
 ] <- 691.3902

# FILL IN MISSING FAMILY IDs WITH SEQUENTIAL VALUES 
# get number missing 
n_missing <- sum(is.na(NUC_FAMS$FamilyID))

# Fill NAs with sequential values after the max
NUC_FAMS$FamilyID[is.na(NUC_FAMS$FamilyID)] <-
  seq(from = max(NUC_FAMS$FamilyID, na.rm = TRUE) + 1,
      length.out = n_missing)


FILENAME <- "ABCD_SAPPENFIELD_FATHERS_RR_FINAL"

write.csv(NUC_FAMS, paste0(FILENAME, ".csv"), row.names=FALSE, na="")

#prepareMplusData(ALL_WAVES_REORDER,"ABCD_HORM_METH_12.29.25.dat")


