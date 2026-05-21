#Library
library(dplyr)
library(psych)
library(stringr)
library(purrr)
library(MplusAutomation)

########## SET WORKING DIRECTORY, OUTPUT DIRECTORY, AND FILE NAME ########## 

setwd("/home/cjh37695/ABCD_PROJECTS/")

OUT_DIR <- file.path(getwd(), "SUICIDALITY_SOCIAL_FACTORS")

# create base filename
FILENAME <- "ABCD_SSP"

# get today's date appended
DATE_STAMP <- format(Sys.Date(), "%m.%d.%y")


########## LOAD IN DATA AND GCDS NAMING DICTIONARY ########## 

df <- read.csv("SUICIDALITY_SOCIAL_FACTORS/RAW_DATA/dataset.csv", na.strings = c("999", "555", "888", "777", "-999"))

NAME_DIC <- read.csv("G-CDS_ABCD.VARIABLE.NAME_DICTIONARY.csv")


################################################################################
#####################  RENAME THE VARIABLES 

#### FIND ANY MISSING VARIABLES FROM REQUEST MISSING IN DICTIONARY ####

not_found <- setdiff(names(df), NAME_DIC$OG_NAME_6on)

# PRINT ANY VARIABLES NOT FOUND IN THE DATA DICTIONARY 

if (length(not_found) > 0) {
  cat("The following variable names are not found in the Center Data Dictionary\n")
  print(not_found)
} else {
  cat("All variable names in the data frame are found in the Center Data Dictionary\n")
}

## PROCEED TO RE-NAMING VARIABLES WHEN ALL ARE FOUND 

# CREATE NAMING DICTIONARY 
name_mapping <- setNames(NAME_DIC$G.CDS_NAME, NAME_DIC$OG_NAME_6on)

#RENAME BASED ON MATCHING IN THE DICTIONARY 
names(df) <- ifelse(names(df) %in% names(name_mapping), 
                                        name_mapping[names(df)], names(df))

## VIEW RESULTS 
names(df)

### RECODE VISIT DATE TO BE USEFUL (if present) 

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

#### Resting State functional connectivity #####

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

##### DIFFUSION TENSOR IMAGING  #####

if (any(grepl("^dti", colnames(df)))) {
  
  # Rows explicitly failing QC (not included, and not missing)
  affected_rows_dti <- df$dti_INC != 1 & !is.na(df$dti_INC)
  
  # Log excluded IDs + total count
  qc_na_subids$dmri <- df$subID[affected_rows_dti]
  qc_na_counts$dmri <- sum(affected_rows_dti)
  
  # NEW: excluded count per wave
  qc_na_counts_by_wave$dmri <- table(df$eventname[affected_rows_dti])
  
  # Add exclusion flag + total count column(s)
  df <- df %>%
    mutate(
      dti_excluded   = affected_rows_dti,
      dti_excluded_n = qc_na_counts$dmri
    ) %>%
    mutate(across(
      starts_with("dti"),
      ~case_when(
        dti_INC != 1 ~ NA_real_,
        TRUE ~ .
      )
    ))
}

# TODO: replicate the same imaging inclusion pattern for NBACK / MID once those blocks are implemented


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

names(df)


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
W7_5 <- subset(df, eventname == "ses-07A") 

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

##################### REMOVE TIME-INVARIANT & IRRELEVANT VARIABLES IN WAVE DFs

# PULL DFs INTO A LIST 
all_dats <- lapply(WAVE_DFS, get, envir = .GlobalEnv)

# REMOVE ANY VARIABLES THAT SHOULD BE EXCLUDED ENTIRELY (FROM ALL WAVES)
all_dats <- lapply(all_dats, function(df) {
  dplyr::select(df, -any_of("eventname")) # ADD MORE HERE IF NECESSARY
})

# REMOVE INVARIANT VARIABLES FROM LATTER WAVES 

## LIST THEM 
vars_to_remove <- c("SiteID", "Y_SEX", "FamilyID", "ppensity", 
                    "Y_HISP", "Y_RACE", "HiParEdu", "INCOME6L", "MarWrkSt") # ADD VARIABLES AS NECESSARY

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


##################### ADD WAVE NUMBER AS VARIABLE SUFFIX

# LIST OF VARIABLES THAT DO NOT NEED WAVE SPECIFIERS 
no_wave_vars <- c("subID", vars_to_remove)

# MAP WAVES TO NUMBERS 
wave_map <- c(
  W1=1, W1_5=2, W2=3, W2_5=4, W3=5, W3_5=6, W4=7, W4_5=8,
  W5=9, W5_5=10, W6=11, W6_5=12, W7=13, W7_5 = 14, SRNR=0
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

# KEEP DFs WITH MORE THAN A SPECIFIED NUMBER (ADJUST THRESHOLD AS NEEDED)
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


###################### REMOVE WAVE 7_5 (_14, IT IS A PARTIAL WAVE - BE INTENTIAL WHEN KEPT)

MERGED_WAVES_RED <- dplyr::select(MERGED_WAVES, -dplyr::ends_with("_14"))


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


################################################################################
##################### SAVE DATA & CODEBOOK

################## handling some ABCD nuisances 

# IMPUTE MISSING PPESNITY VALUE WITH THE MEAN 
ALL_WAVES_REORDER$ppensity[ALL_WAVES_REORDER$subID == "sub-P2G0PXCM"
] <- 691.3902

# FILL IN ANY MISSING FAMILY IDs WITH SEQUENTIAL VALUES 

# get number missing 
n_missing <- sum(is.na(ALL_WAVES_REORDER$FamilyID))

# Fill NAs with sequential values after the max
ALL_WAVES_REORDER$FamilyID[is.na(ALL_WAVES_REORDER$FamilyID)] <-
  seq(from = max(ALL_WAVES_REORDER$FamilyID, na.rm = TRUE) + 1,
      length.out = n_missing)


##################  CREATE FINAL CODEBOOK FOR FULL DATASET

# detect wave suffixes at the end of variable names
has_wave_suffix <- grepl("_[0-9]+$", names(ALL_WAVES_REORDER))

# Build variable-level summary from final FULL names
final_var_summary <- tibble::tibble(
  FULL_NAME = names(ALL_WAVES_REORDER),
  GCDS_NAME = gsub("_[0-9]+$", "", FULL_NAME),
  WAVE_SUFFIX = ifelse(
    grepl("_[0-9]+$", FULL_NAME),
    sub(".*(_[0-9]+)$", "\\1", FULL_NAME),
    NA_character_
  )
) %>%
  dplyr::group_by(GCDS_NAME) %>%
  dplyr::summarise(
    WAVES = ifelse(
      all(is.na(WAVE_SUFFIX)),
      "invariant",
      paste(sort(unique(stats::na.omit(WAVE_SUFFIX))), collapse = ", ")
    ),
    .groups = "drop"
  )

# Clean dictionary down to the fields needed for the codebook
name_dic_codebook <- NAME_DIC %>%
  dplyr::transmute(
    ABCD_NAME = OG_NAME_6on,
    GCDS_NAME = G.CDS_NAME,
    LABEL = LABEL
  ) %>%
  dplyr::distinct(GCDS_NAME, .keep_all = TRUE)

# Join final variable list to the dictionary
FINAL_CODEBOOK <- final_var_summary %>%
  dplyr::left_join(name_dic_codebook, by = "GCDS_NAME") %>%
  dplyr::select(ABCD_NAME, GCDS_NAME, WAVES, LABEL) %>%
  dplyr::arrange(GCDS_NAME)

# Check for any final variables not found in NAME_DIC
missing_from_codebook <- FINAL_CODEBOOK %>%
  dplyr::filter(is.na(ABCD_NAME))

if (nrow(missing_from_codebook) > 0) {
  cat("\nThe following variables are in ALL_WAVES_REORDER but were not found in NAME_DIC:\n")
  print(missing_from_codebook$GCDS_NAME)
} else {
  cat("\nAll variables in ALL_WAVES_REORDER were found in NAME_DIC and included in the codebook.\n")
}


################## SAVE THE FINAL TRANSFORMED DATA AND CODEBOOK 

# save dataframe
write.csv(
  ALL_WAVES_REORDER, 
  file.path(OUT_DIR, paste0(FILENAME, "_BASE_DATA_", DATE_STAMP, ".csv")), 
  row.names=FALSE, na="")

# Save codebook
write.csv(
  FINAL_CODEBOOK,
  file.path(OUT_DIR, paste0(FILENAME, "_CODEBOOK_", DATE_STAMP, ".csv")),
  row.names = FALSE,
  na = ""
)



