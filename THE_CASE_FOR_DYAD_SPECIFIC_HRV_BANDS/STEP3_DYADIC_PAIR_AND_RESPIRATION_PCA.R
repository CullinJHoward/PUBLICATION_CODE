### LIBRARY
library(dplyr)
library(ggplot2)

# Set working directory
setwd("/home/cjh37695/WAVELET_WORKING")
work_dir <- "/home/cjh37695/WAVELET_WORKING"

################################################################################
####################### LOAD IN IBI DATA

# Define the subfolder path within work_dir
folder <- file.path(work_dir, "ALL_CLEAN_IBI_4HZ")

# Get the list of .csv files in the subfolder
csv_files <- list.files(path = folder, pattern = "_CLEAN.csv$", full.names = TRUE)

# Function to load, clean, and rename CSV files
load_and_clean_csv <- function(file_name) {
  # Extract the base file name (without the full path)
  base_name <- basename(file_name)
  
  # Extract the prefix (C_ or P_) and retain underscore
  prefix <- sub("(_.*)", "_", base_name)  # Keeps "C_" or "P_"
  
  # Extract the family number (digits before "_CLEAN.csv")
  family_number <- sub(".*_(\\d+)_CLEAN.csv$", "\\1", base_name)
  
  # Construct the new variable name with underscore
  new_name <- paste0(prefix, family_number)  # Example: "C_0001" or "P_0002"
  
  # Read the CSV file
  data <- read.csv(file_name, stringsAsFactors = FALSE)
  
  # Remove columns that are entirely NA
  data <- data[, colSums(!is.na(data)) > 0]
  
  # Assign the cleaned data frame to the new name in the environment
  assign(new_name, data, envir = .GlobalEnv)
}

# Apply the function to all CSV files
lapply(csv_files, load_and_clean_csv)

################################################################################
################### IDENTIFY DYADS AND INDIVIDUAL DATA

# Get all object names in the environment that are data frames
all_dfs <- ls()[sapply(ls(), function(x) is.data.frame(get(x)))]

# Extract family IDs and role (Parent or Child)
df_info <- data.frame(
  df_name = all_dfs,
  role = sub("_(\\d{4})$", "", all_dfs),  # Extract P_ or C_
  family_id = sub("^[PC]_(\\d{4}).*$", "\\1", all_dfs),  # Extract the 4-digit family ID
  stringsAsFactors = FALSE
)

# Identify unique family IDs
unique_families <- unique(df_info$family_id)

# Create lists for dyads and individuals
dyads <- c()
individuals <- c()

# Loop through each family ID
for (fam_id in unique_families) {
  # Find all data frames related to this family
  family_dfs <- df_info[df_info$family_id == fam_id, ]
  
  # Check if both a parent and child exist
  has_parent <- any(grepl("^P_", family_dfs$df_name))
  has_child  <- any(grepl("^C_", family_dfs$df_name))
  
  if (has_parent & has_child) {
    dyads <- c(dyads, fam_id)  # Add to dyads list
  } else {
    individuals <- c(individuals, fam_id)  # Add to individuals list
  }
}

## DYADS
print(dyads)

## INDIVIDUALS
print(individuals)

###################### REMOVE NON-DYAD DATA

# Remove data frames associated with families in the individuals list
for (fam_id in individuals) {
  # Find all data frames that match this family ID
  dfs_to_remove <- ls()[grepl(paste0("^[PC]_?", fam_id), ls())]
  
  # Remove them from the environment
  rm(list = dfs_to_remove, envir = .GlobalEnv)
}

################################################################################
############### ENSURE ROW NUMBER MATCHES BETWEEN DYAD MEMBERS

# Loop through each family ID in the dyads list
for (family_id in dyads) {
  parent_df_name <- paste0("P_", family_id)
  child_df_name  <- paste0("C_", family_id)
  
  if (exists(parent_df_name, envir = .GlobalEnv) && exists(child_df_name, envir = .GlobalEnv)) {
    df_parent <- get(parent_df_name, envir = .GlobalEnv)
    df_child  <- get(child_df_name, envir = .GlobalEnv)
    
    min_rows <- min(nrow(df_parent), nrow(df_child))
    df_parent <- df_parent[1:min_rows, , drop = FALSE]
    df_child  <- df_child[1:min_rows, , drop = FALSE]
    
    assign(parent_df_name, df_parent, envir = .GlobalEnv)
    assign(child_df_name,  df_child,  envir = .GlobalEnv)
  }
}

################################################################################
################# MERGE EACH PERSONS IBI AND RESP DATAFRAMES

merge_dyad_data <- function() {
  for (family_id in dyads) {
    parent_df_name <- paste0("P_", family_id)
    child_df_name  <- paste0("C_", family_id)
    
    if (exists(parent_df_name, envir = .GlobalEnv) && exists(child_df_name, envir = .GlobalEnv)) {
      df_parent <- get(parent_df_name, envir = .GlobalEnv)
      df_child  <- get(child_df_name,  envir = .GlobalEnv)
      
      min_rows <- min(nrow(df_parent), nrow(df_child))
      df_parent <- df_parent[1:min_rows, , drop = FALSE]
      df_child  <- df_child[1:min_rows,  , drop = FALSE]
      
      # Remove redundant columns from child
      if ("ID" %in% colnames(df_parent) && "ID" %in% colnames(df_child) && all(df_parent$ID == df_child$ID)) {
        df_child <- df_child[, !(colnames(df_child) == "ID")]
      }
      if ("TIME_SEC" %in% colnames(df_parent) && "TIME_SEC" %in% colnames(df_child) && all(df_parent$TIME_SEC == df_child$TIME_SEC)) {
        df_child <- df_child[, !(colnames(df_child) == "TIME_SEC")]
      }
      
      merged_df <- cbind(df_parent, df_child)
      new_name <- paste0(family_id, "_DYAD")
      assign(new_name, merged_df, envir = .GlobalEnv)
    }
  }
}

merge_dyad_data()

################################################################################
############### DYAD-LEVEL RESP PCA (PC1) + SUMMARY + PLOTS

# --- settings ---
pca_plot_dir <- file.path(work_dir, "RESPIRATORY_PCA_PLOTS")
if (!dir.exists(pca_plot_dir)) dir.create(pca_plot_dir, recursive = TRUE)

# Initialize summary df
RESP_PCA_SUMMARY <- data.frame(
  DYAD_ID = character(),
  N_TIMEPOINTS = integer(),
  PC1_VAR_EXPLAINED = numeric(),
  PC1_LOADING_P_RESP = numeric(),
  PC1_LOADING_C_RESP = numeric(),
  PC1_COR_P_RESP = numeric(),
  PC1_COR_C_RESP = numeric(),
  stringsAsFactors = FALSE
)

# Helper: safe z-score (keeps NAs)
zscore <- function(x) {
  mu <- mean(x, na.rm = TRUE)
  sdv <- sd(x, na.rm = TRUE)
  if (is.na(sdv) || sdv == 0) return(rep(NA_real_, length(x)))
  (x - mu) / sdv
}

# Loop through each dyad DF
dyad_dfs <- ls()[grepl("_DYAD$", ls())]

for (dyad_name in dyad_dfs) {
  
  df <- get(dyad_name)
  
  # Find the parent/child resp columns in the merged dyad df
  # Accept either raw or detrended naming variants:
  p_resp_col <- grep("^P_.*RESP(_D)?$", names(df), value = TRUE)
  c_resp_col <- grep("^C_.*RESP(_D)?$", names(df), value = TRUE)
  
  # If multiple candidates exist, prefer detrended (_RESP_D)
  if (length(p_resp_col) > 1) {
    if (any(grepl("RESP_D$", p_resp_col))) p_resp_col <- p_resp_col[grepl("RESP_D$", p_resp_col)][1] else p_resp_col <- p_resp_col[1]
  }
  if (length(c_resp_col) > 1) {
    if (any(grepl("RESP_D$", c_resp_col))) c_resp_col <- c_resp_col[grepl("RESP_D$", c_resp_col)][1] else c_resp_col <- c_resp_col[1]
  }
  
  if (length(p_resp_col) != 1 || length(c_resp_col) != 1) {
    message("Skipping ", dyad_name, " (could not uniquely identify P/C respiration columns).")
    next
  }
  
  if (!"TIME_SEC" %in% names(df)) {
    message("Skipping ", dyad_name, " (TIME_SEC not found).")
    next
  }
  
  # Build PCA input with complete cases
  resp_p <- df[[p_resp_col]]
  resp_c <- df[[c_resp_col]]
  
  # Standardize each resp series before PCA (per dyad)
  resp_p_z <- zscore(resp_p)
  resp_c_z <- zscore(resp_c)
  
  pca_input <- data.frame(P_RESP_Z = resp_p_z, C_RESP_Z = resp_c_z)
  
  complete_idx <- complete.cases(pca_input)
  
  # Need enough points
  if (sum(complete_idx) < 10) {
    message("Skipping ", dyad_name, " (too few complete respiration observations for PCA).")
    next
  }
  
  # Run PCA on standardized respiration
  pca_fit <- prcomp(pca_input[complete_idx, ], center = FALSE, scale. = FALSE)
  
  # Extract PC1 time series for complete rows, then reinsert NAs
  pc1 <- rep(NA_real_, nrow(df))
  pc1[complete_idx] <- pca_fit$x[, 1]
  
  # Add PC1 to dyad df
  df$RESP_PC1 <- pc1
  
  # Summarize
  var_expl <- (pca_fit$sdev[1]^2) / sum(pca_fit$sdev^2)
  loadings <- pca_fit$rotation[, 1]
  loading_p <- unname(loadings["P_RESP_Z"])
  loading_c <- unname(loadings["C_RESP_Z"])
  
  # Correlations (using available points)
  cor_p <- suppressWarnings(cor(pc1, resp_p_z, use = "pairwise.complete.obs"))
  cor_c <- suppressWarnings(cor(pc1, resp_c_z, use = "pairwise.complete.obs"))
  
  dyad_id <- sub("_DYAD$", "", dyad_name)
  
  RESP_PCA_SUMMARY <- rbind(
    RESP_PCA_SUMMARY,
    data.frame(
      DYAD_ID = dyad_id,
      N_TIMEPOINTS = nrow(df),
      PC1_VAR_EXPLAINED = var_expl,
      PC1_LOADING_P_RESP = loading_p,
      PC1_LOADING_C_RESP = loading_c,
      PC1_COR_P_RESP = cor_p,
      PC1_COR_C_RESP = cor_c,
      stringsAsFactors = FALSE
    )
  )
  
  # Save plot: P resp, C resp, and PC1 (scaled for overlay)
  plot_df <- data.frame(
    TIME_SEC = df$TIME_SEC,
    P_RESP_Z = resp_p_z,
    C_RESP_Z = resp_c_z,
    RESP_PC1 = df$RESP_PC1
  )
  
  # For overlay readability, put everything on z-scale:
  # PC1 is already in z-like units; keep as-is.
  p <- ggplot(plot_df, aes(x = TIME_SEC)) +
    geom_line(aes(y = P_RESP_Z, color = "P_RESP (z)"), linewidth = 0.5, alpha = 0.8) +
    geom_line(aes(y = C_RESP_Z, color = "C_RESP (z)"), linewidth = 0.5, alpha = 0.8) +
    geom_line(aes(y = RESP_PC1,  color = "RESP_PC1"), linewidth = 0.8) +
    labs(
      title = paste0("Respiration + PC1: Dyad ", dyad_id),
      x = "Time (sec)",
      y = "Z-scored units",
      color = "Series"
    ) +
    theme_minimal()
  
  plot_path <- file.path(pca_plot_dir, paste0("RESP_PC1_", dyad_id, ".png"))
  ggsave(filename = plot_path, plot = p, width = 10, height = 4.5, dpi = 300)
  
  # Update dyad df back into environment
  assign(dyad_name, df, envir = .GlobalEnv)
  
  message("PCA + plot saved for dyad: ", dyad_id)
}

################################################################################
############# REMOVE WORKING DATAFRAMES (keep only DYAD dfs)

all_dfs <- ls()[sapply(ls(), function(x) is.data.frame(get(x)))]
dfs_to_remove <- all_dfs[!grepl("_DYAD$", all_dfs) & all_dfs != "RESP_PCA_SUMMARY"]
invisible(lapply(dfs_to_remove, function(df_name) rm(list = df_name, envir = .GlobalEnv)))


################################################################################
################ REMOVE RESP-REGRESSED IBI VARIABLES FROM DYAD DFS

# Identify all dyad-level data frames
dyad_dfs <- ls()[grepl("_DYAD$", ls())]

# Loop through each dyad df
for (df_name in dyad_dfs) {
  
  df <- get(df_name)
  
  # Identify columns to remove if present
  cols_to_remove <- intersect(
    c("P_IBI_D_R", "C_IBI_D_R"),
    names(df)
  )
  
  # Remove columns
  if (length(cols_to_remove) > 0) {
    df <- df[, !names(df) %in% cols_to_remove, drop = FALSE]
    assign(df_name, df, envir = .GlobalEnv)
    message("Removed ", paste(cols_to_remove, collapse = ", "),
            " from ", df_name)
  } else {
    message("No resp-regressed IBI variables found in ", df_name)
  }
}

rm(df)

################################################################################
############################ SAVE THE DYAD DATA FRAMES (+ PC1)

output_dir <- file.path(work_dir, "TEST_RESP/MERGED_DYADIC_IBI_4HZ")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

dyad_dfs <- ls()[grepl("_DYAD$", ls())]

for (df_name in dyad_dfs) {
  df <- get(df_name)
  file_path <- file.path(output_dir, paste0(df_name, ".csv"))
  write.csv(df, file_path, row.names = FALSE)
  message(paste("Saved:", file_path))
}

# Save PCA summary df
summary_path <- file.path(output_dir, "../RESP_PCA_SUMMARY.csv")
write.csv(RESP_PCA_SUMMARY, summary_path, row.names = FALSE)
message(paste("Saved:", summary_path))

