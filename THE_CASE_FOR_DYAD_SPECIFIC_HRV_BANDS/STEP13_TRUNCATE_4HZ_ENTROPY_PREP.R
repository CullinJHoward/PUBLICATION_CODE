### LIBRARY
library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
library(RespirAnalyzer)
library(ggplot2)

## We need to test entropy on the non-downsampled data. This script accomplishes 
## by truncating the oringal full 4Hz data so entropy can be computed without 
## needlessly smoothening it out via downsampling means. 

################################################################################
############################ SET DIRECTORY & LOAD FOB SUMMARY

work_dir <- "/home/cjh37695/WAVELET_WORKING/TEST_RESP"
setwd(work_dir)

FOB_SUM <- read.csv("FINAL_FREQUENCY_OPTIMIZED_BAND_SELECTION.csv")


################################################################################
##################### INPUT DIRECTORIES (FULL; NO DOWNSAMPLING)

## PAIRED DATA 
FOB_COH <- file.path("FOB_COHERENCE_FULL/")
FOB_PHS <- file.path("FOB_PHASE_FULL/")
RSB_COH <- file.path("RSB_COHERENCE_FULL/")
RSB_PHS <- file.path("RSB_PHASE_FULL/")
HF_COH  <- file.path("HF_COHERENCE_FULL/")
HF_PHS  <- file.path("HF_PHASE_FULL/")
LF_COH  <- file.path("LF_COHERENCE_FULL/")
LF_PHS  <- file.path("LF_PHASE_FULL/")

## RANDOM DYADS 
CV_FOB_COH <- file.path("RESAMPLING_VALIDATION/CV_FOB_COHERENCE_FULL/")
CV_FOB_PHS <- file.path("RESAMPLING_VALIDATION/CV_FOB_PHASE_FULL/")
CV_HF_COH  <- file.path("RESAMPLING_VALIDATION/CV_HF_COHERENCE_FULL/")
CV_HF_PHS  <- file.path("RESAMPLING_VALIDATION/CV_HF_PHASE_FULL/")
CV_LF_COH  <- file.path("RESAMPLING_VALIDATION/CV_LF_COHERENCE_FULL/")
CV_LF_PHS  <- file.path("RESAMPLING_VALIDATION/CV_LF_PHASE_FULL/")

################################################################################
##################### OUTPUT DIRECTORIES / FILENAMES

OUT_COH_FILE <- "SAMPLE_TRUNC_COHERENCE_FULL_4Hz_SERIES.csv"
OUT_PHS_FILE <- "SAMPLE_TRUNC_PHASE_FULL_4Hz_SERIES.csv"
OUT_CV_COH_FILE <- "RESAMPLING_VALIDATION/CV_SAMPLE_TRUNC_COHERENCE_FULL_4Hz_SERIES.csv"
OUT_CV_PHS_FILE <- "RESAMPLING_VALIDATION/CV_SAMPLE_TRUNC_PHASE_FULL_4Hz_SERIES.csv"

################################################################################
##################### HELPERS


extract_id <- function(filepath) {
  fn <- basename(filepath)
  
  # Extract all digit sequences anywhere in the filename
  id <- stringr::str_extract(fn, "\\d+")
  
  if (is.na(id)) {
    warning("No numeric ID found in filename: ", filepath)
    return(NA_character_)
  }
  
  # Pad to 4 digits
  sprintf("%04d", as.integer(id))
}

read_csv_fast <- function(path) {
  if (requireNamespace("readr", quietly = TRUE)) {
    readr::read_csv(path, show_col_types = FALSE)
  } else {
    read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  }
}

write_csv_fast <- function(df, path) {
  if (requireNamespace("readr", quietly = TRUE)) {
    readr::write_csv(df, path)
  } else {
    write.csv(df, path, row.names = FALSE)
  }
}

# NA-adaptive arithmetic mean by column (returns NA if all NA in that column)
col_mean_na <- function(M) {
  out <- colMeans(M, na.rm = TRUE)
  out[is.nan(out)] <- NA_real_
  out
}

# NA-adaptive circular mean by column (radians; returns NA if all NA)
col_circmean_na <- function(M) {
  M <- as.matrix(M)
  mode(M) <- "numeric"
  out <- rep(NA_real_, ncol(M))
  for (j in seq_len(ncol(M))) {
    v <- M[, j]
    v <- v[!is.na(v)]
    if (length(v) == 0) {
      out[j] <- NA_real_
    } else if (length(v) == 1) {
      out[j] <- v
    } else {
      out[j] <- atan2(mean(sin(v)), mean(cos(v)))
    }
  }
  out
}

# Pull time axis from column names like "t_0.625"; fallback to 1..N if absent
get_time_axis <- function(time_names) {
  x <- suppressWarnings(as.numeric(sub("^t[_\\.]", "", time_names)))
  if (all(is.na(x))) seq_along(time_names) else x
}

################################################################################
##################### BAND META 

HF_HZ <- c(0.15, 0.40)
LF_HZ <- c(0.04, 0.15)

get_band_meta <- function(id, band) {
  band <- toupper(band)
  
  # helper: NA meta row
  na_meta <- data.frame(
    period_start = NA_real_,
    period_end   = NA_real_,
    HZ_start     = NA_real_,
    HZ_end       = NA_real_
  )
  
  if (band == "FOB") {
    row <- FOB_SUM[FOB_SUM$ID == id, , drop = FALSE]
    if (nrow(row) == 0) return(na_meta)
    p1 <- as.numeric(row$FOB_period_start[1])
    p2 <- as.numeric(row$FOB_period_end[1])
    
  } else if (band == "RSB") {
    row <- FOB_SUM[FOB_SUM$ID == id, , drop = FALSE]
    if (nrow(row) == 0) return(na_meta)
    p1 <- as.numeric(row$RSB_period_start[1])
    p2 <- as.numeric(row$RSB_period_end[1])
    
  } else if (band == "HF") {
    p1 <- 1 / HF_HZ[2]
    p2 <- 1 / HF_HZ[1]
    
  } else if (band == "LF") {
    p1 <- 1 / LF_HZ[2]
    p2 <- 1 / LF_HZ[1]
    
  } else {
    return(na_meta)
  }
  
  period_start <- min(p1, p2, na.rm = TRUE)
  period_end   <- max(p1, p2, na.rm = TRUE)
  
  # If p1/p2 were both NA, min/max become Inf/-Inf; guard that:
  if (!is.finite(period_start) || !is.finite(period_end)) return(na_meta)
  
  HZ_start <- 1 / period_end
  HZ_end   <- 1 / period_start
  
  data.frame(
    period_start = period_start,
    period_end   = period_end,
    HZ_start     = HZ_start,
    HZ_end       = HZ_end
  )
}


################################################################################
##################### CORE: process one FULL file -> mean row + plot

process_one_file_mean_full <- function(f, band, kind = c("COHERENCE", "PHASE"), plot_dir) {
  kind <- match.arg(kind)
  id <- extract_id(f)
  if (is.na(id)) {
    warning("No 4-digit ID found in filename: ", f)
    return(NULL)
  }
  
  df <- read_csv_fast(f)
  if (ncol(df) < 3) {
    warning("Too few columns in: ", f)
    return(NULL)
  }
  
  Period <- df[[1]]
  X <- as.matrix(df[, -1, drop = FALSE])
  mode(X) <- "numeric"
  
  time_names <- colnames(df)[-1]
  if (is.null(time_names) || any(time_names == "")) {
    time_names <- paste0("t", seq_len(ncol(X)))
    colnames(X) <- time_names
  }
  
  time_axis <- get_time_axis(time_names)
  
  # Mean series across period-rows (column-wise)
  mean_series <- if (kind == "COHERENCE") col_mean_na(X) else col_circmean_na(X)
  
  meta <- get_band_meta(id, band)
  if (all(is.na(meta))) {
    warning("No metadata found for ID ", id, " band ", band, " (keeping row; metadata = NA).")
  }
  
  
  out_row <- data.frame(
    ID = id,
    BAND = band,
    period_start = meta$period_start,
    period_end   = meta$period_end,
    HZ_start     = meta$HZ_start,
    HZ_end       = meta$HZ_end,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  
  ts_df <- as.data.frame(t(mean_series))
  colnames(ts_df) <- time_names
  out_row <- cbind(out_row, ts_df)
  
  # Plot: all period rows (grey) + mean (colored)
  n_period <- nrow(X)
  n_time   <- ncol(X)
  
  orig_long <- data.frame(
    Period = rep(Period, times = n_time),
    Time   = rep(time_axis, each = n_period),
    Value  = as.vector(X)
  )
  
  mean_line <- data.frame(Time = time_axis, Value = mean_series)
  
  p <- ggplot() +
    geom_line(
      data = orig_long,
      aes(x = Time, y = Value, group = Period),
      color = "grey70",
      linewidth = 0.35,
      alpha = 0.7,
      na.rm = TRUE
    ) +
    geom_line(
      data = mean_line,
      aes(x = Time, y = Value),
      linewidth = 1.1,
      na.rm = TRUE,
      inherit.aes = FALSE
    ) +
    labs(
      x = "Time",
      y = ifelse(kind == "COHERENCE", "Coherence", "Phase (radians)"),
      title = paste0(band, " ", kind, " — period rows + mean (FULL; ID ", id, ")")
    ) +
    theme_minimal(base_family = "Arial")
  
  dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
  plot_path <- file.path(plot_dir, paste0(band, "_", kind, "_", id, "_MEAN_FULL_4Hz_PLOT_FULL.png"))
  ggsave(plot_path, plot = p, width = 11, height = 5, dpi = 300, bg = "white")
  
  out_row
}

################################################################################
##################### BATCH: process one FULL directory -> stacked df

process_dir_mean_full <- function(in_dir, band, kind = c("COHERENCE", "PHASE")) {
  kind <- match.arg(kind)
  
  plot_dir <- paste0(band, "_", kind, "_FULL_4Hz_MEAN_PLOTS_FULL/")
  files <- list.files(in_dir, pattern = "\\.csv$", full.names = TRUE)
  
  if (length(files) == 0) {
    warning("No CSV files found in: ", in_dir)
    return(NULL)
  }
  
  rows <- vector("list", length(files))
  for (i in seq_along(files)) {
    rows[[i]] <- process_one_file_mean_full(
      f = files[i],
      band = band,
      kind = kind,
      plot_dir = plot_dir
    )
  }
  
  rows <- rows[!vapply(rows, is.null, logical(1))]
  if (length(rows) == 0) return(NULL)
  
  dplyr::bind_rows(rows)
}

################################################################################
##################### RUN ALL 8 FOR PAIRED AND RANDOM DYADS AND SAVE

## PAIRED DYADS
# --- COHERENCE (arithmetic mean) ---
coh_FOB <- process_dir_mean_full(FOB_COH, band = "FOB", kind = "COHERENCE")
print("Done collapsing FOB COHERENCE (FULL)")
coh_RSB <- process_dir_mean_full(RSB_COH, band = "RSB", kind = "COHERENCE")
print("Done collapsing RSB COHERENCE (FULL)")
coh_HF  <- process_dir_mean_full(HF_COH,  band = "HF",  kind = "COHERENCE")
print("Done collapsing HF COHERENCE (FULL)")
coh_LF  <- process_dir_mean_full(LF_COH,  band = "LF",  kind = "COHERENCE")
print("Done collapsing LF COHERENCE (FULL)")

SAMPLE_COH_FULL <- dplyr::bind_rows(coh_FOB, coh_RSB, coh_HF, coh_LF)
write_csv_fast(SAMPLE_COH_FULL, OUT_COH_FILE)
print("Done collapsing ALL COHERENCE data")

# --- PHASE (circular mean) ---
phs_FOB <- process_dir_mean_full(FOB_PHS, band = "FOB", kind = "PHASE")
print("Done collapsing FOB PHASE (FULL)")
phs_RSB <- process_dir_mean_full(RSB_PHS, band = "RSB", kind = "PHASE")
print("Done collapsing RSB PHASE (FULL)")
phs_HF  <- process_dir_mean_full(HF_PHS,  band = "HF",  kind = "PHASE")
print("Done collapsing HF PHASE (FULL)")
phs_LF  <- process_dir_mean_full(LF_PHS,  band = "LF",  kind = "PHASE")
print("Done collapsing LF PHASE (FULL)")

SAMPLE_PHS_FULL <- dplyr::bind_rows(phs_FOB, phs_RSB, phs_HF, phs_LF)
write_csv_fast(SAMPLE_PHS_FULL, OUT_PHS_FILE)
print("Done collapsing ALL PHASE data")

## RANDOM DYADS 

# --- COHERENCE (arithmetic mean) ---
CV_coh_FOB <- process_dir_mean_full(CV_FOB_COH, band = "FOB", kind = "COHERENCE")
print("Done collapsing FOB COHERENCE (FULL)")
CV_coh_HF  <- process_dir_mean_full(CV_HF_COH,  band = "HF",  kind = "COHERENCE")
print("Done collapsing HF COHERENCE (FULL)")
CV_coh_LF  <- process_dir_mean_full(CV_LF_COH,  band = "LF",  kind = "COHERENCE")
print("Done collapsing LF COHERENCE (FULL)")

CV_SAMPLE_COH_FULL <- dplyr::bind_rows(CV_coh_FOB, CV_coh_HF, CV_coh_LF)
write_csv_fast(CV_SAMPLE_COH_FULL, OUT_CV_COH_FILE)
print("Done collapsing ALL Random Paired COHERENCE data")

# --- PHASE (circular mean) ---
CV_phs_FOB <- process_dir_mean_full(CV_FOB_PHS, band = "FOB", kind = "PHASE")
print("Done collapsing FOB PHASE (FULL)")
CV_phs_HF  <- process_dir_mean_full(CV_HF_PHS,  band = "HF",  kind = "PHASE")
print("Done collapsing HF PHASE (FULL)")
CV_phs_LF  <- process_dir_mean_full(CV_LF_PHS,  band = "LF",  kind = "PHASE")
print("Done collapsing LF PHASE (FULL)")

CV_SAMPLE_PHS_FULL <- dplyr::bind_rows(CV_phs_FOB, CV_phs_HF, CV_phs_LF)
write_csv_fast(CV_SAMPLE_PHS_FULL, OUT_CV_PHS_FILE)
print("Done collapsing ALL Random Paired PHASE data")
