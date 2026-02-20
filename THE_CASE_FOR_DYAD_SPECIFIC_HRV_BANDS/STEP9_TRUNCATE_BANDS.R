################################################################################
######################### LIBRARY SHTUFFF
library(dplyr)
library(DescTools)
library(ggplot2)
library(tidyr)


################################################################################
############################ SET DIRECTORY & LOAD FOB SUMMARY 

work_dir <- "/home/cjh37695/WAVELET_WORKING/TEST_RESP"
setwd(work_dir)

# SUMMARY OF BAND SELECTION 
FOB_SUM <- read.csv("FINAL_FREQUENCY_OPTIMIZED_BAND_SELECTION.csv")

################################################################################
##################### SET DIRECTORIES FOR ALL FILES 

FOB_COH <- file.path("FOB_COHERENCE_FULL/")
FOB_PHS <- file.path("FOB_PHASE_FULL/")
RSB_COH <- file.path("RSB_COHERENCE_FULL/")
RSB_PHS <- file.path("RSB_PHASE_FULL/")
HF_COH <- file.path("HF_COHERENCE_FULL/")
HF_PHS <- file.path("HF_PHASE_FULL/")
LF_COH <- file.path("LF_COHERENCE_FULL/")
LF_PHS <- file.path("LF_PHASE_FULL/")



################################################################################
##################### HELPERS

extract_id4 <- function(filepath) {
  fn <- basename(filepath)
  id <- regmatches(fn, regexpr("\\d{4}", fn))
  if (length(id) == 0 || id == "") NA_character_ else id
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

# Downsample factor for 4 Hz -> 1.25s bins:
# 1.25 seconds * 4 Hz = 5 samples per bin
## This was chosen because the lowest period is 2.5, so to detect
## changes at that level it should be 1/2 speed, i.e., 1.25 second bins
BIN_SIZE <- 5L

# Arithmetic mean per bin with NA-adaptive behavior (NA if all NA in bin)
bin_mean_na <- function(x, bin_size = BIN_SIZE) {
  x <- as.numeric(x)
  n <- length(x)
  if (n == 0) return(numeric(0))
  n_bins <- n %/% bin_size  # for 1200, with bin=5 => 240
  if (n_bins < 1) stop("Not enough samples to bin.")
  x <- x[seq_len(n_bins * bin_size)]
  mat <- matrix(x, nrow = bin_size, ncol = n_bins)
  # colMeans with na.rm=TRUE gives NaN when all NA; convert to NA
  out <- colMeans(mat, na.rm = TRUE)
  out[is.nan(out)] <- NA_real_
  out
}

# Circular mean per bin (radians) with NA-adaptive behavior
# Returns NA if all NA in bin; if 1 value, returns that value.
bin_circmean_na <- function(theta, bin_size = BIN_SIZE) {
  theta <- as.numeric(theta)
  n <- length(theta)
  if (n == 0) return(numeric(0))
  n_bins <- n %/% bin_size
  if (n_bins < 1) stop("Not enough samples to bin.")
  theta <- theta[seq_len(n_bins * bin_size)]
  mat <- matrix(theta, nrow = bin_size, ncol = n_bins)
  
  out <- rep(NA_real_, n_bins)
  for (j in seq_len(n_bins)) {
    v <- mat[, j]
    v <- v[!is.na(v)]
    if (length(v) == 0) {
      out[j] <- NA_real_
    } else if (length(v) == 1) {
      out[j] <- v
    } else {
      # circular mean via atan2(mean(sin), mean(cos))
      out[j] <- atan2(mean(sin(v)), mean(cos(v)))
    }
  }
  out
}

# Downsample a whole wavelet df: keep Period column, downsample columns 2:ncol
downsample_wavelet_df <- function(df, mode = c("coherence", "phase"), bin_size = BIN_SIZE) {
  mode <- match.arg(mode)
  
  if (ncol(df) < 3) stop("DF has too few columns: expected Period + >=2 time columns.")
  period_col <- df[[1]]
  X <- df[, -1, drop = FALSE]
  
  # Ensure we have 1200 samples (or 1201 etc). We'll truncate to full bins.
  n_time <- ncol(X)
  n_bins <- n_time %/% bin_size
  if (n_bins < 1) stop("Not enough timepoints to bin with bin_size=", bin_size)
  
  X <- X[, seq_len(n_bins * bin_size), drop = FALSE]  # preserve time ordering
  # Downsample row-wise (period-wise)
  out_mat <- matrix(NA_real_, nrow = nrow(X), ncol = n_bins)
  
  if (mode == "coherence") {
    for (i in seq_len(nrow(X))) {
      out_mat[i, ] <- bin_mean_na(as.numeric(X[i, ]), bin_size = bin_size)
    }
  } else { # phase
    for (i in seq_len(nrow(X))) {
      out_mat[i, ] <- bin_circmean_na(as.numeric(X[i, ]), bin_size = bin_size)
    }
  }
  
  # Make new time column names (bin centers in seconds)
  # If original is 4Hz over 300s, typical time step is 0.25s.
  # Bin size=5 => 1.25s bins. We'll label by bin center: 0.625, 1.875, ...
  t_step <- 0.25
  bin_width <- bin_size * t_step              # 1.25
  centers <- (seq_len(n_bins) - 0.5) * bin_width
  colnames(out_mat) <- sprintf("t_%0.3f", centers)
  
  out <- data.frame(Period = period_col, out_mat, check.names = FALSE)
  out
}


################################################################################
##################### BATCH PROCESSOR (one directory -> one output directory)

process_dir_downsample <- function(in_dir, out_dir, kind = c("coherence", "phase"), prefix) {
  kind <- match.arg(kind)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  files <- list.files(in_dir, pattern = "\\.csv$", full.names = TRUE)
  if (length(files) == 0) {
    warning("No CSV files found in: ", in_dir)
    return(invisible(NULL))
  }
  
  for (f in files) {
    id <- extract_id4(f)
    if (is.na(id)) {
      warning("No 4-digit ID found in filename: ", f)
      next
    }
    
    df <- read_csv_fast(f)
    
    df_ds <- downsample_wavelet_df(df, mode = kind, bin_size = BIN_SIZE)
    
    out_name <- paste0(prefix, id, "_DOWNSAMPLE.csv")
    out_path <- file.path(out_dir, out_name)
    write_csv_fast(df_ds, out_path)
  }
  
  invisible(TRUE)
}


################################################################################
##################### OUTPUT DIRECTORIES (create as needed)

FOB_COH_DS <- "FOB_COHERENCE_DOWNSAMPLED/"
FOB_PHS_DS <- "FOB_PHASE_DOWNSAMPLED/"
RSB_COH_DS <- "RSB_COHERENCE_DOWNSAMPLED/"
RSB_PHS_DS <- "RSB_PHASE_DOWNSAMPLED/"
HF_COH_DS  <- "HF_COHERENCE_DOWNSAMPLED/"
HF_PHS_DS  <- "HF_PHASE_DOWNSAMPLED/"
LF_COH_DS  <- "LF_COHERENCE_DOWNSAMPLED/"
LF_PHS_DS  <- "LF_PHASE_DOWNSAMPLED/"


################################################################################
##################### RUN ALL 8 BATCHES

# Coherence (arithmetic mean)
process_dir_downsample(FOB_COH, FOB_COH_DS, kind = "coherence", prefix = "FOB_COHERENCE_")
print("Done downsampling FOB Coherence.")
process_dir_downsample(RSB_COH, RSB_COH_DS, kind = "coherence", prefix = "RSB_COHERENCE_")
print("Done downsampling RSB Coherence.")
process_dir_downsample(HF_COH,  HF_COH_DS,  kind = "coherence", prefix = "HF_COHERENCE_")
print("Done downsampling HF Coherence.")
process_dir_downsample(LF_COH,  LF_COH_DS,  kind = "coherence", prefix = "LF_COHERENCE_")
print("Done downsampling LF Coherence.")

# Phase (circular mean)
process_dir_downsample(FOB_PHS, FOB_PHS_DS, kind = "phase", prefix = "FOB_PHASE_")
print("Done downsampling FOB Phase")
process_dir_downsample(RSB_PHS, RSB_PHS_DS, kind = "phase", prefix = "RSB_PHASE_")
print("Done downsampling RSB Phase")
process_dir_downsample(HF_PHS,  HF_PHS_DS,  kind = "phase", prefix = "HF_PHASE_")
print("Done downsampling HF Phase")
process_dir_downsample(LF_PHS,  LF_PHS_DS,  kind = "phase", prefix = "LF_PHASE_")
print("Done downsampling LF Phase")

print("Done downsampling HRV Wavelet Frequency Bands.")
################################################################################
##################### REDUCE TO A SINGLE TIME SERIES VIA MEANS

################################################################################
##################### COLLAPSE PERIOD-ROWS -> SINGLE MEAN SERIES (COH + PHASE)
# - Input: *DOWNSAMPLED* directories (FOB/RSB/HF/LF x COH/PHASE)
# - Output:
#   1) Plots for every file: all period-row series (grey) + mean series (colored),
#      saved to: [BAND]_[COHERENCE|PHASE]_MEAN_PLOTS/
#      filename: [BAND]_[COHERENCE|PHASE]_[ID]_MEAN_PLOT.png
#   2) Two sample-wide dataframes (one coherence, one phase) with one row per ID x BAND:
#        ID, BAND, period_start, period_end, HZ_start, HZ_end, [time series columns...]
#      saved as:
#        SAMPLE_TRUNC_COHERENCE_DW_SERIES.csv
#        SAMPLE_TRUNC_PHASE_DW_SERIES.csv
#
# Notes:
# - Handles COI edge NAs: means are NA only if an entire column is NA.
# - Time columns are kept exactly as in the downsampled files (dynamic; 240, 420, etc).
################################################################################

# HF AND LOW BANDS (Hz)
HF_HZ <- c(0.15, 0.40)
LF_HZ <- c(0.04, 0.15)

################################################################################
##################### INPUT DIRECTORIES (DOWNSAMPLED)

FOB_COH_DS <- "FOB_COHERENCE_DOWNSAMPLED/"
FOB_PHS_DS <- "FOB_PHASE_DOWNSAMPLED/"
RSB_COH_DS <- "RSB_COHERENCE_DOWNSAMPLED/"
RSB_PHS_DS <- "RSB_PHASE_DOWNSAMPLED/"
HF_COH_DS  <- "HF_COHERENCE_DOWNSAMPLED/"
HF_PHS_DS  <- "HF_PHASE_DOWNSAMPLED/"
LF_COH_DS  <- "LF_COHERENCE_DOWNSAMPLED/"
LF_PHS_DS  <- "LF_PHASE_DOWNSAMPLED/"

################################################################################
##################### HELPERS 

extract_id4 <- function(filepath) {
  fn <- basename(filepath)
  id <- regmatches(fn, regexpr("\\d{4}", fn))
  if (length(id) == 0 || id == "") NA_character_ else id
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
    readr::write_csv(df)
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
names(FOB_SUM)
# Metadata builder: period_start/end in seconds; HZ_start/end in Hz
get_band_meta <- function(id, band) {
  band <- toupper(band)
  
  if (band == "FOB") {
    row <- FOB_SUM[FOB_SUM$ID == as.integer(id) | FOB_SUM$ID == id, , drop = FALSE]
    if (nrow(row) == 0) return(NULL)
    p1 <- as.numeric(row$FOB_period_start[1])
    p2 <- as.numeric(row$FOB_period_end[1])
  } else if (band == "RSB") {
    row <- FOB_SUM[FOB_SUM$ID == as.integer(id) | FOB_SUM$ID == id, , drop = FALSE]
    if (nrow(row) == 0) return(NULL)
    p1 <- as.numeric(row$RSB_period_start[1])
    p2 <- as.numeric(row$RSB_period_end[1])
  } else if (band == "HF") {
    # Convert Hz band to Period band in seconds:
    # higher Hz = smaller period
    p1 <- 1 / HF_HZ[2]   # smaller period
    p2 <- 1 / HF_HZ[1]   # larger period
  } else if (band == "LF") {
    p1 <- 1 / LF_HZ[2]
    p2 <- 1 / LF_HZ[1]
  } else {
    return(NULL)
  }
  
  # ensure period_start <= period_end
  period_start <- min(p1, p2, na.rm = TRUE)
  period_end   <- max(p1, p2, na.rm = TRUE)
  
  # Hz boundaries (low Hz corresponds to larger period)
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
##################### CORE: process one file -> mean row + plot saved

process_one_file_mean <- function(f, band, kind = c("COHERENCE", "PHASE"), plot_dir) {
  kind <- match.arg(kind)
  id <- extract_id4(f)
  if (is.na(id)) {
    warning("No 4-digit ID found in filename: ", f)
    return(NULL)
  }
  
  df <- read_csv_fast(f)
  if (ncol(df) < 3) {
    warning("Too few columns in: ", f)
    return(NULL)
  }
  
  # Period is col 1 (rows represent periods)
  Period <- df[[1]]
  X <- as.matrix(df[, -1, drop = FALSE])  # periods x time
  mode(X) <- "numeric"
  
  time_names <- colnames(df)[-1]
  if (is.null(time_names) || any(time_names == "")) {
    time_names <- paste0("t", seq_len(ncol(X)))
    colnames(X) <- time_names
  }
  
  time_axis <- get_time_axis(time_names)
  
  # Compute band mean series (unweighted)
  if (kind == "COHERENCE") {
    mean_series <- col_mean_na(X)
  } else {
    mean_series <- col_circmean_na(X)
  }
  
  # Build 1-row output with dynamic number of time columns
  meta <- get_band_meta(id, band)
  if (is.null(meta)) {
    warning("No metadata found for ID ", id, " band ", band, " (skipping).")
    return(NULL)
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
  
  # Append time series columns
  ts_df <- as.data.frame(t(mean_series))
  colnames(ts_df) <- time_names
  out_row <- cbind(out_row, ts_df)
  
  # ---- Plot: all rows + mean overlay ----
  n_period <- nrow(X)
  n_time   <- ncol(X)
  
  orig_long <- data.frame(
    Period = rep(Period, times = n_time),
    Time   = rep(time_axis, each = n_period),
    Value  = as.vector(X)
  )
  
  mean_line <- data.frame(
    Time  = time_axis,
    Value = mean_series
  )
  
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
      x = "Time (downsampled bins)",
      y = ifelse(kind == "COHERENCE", "Coherence", "Phase (radians)"),
      title = paste0(band, " ", kind, " — period rows + mean (ID ", id, ")")
    ) +
    theme_minimal(base_family = "Arial")
  
  dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
  plot_path <- file.path(plot_dir, paste0(band, "_", kind, "_", id, "_MEAN_PLOT.png"))
  ggsave(plot_path, plot = p, width = 11, height = 5, dpi = 300, bg = "white")
  
  out_row
}

################################################################################
##################### BATCH: process a directory -> returns stacked df

process_dir_mean <- function(in_dir, band, kind = c("COHERENCE", "PHASE")) {
  kind <- match.arg(kind)
  
  plot_dir <- paste0(band, "_", kind, "_MEAN_PLOTS/")
  files <- list.files(in_dir, pattern = "\\.csv$", full.names = TRUE)
  
  if (length(files) == 0) {
    warning("No CSV files found in: ", in_dir)
    return(NULL)
  }
  
  rows <- vector("list", length(files))
  for (i in seq_along(files)) {
    rows[[i]] <- process_one_file_mean(
      f = files[i],
      band = band,
      kind = kind,
      plot_dir = plot_dir
    )
  }
  
  rows <- rows[!vapply(rows, is.null, logical(1))]
  if (length(rows) == 0) return(NULL)
  
  # Bind with fill if time columns differ (they shouldn't, but just in case)
  out <- dplyr::bind_rows(rows)
  out
}

################################################################################
##################### RUN ALL 8 (build two sample-wide dfs)

# --- COHERENCE ---
coh_FOB <- process_dir_mean(FOB_COH_DS, band = "FOB", kind = "COHERENCE")
print("Done collapsing FOB COHERENCE")
coh_RSB <- process_dir_mean(RSB_COH_DS, band = "RSB", kind = "COHERENCE")
print("Done collapsing RSB COHERENCE")
coh_HF  <- process_dir_mean(HF_COH_DS,  band = "HF",  kind = "COHERENCE")
print("Done collapsing HF COHERENCE")
coh_LF  <- process_dir_mean(LF_COH_DS,  band = "LF",  kind = "COHERENCE")
print("Done collapsing LF COHERENCE")

SAMPLE_COH <- dplyr::bind_rows(coh_FOB, coh_RSB, coh_HF, coh_LF)
print("Collapsed COHERENCE DF COMBINED")

# --- PHASE ---
phs_FOB <- process_dir_mean(FOB_PHS_DS, band = "FOB", kind = "PHASE")
print("Done collapsing FOB PHASE")
phs_RSB <- process_dir_mean(RSB_PHS_DS, band = "RSB", kind = "PHASE")
print("Done collapsing RSB PHASE")
phs_HF  <- process_dir_mean(HF_PHS_DS,  band = "HF",  kind = "PHASE")
print("Done collapsing HF PHASE")
phs_LF  <- process_dir_mean(LF_PHS_DS,  band = "LF",  kind = "PHASE")
print("Done collapsing LF PHASE")

SAMPLE_PHS <- dplyr::bind_rows(phs_FOB, phs_RSB, phs_HF, phs_LF)
print("Collapsed PHASE DF COMBINED")


################################################################################
##################### SAVE SAMPLE-WIDE DFS

write.csv(SAMPLE_COH, "SAMPLE_TRUNC_COHERENCE_DS_SERIES.csv")
write.csv(SAMPLE_PHS, "SAMPLE_TRUNC_PHASE_DS_SERIES.csv")

cat("DONE.\n",
    "Saved: SAMPLE_TRUNC_COHERENCE_DW_SERIES.csv (rows=", nrow(SAMPLE_COH), ")\n",
    "Saved: SAMPLE_TRUNC_PHASE_DW_SERIES.csv (rows=", nrow(SAMPLE_PHS), ")\n", sep = "")

