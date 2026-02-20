## LIBRARIES SCHTUFF
library(dplyr)
library(DescTools)
library(ggplot2)
library(slider)
library(tidyr)
library(stringr)
library(circular)
library(mgcv)
library(purrr)

################################################################################
############################ SET DIRECTORY & LOAD FOB SUMMARY 

work_dir <- "/home/cjh37695/WAVELET_WORKING/TEST_RESP"
setwd(work_dir)

FOB_SUM <- read.csv("FINAL_FREQUENCY_OPTIMIZED_BAND_SELECTION.csv")

## RENAME & FIX ID TO MATCH CROSS-DF FORMAT 

FOB_SUM <- FOB_SUM %>%
  rename(ID = "fam_id")

FOB_SUM$ID <- sprintf("%04d", FOB_SUM$ID)

## LIST DYADIC DATA TO PROCESS 

COHERENCE_FILES <- list.files(
  path = "/home/cjh37695/WAVELET_WORKING/TEST_RESP/COHERENCE_FILTERED",
  pattern = "\\.csv$",
  full.names = TRUE
)

PHASE_FILES <- list.files(
  path = "/home/cjh37695/WAVELET_WORKING/TEST_RESP/PHASE_FILTERED",
  pattern = "\\.csv$",
  full.names = TRUE
)


################################################################################
############################ DEFINE HELPER FUNCTIONS FOR EXTRACTION 


extract_id4 <- function(filepath) {
  fn <- basename(filepath)
  id <- regmatches(fn, regexpr("\\d{4}", fn))
  if (length(id) == 0 || id == "") NA_character_ else id
}

read_csv_fast <- function(path) {
  # readr is nicer if available; otherwise base R
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

################################################################################
############################ EXTRACT FULL FREQUENCY OPTIMZED BANDS


process_fob_list <- function(files, out_dir, out_prefix, FOB_SUM) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Ensure ID formats match (4-digit character)
  FOB_SUM$ID <- sprintf("%04d", as.integer(as.character(FOB_SUM$ID)))
  
  for (f in files) {
    id <- extract_id4(f)
    if (is.na(id)) {
      warning("No 4-digit ID found in filename: ", f)
      next
    }
    
    # Match row in FOB_SUM
    idx <- match(id, FOB_SUM$ID)
    if (is.na(idx)) {
      warning("ID ", id, " not found in FOB_SUM$ID; skipping: ", f)
      next
    }
    
    p_start <- FOB_SUM$period_start[idx]
    p_end   <- FOB_SUM$period_end[idx]
    
    # Read the CSV
    df <- read_csv_fast(f)
    
    # Get Period column
    if ("Period" %in% names(df)) {
      period <- df$Period
    } else {
      period <- df[[1]]
      # Optional: name it for clarity
      # names(df)[1] <- "Period"
    }
    
    # Coerce to numeric just in case
    period <- suppressWarnings(as.numeric(period))
    
    # Trim rows: >= start and <= end
    keep <- !is.na(period) & period >= p_start & period <= p_end
    df_trim <- df[keep, , drop = FALSE]
    
    # Write out
    out_path <- file.path(out_dir, paste0(out_prefix, id, "_FULL.csv"))
    write_csv_fast(df_trim, out_path)
  }
  
  invisible(TRUE)
}

# ---------- run ----------
process_fob_list(
  files = COHERENCE_FILES,
  out_dir = "/home/cjh37695/WAVELET_WORKING/TEST_RESP/FOB_COHERENCE_FULL",
  out_prefix = "FOB_COHERENCE_",
  FOB_SUM = FOB_SUM
)

process_fob_list(
  files = PHASE_FILES,
  out_dir = "/home/cjh37695/WAVELET_WORKING/TEST_RESP/FOB_PHASE_FULL",
  out_prefix = "FOB_PHASE_",
  FOB_SUM = FOB_SUM
)


################################################################################
#################### EXTRACT FULL HIGH AND LOW FREQUENCY BANDS 

# ---------- USER SETTINGS: define your fixed bands in Hz ----------
HF_BAND <- c(0.15, 0.40)  # CHANGE AS NEEDED 
LF_BAND <- c(0.04, 0.15)  # CHANGE AS NEEDED 


get_period_seconds <- function(df) {
  if ("Period" %in% names(df)) df$Period else df[[1]]
}

subset_by_band_hz <- function(df, band_hz, add_hz_col = TRUE) {
  period_s <- suppressWarnings(as.numeric(get_period_seconds(df)))
  hz <- 1 / period_s
  
  keep <- !is.na(hz) & is.finite(hz) & hz >= band_hz[1] & hz <= band_hz[2]
  out <- df[keep, , drop = FALSE]
  
  if (add_hz_col) {
    # Add Hz as a second column (keeps Period as col 1)
    out_hz <- hz[keep]
    out <- cbind(out[, 1, drop = FALSE], Hz = out_hz, out[, -1, drop = FALSE])
  }
  
  out
}

process_HFLF_files <- function(files, out_dir_hf, out_dir_lf,
                               prefix_hf, prefix_lf,
                               HF_BAND, LF_BAND) {
  dir.create(out_dir_hf, recursive = TRUE, showWarnings = FALSE)
  dir.create(out_dir_lf, recursive = TRUE, showWarnings = FALSE)
  
  for (f in files) {
    id <- extract_id4(f)
    if (is.na(id)) {
      warning("No 4-digit ID found in filename: ", f)
      next
    }
    
    df <- read_csv_fast(f)
    
    df_hf <- subset_by_band_hz(df, HF_BAND, add_hz_col = TRUE)
    df_lf <- subset_by_band_hz(df, LF_BAND, add_hz_col = TRUE)
    
    out_hf <- file.path(out_dir_hf, paste0(prefix_hf, id, "_FULL.csv"))
    out_lf <- file.path(out_dir_lf, paste0(prefix_lf, id, "_FULL.csv"))
    
    write_csv_fast(df_hf, out_hf)
    write_csv_fast(df_lf, out_lf)
  }
  
  invisible(TRUE)
}

# ---------- OUTPUT DIRECTORIES ----------
base_dir <- "/home/cjh37695/WAVELET_WORKING/TEST_RESP"

HF_COH_DIR <- file.path(base_dir, "HF_COHERENCE_FULL")
HF_PHA_DIR <- file.path(base_dir, "HF_PHASE_FULL")
LF_COH_DIR <- file.path(base_dir, "LF_COHERENCE_FULL")
LF_PHA_DIR <- file.path(base_dir, "LF_PHASE_FULL")

# ---------- run for COHERENCE ----------
process_HFLF_files(
  files = COHERENCE_FILES,
  out_dir_hf = HF_COH_DIR,
  out_dir_lf = LF_COH_DIR,
  prefix_hf = "HF_COHERENCE_",
  prefix_lf = "LF_COHERENCE_",
  HF_BAND = HF_BAND,
  LF_BAND = LF_BAND
)

# ---------- run for PHASE ----------
process_HFLF_files(
  files = PHASE_FILES,
  out_dir_hf = HF_PHA_DIR,
  out_dir_lf = LF_PHA_DIR,
  prefix_hf = "HF_PHASE_",
  prefix_lf = "LF_PHASE_",
  HF_BAND = HF_BAND,
  LF_BAND = LF_BAND
)


print("done extracfting FOB & HF/LF bands.")

################################################################################
#################### EXTRACT FULL RANDOM SURROGATE BANDS 

## ALL POSSIBLE PERIODS 
## ============================================================
## ROUTE B: match surrogate band by k (number of rows / scales)
## Prefer NO OVERLAP; if impossible, choose MINIMUM OVERLAP
## Plot on LOG-frequency axis; plot FOB using GRID-SNAPPED bounds
## ============================================================

# ---- Dynamically derive AVAILABLE_PERIODS from one example file ----
# Prefer coherence file if available; otherwise phase file.
example_file <- dplyr::coalesce(COHERENCE_FILES[1], PHASE_FILES[1])

if (is.na(example_file) || is.null(example_file) || !file.exists(example_file)) {
  stop("Could not find an example coherence/phase CSV to derive AVAILABLE_PERIODS.")
}

df_ex <- read_csv_fast(example_file)

# Pull Period column robustly (named 'Period' or first column)
if ("Period" %in% names(df_ex)) {
  AVAILABLE_PERIODS <- df_ex$Period
} else {
  AVAILABLE_PERIODS <- df_ex[[1]]
}

AVAILABLE_PERIODS <- suppressWarnings(as.numeric(AVAILABLE_PERIODS))
AVAILABLE_PERIODS <- AVAILABLE_PERIODS[!is.na(AVAILABLE_PERIODS)]

# Ensure unique + sorted (typical for period grids)
AVAILABLE_PERIODS <- sort(unique(AVAILABLE_PERIODS))

# Basic sanity checks
if (length(AVAILABLE_PERIODS) < 5) {
  stop("AVAILABLE_PERIODS looks too short—check that the example CSV has a valid Period column.")
}

# (Optional) message for QC
message("Derived AVAILABLE_PERIODS from: ", basename(example_file),
        " (n = ", length(AVAILABLE_PERIODS), ")")


##ENSURE THESE MATCH wavelet output!
print(AVAILABLE_PERIODS)


set.seed(123)

P  <- AVAILABLE_PERIODS
nP <- length(P)

# Make sure band_k is integer
FOB_SUM$band_k <- as.integer(FOB_SUM$band_k)

if (any(is.na(FOB_SUM$band_k))) stop("FOB_SUM$band_k has NA(s).")
if (any(FOB_SUM$band_k < 1)) stop("FOB_SUM$band_k has values < 1.")
if (any(FOB_SUM$band_k > nP)) stop("Some band_k exceed available period grid length.")

# Snap rounded FOB bounds to nearest grid index
period_to_idx_nearest <- function(x, P) {
  if (is.na(x)) return(NA_integer_)
  which.min(abs(P - x))
}

FOB_SUM$fob_idx_lower <- vapply(FOB_SUM$period_start, period_to_idx_nearest, integer(1), P = P)
FOB_SUM$fob_idx_upper <- vapply(FOB_SUM$period_end,   period_to_idx_nearest, integer(1), P = P)

# Correct swap (temp variable)
swap <- FOB_SUM$fob_idx_lower > FOB_SUM$fob_idx_upper
tmp <- FOB_SUM$fob_idx_lower[swap]
FOB_SUM$fob_idx_lower[swap] <- FOB_SUM$fob_idx_upper[swap]
FOB_SUM$fob_idx_upper[swap] <- tmp

# Store snapped FOB periods (use these for plotting)
FOB_SUM$fob_period_lower_grid <- P[FOB_SUM$fob_idx_lower]
FOB_SUM$fob_period_upper_grid <- P[FOB_SUM$fob_idx_upper]

# ---- Generate RSB: contiguous window of SAME k, no overlap if possible; else minimal overlap ----
rsb_start_idx <- rep(NA_integer_, nrow(FOB_SUM))
rsb_overlap_k <- rep(NA_integer_, nrow(FOB_SUM))  # diagnostic: overlap rows

for (r in seq_len(nrow(FOB_SUM))) {
  k  <- FOB_SUM$band_k[r]
  fL <- FOB_SUM$fob_idx_lower[r]
  fU <- FOB_SUM$fob_idx_upper[r]
  
  if (any(is.na(c(k, fL, fU)))) next
  
  starts <- seq_len(nP - k + 1L)
  ends   <- starts + k - 1L
  
  # overlap in rows for each candidate window
  overlap_len <- pmax(0L, pmin(ends, fU) - pmax(starts, fL) + 1L)
  
  # prefer no overlap; else minimal overlap
  if (any(overlap_len == 0L)) {
    best_starts <- starts[overlap_len == 0L]
  } else {
    best_starts <- starts[overlap_len == min(overlap_len)]
  }
  
  s <- sample(best_starts, 1L)
  rsb_start_idx[r] <- s
  rsb_overlap_k[r] <- overlap_len[starts == s][1]
}

FOB_SUM$rsb_lower <- ifelse(is.na(rsb_start_idx), NA_real_, P[rsb_start_idx])
FOB_SUM$rsb_upper <- ifelse(is.na(rsb_start_idx), NA_real_, P[rsb_start_idx + FOB_SUM$band_k - 1L])

FOB_SUM$rsb_overlap_rows <- rsb_overlap_k  # optional diagnostic

# ---- prep plotting df ----
plot_df <- FOB_SUM
plot_df <- plot_df[
  !is.na(plot_df$fob_period_lower_grid) & !is.na(plot_df$fob_period_upper_grid) &
    !is.na(plot_df$rsb_lower) & !is.na(plot_df$rsb_upper),
]

# Dyad index (hide IDs)
plot_df$dyad_index <- seq_len(nrow(plot_df))

# Convert to Hz using GRID-SNAPPED FOB bounds and grid RSB bounds
plot_df$fob_hz_lo <- 1 / pmax(plot_df$fob_period_lower_grid, plot_df$fob_period_upper_grid)
plot_df$fob_hz_hi <- 1 / pmin(plot_df$fob_period_lower_grid, plot_df$fob_period_upper_grid)

plot_df$rsb_hz_lo <- 1 / pmax(plot_df$rsb_lower, plot_df$rsb_upper)
plot_df$rsb_hz_hi <- 1 / pmin(plot_df$rsb_lower, plot_df$rsb_upper)

# Stack FOB/RSB per dyad
offset <- 0.15
plot_df$y_fob <- plot_df$dyad_index + offset
plot_df$y_rsb <- plot_df$dyad_index - offset

# ---- plot (LOG frequency axis is the key for Route B) ----
p <- ggplot() +
  geom_segment(
    data = plot_df,
    aes(x = fob_hz_lo, xend = fob_hz_hi, y = y_fob, yend = y_fob, color = "FOB"),
    linewidth = 0.7, lineend = "round"
  ) +
  geom_segment(
    data = plot_df,
    aes(x = rsb_hz_lo, xend = rsb_hz_hi, y = y_rsb, yend = y_rsb, color = "RSB"),
    linewidth = 0.7, lineend = "round"
  ) +
  scale_color_manual(values = c("FOB" = "#1f77b4", "RSB" = "#ff7f0e")) +
  scale_x_log10() +
  labs(x = "Frequency (Hz, log scale)", y = "Dyad", color = NULL) +
  theme_minimal(base_family = "Arial") +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.position = "bottom"
  )

p


## SAVE THE PLOT 

ggsave(
  filename = "RSB_FOB_PAIRING.png",
  plot = p,
  width = 10,
  height = 8,
  units = "in",
  dpi = 300,
  bg = "white"
)

## SAVE THE RSB UPDATED DF 
names(FOB_SUM)

FOB_SUM_RED <- FOB_SUM %>%
  select(c("ID", "band_k", "period_start", 
           "period_end", "sig_scan_mean", "hz_start",
           "hz_end", "sample",
           "rsb_lower", "rsb_upper"))

FOB_SUM_RED_NAMED <- FOB_SUM_RED %>%
  rename(FOB_BAND_K = "band_k",
         FOB_period_start = "period_start",
         FOB_period_end = "period_end",
         FOB_BAND_SIG = "sig_scan_mean",
         FOB_hz_start = "hz_start",
         FOB_hz_end = "hz_end",
         RSB_period_start = "rsb_lower", 
         RSB_period_end = "rsb_upper")
   
write.csv(FOB_SUM_RED_NAMED, "FINAL_FREQUENCY_OPTIMIZED_BAND_SELECTION.csv", row.names = FALSE) 

################## EXTRACT THE RSB FILES 


process_rsb_list <- function(files, out_dir, out_prefix, FOB_SUM_RED_NAMED) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Ensure ID formats match (4-digit character)
  FOB_SUM_RED_NAMED$ID <- sprintf("%04d", as.integer(as.character(FOB_SUM_RED_NAMED$ID)))
  
  for (f in files) {
    id <- extract_id4(f)
    if (is.na(id)) {
      warning("No 4-digit ID found in filename: ", f)
      next
    }
    
    # Match row in FOB_SUM
    idx <- match(id, FOB_SUM_RED_NAMED$ID)
    if (is.na(idx)) {
      warning("ID ", id, " not found in FOB_SUM_RED_NAMED$ID; skipping: ", f)
      next
    }
    names(FOB_SUM_RED_NAMED)
    p_start <- FOB_SUM_RED_NAMED$RSB_period_start[idx]
    p_end   <- FOB_SUM_RED_NAMED$RSB_period_end[idx]
    
    # Read the CSV
    df <- read_csv_fast(f)
    
    # Get Period column
    if ("Period" %in% names(df)) {
      period <- df$Period
    } else {
      period <- df[[1]]
      # Optional: name it for clarity
      # names(df)[1] <- "Period"
    }
    
    # Coerce to numeric just in case
    period <- suppressWarnings(as.numeric(period))
    
    # Trim rows: >= start and <= end
    keep <- !is.na(period) & period >= p_start & period <= p_end
    df_trim <- df[keep, , drop = FALSE]
    
    # Write out
    out_path <- file.path(out_dir, paste0(out_prefix, id, "_FULL.csv"))
    write_csv_fast(df_trim, out_path)
  }
  
  invisible(TRUE)
}

# ---------- run ----------
process_rsb_list(
  files = COHERENCE_FILES,
  out_dir = "/home/cjh37695/WAVELET_WORKING/TEST_RESP/RSB_COHERENCE_FULL",
  out_prefix = "RSB_COHERENCE_",
  FOB_SUM = FOB_SUM_RED_NAMED
)

process_rsb_list(
  files = PHASE_FILES,
  out_dir = "/home/cjh37695/WAVELET_WORKING/TEST_RESP/RSB_PHASE_FULL",
  out_prefix = "RSB_PHASE_",
  FOB_SUM = FOB_SUM_RED_NAMED
)

print("Done extracting the RSB files")
