################################################################################
# BIWAVELET OUTPUT EXTRACTION
# - COI mask
# - Save COI-filtered COHERENCE (rsq) and PHASE (phase): keep ALL values inside COI
# - Build significance islands from obj$signif (threshold surface for rsq)
# - Compute RIDGE as local maxima across period for each timepoint
#     * RIDGE computation is done ONLY within (COI & significant-coherence) cells
# - Save:
#     COHERENCE_FILTERED/COHERENCE_<ID>.csv         (COI only)
#     PHASE_FILTERED/PHASE_<ID>.csv                 (COI only)
#     RIDGE_FILTERED/RIDGE_<ID>.csv                 (ridge from sig rsq within COI)
#     RIDGE_UNFILTERED_WCOI/RIDGE_COI_<ID>.csv      (ridge from raw rsq, no masks)
#     COI_VECTORS/COI_<ID>.csv                      (time + coi period boundary)
#     SIG_RIDGE_BY_PERIOD.csv                       (row-sum of ridge hits per period)
#
# Notes:
# - Your obj$signif is a coherence-threshold surface (NOT p-values), confirmed by:
#     summary(as.vector(obj$signif)) ~ 0.26-0.62
# - Ridge is NOT guaranteed >= 1 per timepoint under local-max + sig gating.
################################################################################

library(dplyr)

# ---------------------------
# WORKING DIRECTORY
# ---------------------------
setwd("/home/cjh37695/WAVELET_WORKING/TEST_RESP")
work_dir <- "/home/cjh37695/WAVELET_WORKING/TEST_RESP"

# ---------------------------
# INPUT: RDS DIRECTORY
# ---------------------------
coherence_dir <- "MATH_PWTC_DATA"
rds_files <- list.files(coherence_dir, pattern = "\\.RDS$", full.names = TRUE)

# ---------------------------
# USER SETTINGS
# ---------------------------
NEIGHBORHOOD <- 1            # local max across +/- NEIGHBORHOOD period bins
ALLOW_EDGE_MAXIMA <- FALSE   # allow first/last period bins to be maxima?

# ---------------------------
# OUTPUT DIRECTORIES
# ---------------------------
dir.create("COHERENCE_FILTERED", showWarnings = FALSE)
dir.create("PHASE_FILTERED", showWarnings = FALSE)
dir.create("RIDGE_FILTERED", showWarnings = FALSE)
dir.create("RIDGE_UNFILTERED_WCOI", showWarnings = FALSE)
dir.create("COI_VECTORS", showWarnings = FALSE)

coherence_output_dir <- "COHERENCE_FILTERED"
phase_output_dir <- "PHASE_FILTERED"
ridge_output_dir <- "RIDGE_FILTERED"
ridge_unfiltered_output_dir <- "RIDGE_UNFILTERED_WCOI"
coi_output_dir <- "COI_VECTORS"

# ---------------------------
# HELPERS
# ---------------------------
extract_family_id <- function(path) {
  base <- basename(path)
  gsub(".*?(\\d{4}).*", "\\1", base)
}

# Build a signif matrix matching rsq dims.
# Handles: scalar, period-vector, or flattened full matrix.
make_signif_matrix <- function(signif_obj, n_period, n_time) {
  if (is.null(signif_obj)) stop("Object missing 'signif' (needed for islands mask).")
  
  # Case 1: scalar threshold
  if (length(signif_obj) == 1) {
    return(matrix(signif_obj, nrow = n_period, ncol = n_time))
  }
  
  # Case 2: period-wise threshold vector
  if (length(signif_obj) == n_period) {
    return(matrix(signif_obj, nrow = n_period, ncol = n_time, byrow = FALSE))
  }
  
  # Case 3: full threshold matrix stored as a vector (flattened)
  if (length(signif_obj) == n_period * n_time) {
    # IMPORTANT: column-major reshape to match R's vectorization
    return(matrix(signif_obj, nrow = n_period, ncol = n_time, byrow = FALSE))
  }
  
  stop(paste0(
    "Unexpected length(obj$signif) = ", length(signif_obj),
    ". Expected 1, nrow(rsq)=", n_period,
    ", or nrow(rsq)*ncol(rsq)=", (n_period * n_time), "."
  ))
}

# Local maxima across period for each timepoint (column)
# strict: value must be > all available neighbors within +/- neighborhood
local_maxima_ridge <- function(mat, neighborhood = 1, allow_edges = FALSE) {
  nr <- nrow(mat); nc <- ncol(mat)
  ridge <- matrix(0, nrow = nr, ncol = nc)
  
  for (j in seq_len(nc)) {
    col <- mat[, j]
    if (all(is.na(col))) next
    
    idx <- seq_len(nr)
    if (!allow_edges) {
      idx <- idx[(idx > neighborhood) & (idx <= (nr - neighborhood))]
      if (length(idx) == 0) next
    }
    
    for (i in idx) {
      if (is.na(col[i])) next
      
      left_idx  <- max(1, i - neighborhood):(i - 1)
      right_idx <- (i + 1):min(nr, i + neighborhood)
      
      neighbors <- c(col[left_idx], col[right_idx])
      neighbors <- neighbors[!is.na(neighbors)]
      if (length(neighbors) == 0) next
      
      if (all(col[i] > neighbors)) ridge[i, j] <- 1
    }
  }
  
  ridge
}

combine_summary_df <- function(summary_list, period, metric_name) {
  df <- do.call(rbind, summary_list)
  df <- as.data.frame(df)
  df$family_id <- rownames(df)
  df <- df[, c("family_id", setdiff(names(df), "family_id"))]
  colnames(df)[-1] <- paste0(metric_name, "_", round(period, 3))
  df
}

# ---------------------------
# CONTAINERS
# ---------------------------
ridge_period_summary <- list()
last_period <- NULL

# ---------------------------
# MAIN LOOP
# ---------------------------
for (file in rds_files) {
  
  obj <- readRDS(file)
  family_id <- extract_family_id(file)
  
  # Required components
  if (is.null(obj$rsq))    stop("Missing obj$rsq in ", file)
  if (is.null(obj$period)) stop("Missing obj$period in ", file)
  if (is.null(obj$t))      stop("Missing obj$t in ", file)
  if (is.null(obj$coi))    stop("Missing obj$coi in ", file)
  if (is.null(obj$signif)) stop("Missing obj$signif in ", file)
  if (is.null(obj$phase))  stop("Missing obj$phase in ", file)
  
  rsq <- obj$rsq              # coherence matrix [period x time]
  phase <- obj$phase          # phase matrix [period x time]
  period <- obj$period        # length = nrow(rsq)
  time <- obj$t               # length = ncol(rsq)
  coi_period <- obj$coi       # length = ncol(rsq), in period units
  
  last_period <- period
  
  # Skip very short recordings (timepoints)
  if (ncol(rsq) < 300) {
    message("Dyad ", family_id, " skipped: fewer than 300 timepoints")
    next
  }
  
  # -------------------------
  # COI MASK (TRUE = inside COI)
  # -------------------------
  coi_mask <- outer(period, coi_period, "<=")
  
  # -------------------------
  # SIGNIFICANCE MASK FOR COHERENCE ISLANDS
  # (coherence significant if rsq >= threshold)
  # -------------------------
  signif_mat <- make_signif_matrix(obj$signif, nrow(rsq), ncol(rsq))
  stopifnot(all(dim(signif_mat) == dim(rsq)))
  
  sig_mask <- (rsq >= signif_mat)
  
  # Keep mask for ridge computation only
  keep_mask <- coi_mask & sig_mask
  
  # -------------------------
  # SAVE COHERENCE (COI ONLY)
  # -------------------------
  rsq_coi <- rsq
  rsq_coi[!coi_mask] <- NA
  
  write.csv(
    data.frame(Period = period, rsq_coi),
    file = file.path(coherence_output_dir, paste0("COHERENCE_", family_id, ".csv")),
    row.names = FALSE
  )
  
  # -------------------------
  # SAVE PHASE (COI ONLY)
  # -------------------------
  phase_coi <- phase
  phase_coi[!coi_mask] <- NA
  
  write.csv(
    data.frame(Period = period, phase_coi),
    file = file.path(phase_output_dir, paste0("PHASE_", family_id, ".csv")),
    row.names = FALSE
  )
  
  # -------------------------
  # SAVE COI VECTOR (time + COI boundary)
  # -------------------------
  write.csv(
    data.frame(Time = time, COI_Period = coi_period),
    file = file.path(coi_output_dir, paste0("COI_", family_id, ".csv")),
    row.names = FALSE
  )
  
  # -------------------------
  # RIDGE (UNFILTERED): local maxima on raw rsq (no COI/sig masking)
  # -------------------------
  ridge_unfiltered <- local_maxima_ridge(
    rsq,
    neighborhood = NEIGHBORHOOD,
    allow_edges = ALLOW_EDGE_MAXIMA
  )
  
  write.csv(
    data.frame(Period = period, ridge_unfiltered),
    file = file.path(ridge_unfiltered_output_dir, paste0("RIDGE_COI_", family_id, ".csv")),
    row.names = FALSE
  )
  
  # -------------------------
  # RIDGE (FILTERED): local maxima only within (COI & sig islands)
  # -------------------------
  rsq_for_ridge <- rsq
  rsq_for_ridge[!keep_mask] <- NA
  
  ridge_filtered <- local_maxima_ridge(
    rsq_for_ridge,
    neighborhood = NEIGHBORHOOD,
    allow_edges = ALLOW_EDGE_MAXIMA
  )
  
  write.csv(
    data.frame(Period = period, ridge_filtered),
    file = file.path(ridge_output_dir, paste0("RIDGE_", family_id, ".csv")),
    row.names = FALSE
  )
  
  # -------------------------
  # PERIOD SUMMARY: count ridge hits across time (row sums)
  # -------------------------
  ridge_counts_by_period <- rowSums(ridge_filtered == 1, na.rm = TRUE)
  ridge_period_summary[[family_id]] <- ridge_counts_by_period
  
  message("Processed Family ID: ", family_id)
}

# ---------------------------
# FINAL SUMMARY CSV
# ---------------------------
if (length(ridge_period_summary) > 0 && !is.null(last_period)) {
  ridge_period_df <- combine_summary_df(ridge_period_summary, last_period, "RIG")
  write.csv(ridge_period_df, "SIG_RIDGE_BY_PERIOD.csv", row.names = FALSE)
} else {
  warning("No dyads processed; SIG_RIDGE_BY_PERIOD.csv not written.")
}

