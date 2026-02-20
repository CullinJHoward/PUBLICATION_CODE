#!/usr/bin/env Rscript

###############################################################################
# process_one_dyad_pwtc.R
# Partial wavelet coherence using biwavelet::pwtc
# y  = C_IBI_D (detrended)
# x1 = P_IBI_D (detrended)
# x2 = RESP_PC1 (shared respiration covariate)
###############################################################################

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("Usage: Rscript process_one_dyad_pwtc.R <input_csv>")
}
input_csv <- args[1]

suppressMessages(library(biwavelet))
library(biwavelet)

# Dyad ID
filename <- basename(input_csv)
dyad_id  <- sub("_DYAD\\.csv$", "", filename)

# Output dirs
out_rds_dir  <- "MATH_PWTC_DATA"
out_plot_dir <- "MATH_PWTC_PLOTS"
if (!dir.exists(out_rds_dir))  dir.create(out_rds_dir,  recursive = TRUE)
if (!dir.exists(out_plot_dir)) dir.create(out_plot_dir, recursive = TRUE)

# Load
dat <- read.csv(input_csv, stringsAsFactors = FALSE)
dat <- dat[, colSums(!is.na(dat)) > 0]

# Drop very short recordings
if (nrow(dat) < 300) {
  message(paste("Dyad", dyad_id, "removed: fewer than 300 rows"))
  quit(status = 0)
}

# Trim to first 1200 rows
dat <- dat[1:min(1200, nrow(dat)), ]

# Required cols (detrended IBI + shared RESP PC1)
required_cols <- c("TIME_SEC", "C_IBI_D", "P_IBI_D", "RESP_PC1")
missing_cols <- setdiff(required_cols, names(dat))
if (length(missing_cols) > 0) {
  stop(paste("Dyad", dyad_id, "missing required columns:", paste(missing_cols, collapse = ", ")))
}

# NA check
na_counts <- sapply(dat[, required_cols], function(x) sum(is.na(x)))
if (any(na_counts > 0)) {
  stop(paste0(
    "Dyad ", dyad_id, " has NAs in required columns: ",
    paste(names(na_counts)[na_counts > 0], "=", na_counts[na_counts > 0], collapse = ", "),
    ". Fix upstream (interpolation/imputation) rather than dropping rows here."
  ))
}

# Regular sampling check
dt_vec <- diff(dat$TIME_SEC)
dt_med <- median(dt_vec)
if (any(abs(dt_vec - dt_med) > 1e-6)) {
  stop(paste0(
    "Dyad ", dyad_id, " has irregular TIME_SEC (dt not constant). ",
    "Upstream interpolation should prevent this—please verify."
  ))
}

# Build biwavelet matrices (n x 2): time, value
y  <- cbind(dat$TIME_SEC, dat$C_IBI_D)
x1 <- cbind(dat$TIME_SEC, dat$P_IBI_D)
x2 <- cbind(dat$TIME_SEC, dat$RESP_PC1)

# Match your WaveletComp resolution & period bounds as closely as possible
dj <- 0.05
lowerPeriod <- 2.5
upperPeriod <- 26

s0 <- lowerPeriod
J1 <- floor(log2(upperPeriod / s0) / dj)

# Partial wavelet coherence
PWTC <- pwtc(
  y = y, x1 = x1, x2 = x2,
  pad = TRUE,
  dj  = dj,
  s0  = s0,
  J1  = J1,
  mother    = "morlet",
  sig.level = 0.99,
  sig.test  = 0,
  nrands    = 2000, #not being used when sig.test = 0 (no monte carlo @ this level)
  quiet     = TRUE
)

# Save RDS
saveRDS(
  PWTC,
  file = file.path(out_rds_dir, paste0(dyad_id, "_MATH_PWTC.RDS"))
)

# Save scalogram with COI + significance contours
pdf(
  file = file.path(out_plot_dir, paste0(dyad_id, "_PartialWTC.pdf")),
  width = 8.5, height = 5.5
)

# ---- Create right-side white space so the colorbar doesn't overlap data ----
old_par <- par(no.readonly = TRUE)
on.exit(par(old_par), add = TRUE)

# mar = c(bottom, left, top, right)
par(mar = c(4.5, 4.5, 3.5, 8.5), xpd = NA)

# --- Native biwavelet plot (original axes unchanged) ---
plot(
  PWTC,
  plot.cb    = TRUE,
  plot.coi   = TRUE,
  plot.sig   = TRUE,
  plot.phase = FALSE,
  xlab = "Seconds",
  ylab = "Period",
  main = paste0("Dyad ", dyad_id, ": Partial WTC (C_IBI_D vs P_IBI_D | RESP_PC1)")
)

# Restore clipping behavior
par(xpd = FALSE)

dev.off()
cat(paste("Completed dyad:", dyad_id, "\n"))

