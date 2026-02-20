################################################################################
# FAMILY-SPECIFIC RIDGE BAND DETECTION PIPELINE (ALTERNATIVE)  — UPDATED
# Person-specific Scan-Null by Band Size (k) using within-dyad Rank-1 windows
# Author: Cullin
# Date: Jan 13 2026
#
# UPDATE (Jan 2026): Stage 5 band selection rule
#   - Keep the top (best) band per family under the existing ordering NO MATTER WHAT.
#   - Then keep ALL additional bands that are:
#       (a) significant on scan_p_mean (<= alpha_scan), AND
#       (b) non-overlapping with already-kept bands.
#
#   - Plots reflect all retained bands per family.
#   - Final CSV now retains multiple bands per dyad if they meet the criteria.
#
# OUTPUTS
#   1) ALL_CANDIDATE_BANDS_MEAN_SUM.csv
#   2) BEST_BANDS_BY_K_PER_FAMILY.csv
#   3) PERSON_SCAN_NULL_THRESHOLDS_BY_K.csv          (one row per fam_id x k)
#   4) PERSON_SCAN_NULL_DISTS_BY_K.rds               (list keyed by fam_id then k)
#   5) BEST_BANDS_WITH_PERSON_SCAN_P.csv
#   6) FREQUENCY_OPTIMIZED_BANDS_PERSON_NULL.csv     (now: primary + sig/non-overlap extras)
#   7) OPTIMIZED_BAND_PLOTS_PERSON_NULL/*.pdf
#   8) FINAL_FREQUENCY_OPTIMIZED_BAND_SELECTION.csv  (now: multiple bands per dyad possible)
################################################################################

## LIBRARIES
library(dplyr)
library(stringr)
library(ggplot2)

## PREP ENVIRONMENT
work_dir <- "/home/cjh37695/WAVELET_WORKING/TEST_RESP"
setwd(work_dir)

df <- read.csv("SIG_RIDGE_BY_PERIOD.csv")

rig_cols <- grep("^RIG_", names(df), value = TRUE)
if (!("family_id" %in% names(df))) stop("SIG_RIDGE_BY_PERIOD.csv must contain a 'family_id' column.")

R <- as.matrix(df[, rig_cols])
periods <- as.numeric(str_remove(rig_cols, "^RIG_"))

n_fam <- nrow(R)
n_per <- ncol(R)

################################################################################
#################### USER SETTINGS

set.seed(1234)

# Band width constraints
min_k <- 3
max_k <- floor(n_per / 2)
band_sizes <- min_k:max_k

# Inference / thresholds (person scan null; rank-1 under permutation; k-specific)
alpha_scan <- 0.05
n_mc_scan <- 10000   # PERSON-level scan null is expensive; start 2000-10000 then increase

# Selection settings (choices made on MEAN p-values)
max_bands_per_family <- Inf
top_n_bands_to_plot <- 10   # bump default because some families may have >5 now

# Plot settings
plot_dir <- file.path(work_dir, "OPTIMIZED_BAND_PLOTS_PERSON_NULL")
plot_width <- 6
plot_height <- 6
y_ticks <- c(4, 8, 16)

# Blue gradient for band ranks (auto-sized to top_n_bands_to_plot)
blue_palette <- colorRampPalette(c("#08306B", "#2171B5", "#6BAED6", "#C6DBEF"))(top_n_bands_to_plot)

################################################################################
#################### HELPER FUNCTIONS

# Upper-tail p-value against SORTED null sample
emp_p_upper <- function(obs, null_sorted) {
  N <- length(null_sorted)
  leq <- findInterval(obs, null_sorted, rightmost.closed = TRUE)  # # <= obs
  geq <- N - leq                                                 # # >= obs
  (geq + 1) / (N + 1)
}

# Fast window sums & means using cumulative sums
window_sums_means <- function(x, k) {
  n <- length(x)
  m <- n - k + 1
  if (m <= 0) return(list(sums = numeric(0), means = numeric(0)))
  cs <- c(0, cumsum(x))
  sums <- cs[(k + 1):(n + 1)] - cs[1:(n - k + 1)]
  means <- sums / k
  list(sums = sums, means = means)
}

# Overlap helper: TRUE if [s,e] overlaps any kept interval
overlaps_any <- function(s, e, kept_intervals) {
  if (length(kept_intervals) == 0) return(FALSE)
  any(vapply(kept_intervals, function(iv) {
    !(e < iv[1] || s > iv[2])
  }, logical(1)))
}

################################################################################
#################### STAGE 1: Enumerate ALL candidate bands (sum + mean only)

message("Stage 1: generating all candidate bands (sum + mean only)...")

band_list <- list()
row_id <- 1

for (i in seq_len(n_fam)) {
  ridge_i <- as.numeric(R[i, ])
  fam_id  <- df$family_id[i]
  
  for (k in band_sizes) {
    ws <- window_sums_means(ridge_i, k)
    sums <- ws$sums
    means <- ws$means
    max_start <- length(means)
    if (max_start == 0) next
    
    for (s in seq_len(max_start)) {
      band_list[[row_id]] <- data.frame(
        band_row_id  = row_id,
        family_id    = fam_id,
        band_k       = k,
        start_idx    = s,
        end_idx      = s + k - 1,
        period_start = periods[s],
        period_end   = periods[s + k - 1],
        band_sum     = sums[s],
        band_mean    = means[s]
      )
      row_id <- row_id + 1
    }
  }
}

candidate_bands <- bind_rows(band_list)
write.csv(candidate_bands, "ALL_CANDIDATE_BANDS_MEAN_SUM.csv", row.names = FALSE)

message("Stage 1 complete: saved ALL_CANDIDATE_BANDS_MEAN_SUM.csv")

################################################################################
#################### STAGE 2: Within-dyad ranking within each k + extract rank-1

message("Stage 2: computing within-dyad within-k ranks (mean & sum) and extracting rank-1 bands...")

candidate_bands_ranked <- candidate_bands %>%
  group_by(family_id, band_k) %>%
  mutate(
    mean_rank_k = percent_rank(band_mean),
    sum_rank_k  = percent_rank(band_sum)
  ) %>%
  ungroup()

best_by_mean <- candidate_bands_ranked %>%
  group_by(family_id, band_k) %>%
  arrange(desc(band_mean), desc(band_sum), start_idx) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(best_type = "best_by_mean")

best_by_sum <- candidate_bands_ranked %>%
  group_by(family_id, band_k) %>%
  arrange(desc(band_sum), desc(band_mean), start_idx) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(best_type = "best_by_sum")

best_bands <- bind_rows(best_by_mean, best_by_sum) %>%
  arrange(family_id, band_k, best_type)

write.csv(best_bands, "BEST_BANDS_BY_K_PER_FAMILY.csv", row.names = FALSE)

message("Stage 2 complete: saved BEST_BANDS_BY_K_PER_FAMILY.csv")

################################################################################
#################### STAGE 3: PERSON-SPECIFIC scan null distributions PER k

message("Stage 3: building PERSON-specific scan null distributions per dyad and k...")

person_scan_null <- vector("list", length = n_fam)
names(person_scan_null) <- as.character(df$family_id)

# Threshold table (one row per fam_id x k)
person_null_table <- expand.grid(
  fam_id = df$family_id,
  band_k = band_sizes,
  stringsAsFactors = FALSE
) %>%
  mutate(
    alpha = alpha_scan,
    mean_thr_95 = NA_real_,
    sum_thr_95  = NA_real_
  )

for (i in seq_len(n_fam)) {
  
  fam_id <- df$family_id[i]
  x_obs  <- as.numeric(R[i, ])
  
  message("  Building person scan null for fam_id = ", fam_id)
  
  fam_list <- vector("list", length(band_sizes))
  names(fam_list) <- as.character(band_sizes)
  
  for (k in band_sizes) {
    
    null_max_mean <- numeric(n_mc_scan)
    null_max_sum  <- numeric(n_mc_scan)
    
    for (mc in seq_len(n_mc_scan)) {
      x_perm <- sample(x_obs, size = n_per, replace = FALSE)
      ws <- window_sums_means(x_perm, k)
      
      null_max_mean[mc] <- max(ws$means, na.rm = TRUE)
      null_max_sum[mc]  <- max(ws$sums,  na.rm = TRUE)
    }
    
    null_max_mean_sorted <- sort(null_max_mean)
    null_max_sum_sorted  <- sort(null_max_sum)
    
    fam_list[[as.character(k)]] <- list(
      max_mean_sorted = null_max_mean_sorted,
      max_sum_sorted  = null_max_sum_sorted
    )
    
    # thresholds for this fam_id x k
    idx_row <- person_null_table$fam_id == fam_id & person_null_table$band_k == k
    person_null_table$mean_thr_95[idx_row] <- quantile(null_max_mean_sorted, 1 - alpha_scan, names = FALSE)
    person_null_table$sum_thr_95[idx_row]  <- quantile(null_max_sum_sorted,  1 - alpha_scan, names = FALSE)
  }
  
  person_scan_null[[as.character(fam_id)]] <- fam_list
}

write.csv(person_null_table, "PERSON_SCAN_NULL_THRESHOLDS_BY_K.csv", row.names = FALSE)
saveRDS(
  list(
    n_mc_scan = n_mc_scan,
    band_sizes = band_sizes,
    alpha_scan = alpha_scan,
    person_scan_null = person_scan_null
  ),
  "PERSON_SCAN_NULL_DISTS_BY_K.rds"
)

message("Stage 3 complete: saved PERSON_SCAN_NULL_THRESHOLDS_BY_K.csv + PERSON_SCAN_NULL_DISTS_BY_K.rds")

################################################################################
#################### STAGE 4: Compute PERSON scan-null p-values for rank-1 bands

message("Stage 4: computing PERSON scan-null p-values for rank-1 candidates (per k)...")

person_obj <- readRDS("PERSON_SCAN_NULL_DISTS_BY_K.rds")
person_scan_null <- person_obj$person_scan_null

best_bands <- read.csv("BEST_BANDS_BY_K_PER_FAMILY.csv")

best_bands <- best_bands %>%
  rowwise() %>%
  mutate(
    scan_p_mean = {
      fam_chr <- as.character(family_id)
      k_chr   <- as.character(band_k)
      null_sorted <- person_scan_null[[fam_chr]][[k_chr]]$max_mean_sorted
      emp_p_upper(band_mean, null_sorted)
    },
    scan_p_sum = {
      fam_chr <- as.character(family_id)
      k_chr   <- as.character(band_k)
      null_sorted <- person_scan_null[[fam_chr]][[k_chr]]$max_sum_sorted
      emp_p_upper(band_sum, null_sorted)
    },
    sig_scan_mean = (scan_p_mean <= alpha_scan),
    sig_scan_sum  = (scan_p_sum  <= alpha_scan)
  ) %>%
  ungroup()

write.csv(best_bands, "BEST_BANDS_WITH_PERSON_SCAN_P.csv", row.names = FALSE)

message("Stage 4 complete: saved BEST_BANDS_WITH_PERSON_SCAN_P.csv")

################################################################################
#################### STAGE 5 (UPDATED): Keep primary band + all sig/non-overlapping additional bands

message("Stage 5: selecting primary band + all significant non-overlapping additional bands (MEAN-based; person null)...")

best_mean_only <- best_bands %>%
  filter(best_type == "best_by_mean")

optimized_select <- function(df_fam) {
  
  df_sorted <- df_fam %>%
    arrange(
      scan_p_mean,
      desc(band_mean),
      start_idx,
      end_idx
    )
  
  
  if (nrow(df_sorted) == 0) return(df_sorted)
  
  kept_rows <- logical(nrow(df_sorted))
  kept_intervals <- list()
  
  # (1) Always keep the best band (row 1) — even if not significant
  kept_rows[1] <- TRUE
  kept_intervals[[1]] <- c(df_sorted$start_idx[1], df_sorted$end_idx[1])
  
  # (2) Additional bands: must be significant AND non-overlapping
  if (nrow(df_sorted) >= 2) {
    for (r in 2:nrow(df_sorted)) {
      
      if (sum(kept_rows) >= max_bands_per_family) break
      
      # must be significant to be eligible after the first
      if (is.na(df_sorted$scan_p_mean[r]) || df_sorted$scan_p_mean[r] > alpha_scan) next
      
      s <- df_sorted$start_idx[r]
      e <- df_sorted$end_idx[r]
      
      if (!overlaps_any(s, e, kept_intervals)) {
        kept_rows[r] <- TRUE
        kept_intervals[[length(kept_intervals) + 1]] <- c(s, e)
      }
    }
  }
  
  df_sorted[kept_rows, , drop = FALSE] %>%
    mutate(
      band_id = row_number(),  # rank among kept bands
      sig_scan_mean = (scan_p_mean <= alpha_scan),
      sig_scan_sum  = (scan_p_sum  <= alpha_scan)
    )
}

optimized_bands <- best_mean_only %>%
  group_by(family_id) %>%
  group_modify(~ optimized_select(.x)) %>%
  ungroup()

optimized_bands_out <- optimized_bands %>%
  transmute(
    fam_id        = family_id,
    band_id       = band_id,
    band_k        = band_k,
    band_size     = band_k,
    start_idx     = start_idx,
    end_idx       = end_idx,
    period_start  = period_start,
    period_end    = period_end,
    
    band_mean     = band_mean,
    band_sum      = band_sum,
    
    mean_rank_k   = mean_rank_k,
    sum_rank_k    = sum_rank_k,
    
    scan_p_mean   = scan_p_mean,
    scan_p_sum    = scan_p_sum,
    sig_scan_mean = sig_scan_mean,
    sig_scan_sum  = sig_scan_sum,
    alpha_scan    = alpha_scan
  ) %>%
  arrange(fam_id, band_id)

write.csv(optimized_bands_out, "FREQUENCY_OPTIMIZED_BANDS_PERSON_NULL.csv", row.names = FALSE)

message("Stage 5 complete: saved FREQUENCY_OPTIMIZED_BANDS_PERSON_NULL.csv")

################################################################################
#################### STAGE 6: VISUALIZATION OUTPUT (one PDF per family)

message("Stage 6: generating per-family band plots (person null)...")

if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

# reload for safety
df_plot <- read.csv("SIG_RIDGE_BY_PERIOD.csv")
rig_cols_plot <- grep("^RIG_", names(df_plot), value = TRUE)
periods_plot <- as.numeric(str_remove(rig_cols_plot, "^RIG_"))
R_plot <- as.matrix(df_plot[, rig_cols_plot])

optimized_bands_out <- read.csv("FREQUENCY_OPTIMIZED_BANDS_PERSON_NULL.csv")

fam_to_row <- setNames(seq_len(nrow(df_plot)), df_plot$family_id)
fam_ids <- unique(optimized_bands_out$fam_id)

for (fam in fam_ids) {
  
  if (!(as.character(fam) %in% names(fam_to_row))) next
  i <- fam_to_row[[as.character(fam)]]
  
  ridge_i <- as.numeric(R_plot[i, ])
  
  fam_bands <- optimized_bands_out %>%
    filter(fam_id == fam) %>%
    arrange(band_id) %>%
    slice_head(n = top_n_bands_to_plot)
  
  # --- zero-pad family ID to 4 digits ---
  fam_pad <- sprintf("%04d", as.integer(as.character(fam)))
  
  out_file <- file.path(plot_dir, paste0("FAM_", fam_pad, "_bands.pdf"))
  pdf(out_file, width = plot_width, height = plot_height)
  
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)
  par(mar = c(5, 4, 4, 3))
  
  xlim <- range(ridge_i, finite = TRUE)
  ylim <- range(periods_plot, finite = TRUE)
  
  plot(
    x = ridge_i, y = periods_plot,
    type = "l",
    xlab = "Ridge magnitude",
    ylab = "Period (low = fast, high = slow)",
    main = paste0("Dyad ", fam_pad, " Optimized Band(s)"),
    xlim = xlim, ylim = ylim,
    lwd = 1,
    yaxt = "n"
  )
  
  y_ticks_in <- y_ticks[y_ticks >= ylim[1] & y_ticks <= ylim[2]]
  axis(2, at = y_ticks_in, labels = y_ticks_in, las = 1)
  
  usr <- par("usr")
  x_left <- usr[1]; x_right <- usr[2]
  
  if (nrow(fam_bands) > 0) {
    
    for (r in seq_len(nrow(fam_bands))) {
      y_bottom <- fam_bands$period_start[r]
      y_top    <- fam_bands$period_end[r]
      
      col_fill <- blue_palette[min(fam_bands$band_id[r], top_n_bands_to_plot)]
      
      rect(
        xleft = x_left,
        xright = x_right,
        ybottom = y_bottom,
        ytop = y_top,
        col = col_fill,
        border = NA
      )
    }
    
    # redraw ridge line on top
    lines(ridge_i, periods_plot, lwd = 1.5)
    
    # star always indicates the PRIMARY band (band_id == 1)
    dom <- fam_bands %>% slice(1)
    y_mid <- (dom$period_start + dom$period_end) / 2
    
    par(xpd = NA)
    x_star <- x_right + 0.03 * (x_right - x_left)
    
    text(
      x = x_star, y = y_mid,
      labels = "*",
      cex = 1.6,
      font = 2
    )
    par(xpd = FALSE)
    
  } else {
    lines(ridge_i, periods_plot, lwd = 1.5)
  }
  
  dev.off()
}

message("Stage 6 complete: plots written to ", plot_dir)
message("Done.")


################################################################################
################# FINAL EXPORT: keep exactly 1 band (FOB) per dyad
# Rule:
#   - Choose the band with the smallest scan_p_mean.
#   - Tie-breakers (deterministic):
#       1) larger band_mean
#       2) earlier start_idx
#   - This returns 1 row per fam_id no matter what.

FOBS <- read.csv("FREQUENCY_OPTIMIZED_BANDS_PERSON_NULL.csv")

FOBS_FINAL <- FOBS %>%
  group_by(fam_id) %>%
  arrange(
    scan_p_mean,
    desc(band_mean),
    start_idx
  ) %>%
  slice(1) %>%
  ungroup() %>%
  transmute(
    fam_id,
    band_k,
    period_start,
    period_end,
    band_mean,
    scan_p_mean,
    sig_scan_mean
  )

# Convert period (seconds) to Hz
FOBS_FINAL <- FOBS_FINAL %>%
  mutate(
    hz_start = 1 / period_start,
    hz_end   = 1 / period_end
  )

# Add sample label (kept from your logic)
FOBS_FINAL$sample <- ifelse(FOBS_FINAL$fam_id > 105, "DORRY", "PDM")

# Save final (1 band per dyad)
write.csv(
  FOBS_FINAL,
  "FINAL_FREQUENCY_OPTIMIZED_BAND_SELECTION.csv",
  row.names = FALSE
)

# OPTIONAL QC: confirm 1 row per dyad
stopifnot(nrow(FOBS_FINAL) == dplyr::n_distinct(FOBS_FINAL$fam_id))



################################################################################
################# PLOT FOB DISTRIBUTIONS (1 band per dyad)

# ---- SIZE & FONT CONTROLS ----
axis_title_size  <- 25
axis_text_size   <- 20
tick_length      <- 0.35   # cm
legend_text_size <- 15
legend_title_size <- 15
font_family <- "Arial"

df_plot2 <- FOBS_FINAL %>%
  mutate(
    hz_low  = pmin(hz_start, hz_end),
    hz_high = pmax(hz_start, hz_end),
    hz_mid  = (hz_low + hz_high) / 2
  ) %>%
  arrange(desc(hz_mid)) %>%              # FASTEST first
  mutate(dyad_rank = row_number())

y_top <- max(df_plot2$dyad_rank)

p <- ggplot(df_plot2) +
  
  # ---- Background frequency regions ----
annotate(
  "rect",
  xmin = .0392, xmax = 0.15, ymin = -Inf, ymax = Inf,
  fill = "lightgreen", alpha = 0.18
) +
  annotate(
    "rect",
    xmin = 0.15, xmax = 0.40, ymin = -Inf, ymax = Inf,
    fill = "lightblue", alpha = 0.18
  ) +
  
  # ---- Region labels ----
annotate(
  "text",
  x = .142, y = y_top,
  label = "Low Frequency",
  angle = 90, hjust = 1, vjust = 1,
  size = 6, family = font_family
) +
  annotate(
    "text",
    x = .391, y = y_top,
    label = "High Frequency",
    angle = 90, hjust = 1, vjust = 1,
    size = 6, family = font_family
  ) +
  
  # ---- Dyad frequency bands ----
geom_rect(
  aes(
    xmin = hz_low,
    xmax = hz_high,
    ymin = dyad_rank - 0.35,
    ymax = dyad_rank + 0.35,
    fill = sample
  ),
  alpha = 0.85,
  color = NA
) +
  
  # ---- Axes ----
scale_x_continuous(
  limits = c(.030, 0.40),
  breaks = c(.040, .075, .15, .275, .40)
) +
  scale_y_continuous(
    breaks = NULL,
    expand = expansion(mult = c(0.01, 0.03))
  ) +
  
  labs(
    x = "Frequency (Hz)",
    y = "Dyad",
    fill = "Sample"
  ) +
  
  # ---- Theme ----
theme_minimal(base_family = font_family) +
  theme(
    axis.title.x = element_text(size = axis_title_size),
    axis.title.y = element_text(size = axis_title_size, margin = margin(r = -51)),
    axis.text.x  = element_text(size = axis_text_size),
    axis.text.y  = element_blank(),  # keep axis, hide labels
    axis.ticks.x = element_line(),
    axis.ticks.length = unit(tick_length, "cm"),
    
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.text  = element_text(size = legend_text_size),
    legend.title = element_text(size = legend_title_size)
  )

ggsave(
  filename = "FOB_DISTRIBUTION.png",
  plot = p,
  width = 10,
  height = 8,
  units = "in",
  dpi = 300,
  bg = "white"
)



