################################################################################
# SAMPLE-LEVEL COHERENCE ANALYSIS (using nlme instead of lme4)
# - Loads DYAD_COHERENCE_COSINE_SUMMARY.csv
# - Visualizes heterogeneity in parameters (mu, a, b, c) across bands
# - Computes random effects models (3 types per parameter) using nlme
# - SAVES random effects results + model comparisons to CSV
# - Creates sample-level plot with reconstructed curves
################################################################################

library(dplyr)
library(tidyr)
library(ggplot2)
library(nlme)
library(stringr)
library(patchwork)
library(cowplot)
###############################################################################
# 0) PATHS + LOAD

work_dir <- "/home/cjh37695/WAVELET_WORKING/TEST_RESP"
setwd(work_dir)

DYAD_COH <- read.csv("DYAD_COHERENCE_COSINE_SUMMARY.csv",
                     stringsAsFactors = FALSE)

DYAD_COH <- DYAD_COH %>%
  dplyr::filter(converged == TRUE) %>%
  dplyr::filter(!is.na(COH_MU), !is.na(COH_A), !is.na(COH_B), !is.na(COH_C))

DYAD_COH$BAND <- factor(toupper(DYAD_COH$BAND), levels = c("FOB", "RSB", "HF", "LF"))
DYAD_COH$ID <- as.factor(DYAD_COH$ID)

message("Loaded ", nrow(DYAD_COH), " converged coherence fits")

###############################################################################
# 1) HETEROGENEITY VISUALIZATION: 2x2 FACET PLOT (MEAN + 95% CI)

# -----------------------------
# Easy knobs to tweak text sizes
# -----------------------------
BASE_FONT_SIZE  <- 20   # overall plot text
TITLE_SIZE      <- 24
STRIP_SIZE      <- 18
AXIS_TITLE_SIZE <- 20
AXIS_TEXT_SIZE  <- 20
LABEL_SIZE      <- 3.8  # mean & 95% text size

# -----------------------------
# Band color mapping (edit here)
# Make sure names match your DYAD_COH$BAND values exactly
# -----------------------------
band_palette <- c(
  "FOB" = "#D55E00",  # orange
  "HF"  = "#6A3D9A",  # purple
  "LF"  = "#0072B2",  # blue
  "RSB" = "#009E73"   # green
)

plot_data <- DYAD_COH %>%
  dplyr::select(ID, BAND, COH_MU, COH_A, COH_B) %>%   # <- dropped COH_C
  tidyr::pivot_longer(
    cols = c(COH_MU, COH_A, COH_B),
    names_to = "parameter",
    values_to = "value"
  ) %>%
  dplyr::mutate(
    parameter = dplyr::case_when(
      parameter == "COH_MU" ~ "Average (μ) [coherence]",
      parameter == "COH_A"  ~ "Amplitude (a) [coherence]",
      parameter == "COH_B"  ~ "Period (b) [seconds]"
    ),
    parameter = factor(
      parameter,
      levels = c("Average (μ) [coherence]",
                 "Amplitude (a) [coherence]",
                 "Period (b) [seconds]")
    )
  )

###############################################################################
# Compute mean + 95% CI per parameter x band (with bottom-of-panel 2-line labels)

summary_stats <- plot_data %>%
  dplyr::group_by(parameter, BAND) %>%
  dplyr::summarise(
    n        = sum(!is.na(value)),
    mean_val = mean(value, na.rm = TRUE),
    sd_val   = sd(value, na.rm = TRUE),
    se_val   = sd_val / sqrt(n),
    tcrit    = stats::qt(0.975, df = n - 1),
    ci_low   = mean_val - tcrit * se_val,
    ci_high  = mean_val + tcrit * se_val,
    .groups  = "drop"
  ) %>%
  dplyr::group_by(parameter) %>%
  dplyr::mutate(
    # facet-specific range (so this works with scales = "free_y")
    y_min = min(plot_data$value[plot_data$parameter == first(parameter)], na.rm = TRUE),
    y_max = max(plot_data$value[plot_data$parameter == first(parameter)], na.rm = TRUE),
    y_rng = y_max - y_min,
    
    # park label in the lower margin area (below data), but still visible via expand()
    label_y = y_min - 0.12 * y_rng,
    
    # 2-line label (skinny): "M = <mean>" then "(LB, UB)"
    label_txt = sprintf("M = %.2f\n(%.2f, %.2f)", mean_val, ci_low, ci_high)
  ) %>%
  dplyr::ungroup()

###############################################################################
# Plot (band-colored means + 95% CI + bottom labels)

p_main <- ggplot(plot_data, aes(x = BAND, y = value)) +
  geom_rect(
    data = data.frame(
      parameter = rep(levels(plot_data$parameter), each = 2),
      xmin = rep(c(0.5, 2.5), times = 3),
      xmax = rep(c(1.5, 3.5), times = 3),
      ymin = -Inf, ymax = Inf
    ),
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = "grey90", alpha = 0.3, inherit.aes = FALSE
  ) +
  geom_jitter(width = 0.2, height = 0, alpha = 0.35, size = 0.8, color = "grey35") +
  geom_errorbar(
    data = summary_stats,
    aes(x = BAND, ymin = ci_low, ymax = ci_high, color = BAND),
    width = 0.25, linewidth = 1.1,
    inherit.aes = FALSE
  ) +
  geom_point(
    data = summary_stats,
    aes(x = BAND, y = mean_val, color = BAND),
    size = 3,
    inherit.aes = FALSE
  ) +
  facet_wrap(~ parameter, ncol = 3, scales = "free_y") +
  scale_color_manual(values = band_palette, drop = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.06))) +
  labs(
    x = NULL,
    y = "Coherence Value",
    color = "Band"
  ) +
  theme_minimal(base_family = "serif", base_size = BASE_FONT_SIZE) +
  theme(
    strip.text   = element_text(face = "bold", size = STRIP_SIZE),
    axis.title.y = element_text(size = AXIS_TITLE_SIZE),
    axis.text    = element_text(size = AXIS_TEXT_SIZE),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    plot.margin = margin(10, 10, 2, 10),
    legend.position = "none"
  )

###############################################################################
# STATS STRIP (NO LEGEND) — BELOW PLOT, ABOVE LEGEND

p_stats <- ggplot(summary_stats, aes(x = BAND, y = 1, label = label_txt, color = BAND)) +
  geom_text(size = LABEL_SIZE, family = "serif", lineheight = 0.95, vjust = 1) +
  facet_wrap(~ parameter, ncol = 3) +
  scale_color_manual(values = band_palette, drop = FALSE) +
  scale_y_continuous(limits = c(0.5, 1.1), expand = c(0, 0)) +
  labs(x = NULL, y = NULL) +
  theme_void(base_family = "serif", base_size = BASE_FONT_SIZE) +
  theme(
    strip.text = element_blank(),
    plot.margin = margin(0, 10, 2, 10),
    legend.position = "none"
  )

###############################################################################
# DUMMY LEGEND PLOT (ONLY LEGEND)

p_leg <- ggplot(summary_stats, aes(x = BAND, y = mean_val, color = BAND)) +
  geom_point(size = 3) +
  scale_color_manual(values = band_palette, drop = FALSE) +
  labs(color = "Band") +
  theme_void(base_family = "serif", base_size = BASE_FONT_SIZE) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = AXIS_TEXT_SIZE),
    legend.text  = element_text(size = AXIS_TEXT_SIZE)
  )

leg <- patchwork::wrap_elements(cowplot::get_legend(p_leg))

###############################################################################
# Combine: main -> stats -> legend

p_heterogeneity <- (p_main / p_stats / leg) +
  plot_layout(heights = c(1, 0.22, 0.10))

ggsave(
  "COHERENCE_PARAMETER_BANDWISE_HETEROGENEITY.png",
  p_heterogeneity, width = 14, height = 7.2, dpi = 300, bg = "white"
)

message("DONE: saved COHERENCE_PARAMETER_BANDWISE_HETEROGENEITY.png")


###############################################################################
# 2) SAMPLE-LEVEL PLOT: RECONSTRUCTED CURVES

coh_file <- "SAMPLE_TRUNC_COHERENCE_DS_SERIES.csv"
COH_WIDE <- read.csv(coh_file, check.names = FALSE, stringsAsFactors = FALSE)
COH_WIDE <- COH_WIDE[, -1]

time_cols <- names(COH_WIDE)[stringr::str_detect(names(COH_WIDE), "^t[_\\.]")]
time_vals <- suppressWarnings(as.numeric(sub("^t[_\\.]", "", time_cols)))
time_vals <- time_vals[!is.na(time_vals)]
t_grid <- sort(unique(time_vals))

y_hat_dyad <- function(t, mu, a, b, c) {
  mu + a * cos((2*pi/b) * (t - c))
}

bands <- c("FOB", "RSB", "HF", "LF")
coh_plot_list <- list()

for (band in bands) {
  
  message("Reconstructing trajectories for band: ", band)
  
  dyad_band <- DYAD_COH %>%
    dplyr::filter(toupper(BAND) == toupper(band))
  
  if (nrow(dyad_band) == 0) next
  
  dyad_trajs <- lapply(seq_len(nrow(dyad_band)), function(i) {
    mu <- dyad_band$COH_MU[i]
    a  <- dyad_band$COH_A[i]
    b  <- dyad_band$COH_B[i]
    c  <- dyad_band$COH_C[i]
    id <- dyad_band$ID[i]
    
    y_fitted <- y_hat_dyad(t_grid, mu, a, b, c)
    y_fitted <- pmax(0, pmin(1, y_fitted))
    
    data.frame(
      ID = id,
      BAND = band,
      t = t_grid,
      y_fitted = y_fitted,
      stringsAsFactors = FALSE
    )
  })
  
  dyad_traj_df <- dplyr::bind_rows(dyad_trajs)
  
  sample_curve <- dyad_traj_df %>%
    dplyr::group_by(BAND, t) %>%
    dplyr::summarise(
      y_sample = mean(y_fitted, na.rm = TRUE),
      .groups = "drop"
    )
  
  coh_plot_list[[band]] <- list(
    dyad_trajs = dyad_traj_df,
    sample_curve = sample_curve
  )
}

all_dyad_trajs <- dplyr::bind_rows(lapply(coh_plot_list, `[[`, "dyad_trajs"))
all_sample_curves <- dplyr::bind_rows(lapply(coh_plot_list, `[[`, "sample_curve"))

p_sample <- ggplot() +
  geom_line(
    data = all_dyad_trajs,
    aes(x = t, y = y_fitted, group = ID),
    linewidth = 0.35,
    alpha = 0.25,
    color = "grey50",
    na.rm = TRUE
  ) +
  geom_line(
    data = all_sample_curves,
    aes(x = t, y = y_sample),
    linewidth = 1.4,
    color = "#D55E00",
    na.rm = TRUE
  ) +
  facet_wrap(~ BAND, ncol = 2) +
  labs(
   # title = "Coherence: sample-level mean from dyad-level cosine fits",
    x = "Time (seconds)",
    y = "Coherence"
  ) +
  ylim(0, 1) +
  theme_minimal(base_family = "Arial") +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "none"
  )

ggsave("COHERENCE_SAMPLE_FROM_DYADS_FACET.png",
       p_sample, width = 12, height = 7, dpi = 300, bg = "white")

message("DONE: saved COHERENCE_SAMPLE_FROM_DYADS_FACET.png")
message("\n=== SAMPLE-LEVEL COHERENCE ANALYSIS COMPLETE ===\n")

