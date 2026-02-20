################################################################################
# SAMPLE-LEVEL PHASE ANALYSIS (using nlme instead of lme4)
# - Loads DYAD_PHASE_COSINE_SUMMARY.csv
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

DYAD_PHS <- read.csv("DYAD_PHASE_COSINE_SUMMARY.csv",
                     stringsAsFactors = FALSE)

DYAD_PHS <- DYAD_PHS %>%
  dplyr::filter(converged == TRUE) %>%
  dplyr::filter(!is.na(PHS_MU), !is.na(PHS_A), !is.na(PHS_B), !is.na(PHS_C))

DYAD_PHS$BAND <- factor(toupper(DYAD_PHS$BAND), levels = c("FOB", "RSB", "HF", "LF"))
DYAD_PHS$ID <- as.factor(DYAD_PHS$ID)

message("Loaded ", nrow(DYAD_PHS), " converged coherence fits")

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

###############################################################################
# Build long plot data (REMOVE c)

plot_data <- DYAD_PHS %>%
  dplyr::select(ID, BAND, PHS_MU, PHS_A, PHS_B) %>%
  tidyr::pivot_longer(
    cols = c(PHS_MU, PHS_A, PHS_B),
    names_to = "parameter",
    values_to = "value"
  ) %>%
  dplyr::mutate(
    parameter = dplyr::case_when(
      parameter == "PHS_MU" ~ "Average (μ) [radians]",
      parameter == "PHS_A"  ~ "Amplitude (a) [radians]",
      parameter == "PHS_B"  ~ "Period (b) [seconds]"
    ),
    parameter = factor(
      parameter,
      levels = c("Average (μ) [radians]", "Amplitude (a) [radians]", "Period (b) [seconds]")
    )
  )

###############################################################################
# Circular helpers for μ (radians)
# - Circular mean: atan2(mean(sin), mean(cos))
# - 95% CI: bootstrap, then compute percentile CI on circular differences
#   (avoids boundary issues at -pi/pi)

circ_mean <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  atan2(mean(sin(x)), mean(cos(x)))
}

circ_ci_boot <- function(x, B = 2000, alpha = 0.05, seed = 1) {
  x <- x[is.finite(x)]
  n <- length(x)
  if (n < 2) return(c(NA_real_, NA_real_))
  
  set.seed(seed)
  m_hat <- circ_mean(x)
  
  boot_m <- replicate(B, {
    xb <- sample(x, size = n, replace = TRUE)
    circ_mean(xb)
  })
  
  # Circular difference around m_hat in (-pi, pi]
  d <- atan2(sin(boot_m - m_hat), cos(boot_m - m_hat))
  
  lo_d <- stats::quantile(d, probs = alpha/2, na.rm = TRUE, names = FALSE)
  hi_d <- stats::quantile(d, probs = 1 - alpha/2, na.rm = TRUE, names = FALSE)
  
  # Put CI back around m_hat
  ci_low  <- m_hat + lo_d
  ci_high <- m_hat + hi_d
  
  c(ci_low, ci_high)
}

###############################################################################
# Summary stats:
# - For μ: circular mean + bootstrap CI
# - For a, b: arithmetic mean + t-based 95% CI

summary_stats <- plot_data %>%
  dplyr::group_by(parameter, BAND) %>%
  dplyr::summarise(
    n = sum(!is.na(value)),
    mean_val = dplyr::if_else(
      first(parameter) == "Average (μ) [radians]",
      circ_mean(value),
      mean(value, na.rm = TRUE)
    ),
    ci_low = {
      if (first(parameter) == "Average (μ) [radians]") {
        circ_ci_boot(value, B = 2000, alpha = 0.05, seed = 1)[1]
      } else {
        sd_val <- sd(value, na.rm = TRUE)
        se_val <- sd_val / sqrt(sum(!is.na(value)))
        tcrit  <- stats::qt(0.975, df = sum(!is.na(value)) - 1)
        mean(value, na.rm = TRUE) - tcrit * se_val
      }
    },
    ci_high = {
      if (first(parameter) == "Average (μ) [radians]") {
        circ_ci_boot(value, B = 2000, alpha = 0.05, seed = 1)[2]
      } else {
        sd_val <- sd(value, na.rm = TRUE)
        se_val <- sd_val / sqrt(sum(!is.na(value)))
        tcrit  <- stats::qt(0.975, df = sum(!is.na(value)) - 1)
        mean(value, na.rm = TRUE) + tcrit * se_val
      }
    },
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    label_txt = sprintf("M = %.2f\n(%.2f, %.2f)", mean_val, ci_low, ci_high)
  )


###############################################################################
# Plot (1 x 3 facet)

###############################################################################
# Plot (1 x 3 facet) — MAIN (NO LEGEND)

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
    y = "Phase Value",
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
    legend.position = "none"   # <- IMPORTANT: kill legend here
  )

###############################################################################
# STATS STRIP (NO LEGEND)

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
    legend.position = "none"  # <- IMPORTANT: kill legend here too
  )

###############################################################################
# DUMMY LEGEND PLOT (THE ONLY LEGEND)

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

# Extract legend as a patchwork element
leg <- wrap_elements(get_legend(p_leg))

###############################################################################
# Combine: main -> stats -> legend (legend truly last)

p_heterogeneity <- (p_main / p_stats / leg) +
  plot_layout(heights = c(1, 0.22, 0.10))

ggsave(
  "PHASE_PARAMETER_BANDWISE_HETEROGENEITY.png",
  p_heterogeneity, width = 14, height = 7.2, dpi = 300, bg = "white"
)

message("DONE: saved PHASE_PARAMETER_BANDWISE_HETEROGENEITY.png")

