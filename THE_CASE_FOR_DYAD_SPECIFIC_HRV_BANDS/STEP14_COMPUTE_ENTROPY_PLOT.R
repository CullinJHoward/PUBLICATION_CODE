### LIBRARY
library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
library(RespirAnalyzer)
library(ggplot2)


################################################################################
##################### COMPUTE BAND ENTROPY 

## SET WORKING DIRECTORY 

work_dir <- "/home/cjh37695/WAVELET_WORKING/TEST_RESP"
setwd(work_dir)

FOB_SUM <- read.csv("FINAL_FREQUENCY_OPTIMIZED_BAND_SELECTION.csv")

###### LOAD IN TRUNCATED DATA 

PAIRED_COH <- read.csv("SAMPLE_TRUNC_COHERENCE_FULL_4Hz_SERIES.csv",
                       stringsAsFactors = FALSE)

NONPAIRED_COH <- read.csv("RESAMPLING_VALIDATION/CV_SAMPLE_TRUNC_COHERENCE_FULL_4Hz_SERIES.csv",
                          stringsAsFactors = FALSE)

PAIRED_PHS <- read.csv("SAMPLE_TRUNC_PHASE_FULL_4Hz_SERIES.csv",
                       stringsAsFactors = FALSE)

NONPAIRED_PHS <- read.csv("RESAMPLING_VALIDATION/CV_SAMPLE_TRUNC_PHASE_FULL_4Hz_SERIES.csv",
                          stringsAsFactors = FALSE)

# ADD PAIRED VS NONPAIRED ID

PAIRED_COH$PAIR <- "FAM"
NONPAIRED_COH$PAIR <- "RAND"
PAIRED_PHS$PAIR <- "FAM"
NONPAIRED_PHS$PAIR <- "RAND"


## FULL DF's 

COHERENCE <- rbind(PAIRED_COH, NONPAIRED_COH)

PHASE <- rbind(PAIRED_PHS, NONPAIRED_PHS)

#QC CHECK
table(PHASE$PAIR)


## REMOVE DYADS FAILING QC CHECKS
BAD_IDS <- c(1035)  # VISUAL INSPECTION

COHERENCE$FAIL_QC <- ifelse(COHERENCE$ID %in% BAD_IDS, 1, 0)
COHERENCE <- subset(COHERENCE, FAIL_QC == 0)
PHASE$FAIL_QC <- ifelse(PHASE$ID %in% BAD_IDS, 1, 0)
PHASE <- subset(PHASE, FAIL_QC == 0)

### REMOVE CLEANING vAR 
COHERENCE <- COHERENCE %>%
  select(-c("FAIL_QC", "period_start", "period_end", "HZ_start", "HZ_end"))

PHASE <- PHASE %>%
  select(-c("FAIL_QC", "period_start", "period_end", "HZ_start", "HZ_end"))


## LOCATE PAIR EARLIER THAN THE TIME SERIES 
COHERENCE <- COHERENCE %>% relocate(PAIR, .after = BAND)
PHASE <- PHASE %>% relocate(PAIR, .after = BAND)

###############################################################################
# PIVOT TO LONG FORMAT 


COHERENCE_LONG <- COHERENCE %>%
  pivot_longer(
    cols      = 4:ncol(COHERENCE),
    names_to  = "TIME",
    values_to = "COHERENCE"
  ) %>%
  arrange(ID, BAND, PAIR, TIME)

PHASE_LONG <- PHASE %>%
  pivot_longer(
    cols      = 4:ncol(PHASE),
    names_to  = "TIME",
    values_to = "PHASE"
  ) %>%
  arrange(ID, BAND, PAIR, TIME)

## CREATE A BETTER NESTED BAND/PAIR VARIABLE 

COHERENCE_LONG$BAND_PAIR <- ifelse(COHERENCE_LONG$BAND == "FOB" & COHERENCE_LONG$PAIR == "FAM",  "FOB_FAM",
                                   ifelse(COHERENCE_LONG$BAND == "FOB" & COHERENCE_LONG$PAIR == "RAND", "FOB_RAND",
                                          ifelse(COHERENCE_LONG$BAND == "RSB" & COHERENCE_LONG$PAIR == "FAM",  "RSB_FAM",
                                                 ifelse(COHERENCE_LONG$BAND == "RSB" & COHERENCE_LONG$PAIR == "RAND", "RSB_RAND",
                                                        ifelse(COHERENCE_LONG$BAND == "HF"  & COHERENCE_LONG$PAIR == "FAM",  "HF_FAM",
                                                               ifelse(COHERENCE_LONG$BAND == "HF"  & COHERENCE_LONG$PAIR == "RAND", "HF_RAND",
                                                                      ifelse(COHERENCE_LONG$BAND == "LF"  & COHERENCE_LONG$PAIR == "FAM",  "LF_FAM",
                                                                             ifelse(COHERENCE_LONG$BAND == "LF"  & COHERENCE_LONG$PAIR == "RAND", "LF_RAND", NA))))))))


PHASE_LONG$BAND_PAIR <- ifelse(PHASE_LONG$BAND == "FOB" & PHASE_LONG$PAIR == "FAM",  "FOB_FAM",
                               ifelse(PHASE_LONG$BAND == "FOB" & PHASE_LONG$PAIR == "RAND", "FOB_RAND",
                                      ifelse(PHASE_LONG$BAND == "RSB" & PHASE_LONG$PAIR == "FAM",  "RSB_FAM",
                                             ifelse(PHASE_LONG$BAND == "RSB" & PHASE_LONG$PAIR == "RAND", "RSB_RAND",
                                                    ifelse(PHASE_LONG$BAND == "HF"  & PHASE_LONG$PAIR == "FAM",  "HF_FAM",
                                                           ifelse(PHASE_LONG$BAND == "HF"  & PHASE_LONG$PAIR == "RAND", "HF_RAND",
                                                                  ifelse(PHASE_LONG$BAND == "LF"  & PHASE_LONG$PAIR == "FAM",  "LF_FAM",
                                                                         ifelse(PHASE_LONG$BAND == "LF"  & PHASE_LONG$PAIR == "RAND", "LF_RAND", NA))))))))


## UPDATE TIME VARIABLE 
COHERENCE_LONG <- COHERENCE_LONG %>%
  group_by(ID, BAND_PAIR) %>%            
  arrange(ID, BAND_PAIR, TIME, .by_group = TRUE) %>%  
  mutate(TIME = (row_number() - 1) * 0.25) %>%
  ungroup()

PHASE_LONG <- PHASE_LONG %>%
  group_by(ID, BAND_PAIR) %>%            
  arrange(ID, BAND_PAIR, TIME, .by_group = TRUE) %>%  
  mutate(TIME = (row_number() - 1) * 0.25) %>%
  ungroup()


###############################################################################
# REMOVE MISSING VALUES  

COHERENCE_LONG_CLEAN <- COHERENCE_LONG %>%
  dplyr::filter(!is.na(COHERENCE))

PHASE_LONG_CLEAN <- PHASE_LONG %>%
  dplyr::filter(!is.na(PHASE))

## CHECK MISSINGNESS (ALL MUST BE GONE!)
COH_TBL <- COHERENCE_LONG_CLEAN %>%
  group_by(ID, BAND_PAIR) %>%
  summarize(
    n_total = n(),
    n_nonNA = sum(!is.na(COHERENCE)),
    .groups = "drop"
  )

COH_TBL$FLAG <- ifelse(COH_TBL$n_total - COH_TBL$n_nonNA != 0, 1, 0)
table(COH_TBL$FLAG)

PHS_TBL <- PHASE_LONG_CLEAN %>%
  group_by(ID, BAND_PAIR) %>%
  summarize(
    n_total = n(),
    n_nonNA = sum(!is.na(PHASE)),
    .groups = "drop"
  )
PHS_TBL$FLAG <- ifelse(PHS_TBL$n_total - PHS_TBL$n_nonNA != 0, 1, 0)
table(PHS_TBL$FLAG)

#### ENSURE WE ONLY LOOK AT DYADS IN THE ORIGINAL PAIRING 
VALID_IDS <- setdiff(unique(FOB_SUM$ID), "1035")

COHERENCE_LONG_CLEAN_RED <- COHERENCE_LONG_CLEAN[COHERENCE_LONG_CLEAN$ID %in% VALID_IDS, ]
PHASE_LONG_CLEAN_RED <- PHASE_LONG_CLEAN[PHASE_LONG_CLEAN$ID %in% VALID_IDS, ]

################################################################################
########################## RUN MSE ENTROPY 

###### COHERENCE 
# ----------------------------
# Settings
# ----------------------------
TAU_VEC <- 1:20
M_VAL   <- 2
R_VAL   <- 0.15

df_in <- COHERENCE_LONG_CLEAN_RED

# ----------------------------
# Helper: run MSE and return a tibble(tau, mse)
# ----------------------------

run_mse_one <- function(x, tau = tau, m = M_VAL, r = R_VAL) {
  out <- RespirAnalyzer::MSE(x, tau = tau, m = m, r = r, I = length(x))
  
  # out is a data.frame: tau, m, rSD, SampEn
  out %>%
    as_tibble() %>%
    transmute(
      tau = tau,
      m = m,
      rSD = rSD,          # the actual tolerance used (r * SD_used)
      SampEn = SampEn     # entropy
    )
}

# ----------------------------
# Step A: difference within each ID × BAND_PAIR series
# ----------------------------
df_dx <- df_in %>%
  group_by(ID, BAND_PAIR) %>%
  arrange(TIME, .by_group = TRUE) %>%
  mutate(dx = COHERENCE - lag(COHERENCE)) %>%
  filter(!is.na(dx)) %>%
  ungroup()

# ----------------------------
# Step B: run MSE per series and build Option-1 summary
#   one row per ID × BAND_PAIR × tau
# ----------------------------

entropy_summary <- df_dx %>%
  group_by(ID, BAND_PAIR) %>%
  summarise(
    BAND = first(BAND),
    PAIR = first(PAIR),
    n_dx = n(),
    sd_dx = sd(dx, na.rm = TRUE),     # SD of differenced series (audit)
    mse_tbl = list(run_mse_one(dx, tau = TAU_VEC, m = M_VAL, r = R_VAL)),
    .groups = "drop"
  ) %>%
  unnest(mse_tbl) %>%
  relocate(ID, BAND_PAIR, BAND, PAIR, tau, SampEn, rSD, m, n_dx, sd_dx)

## ADD FLAG FOR ODD ENTROPY
entropy_summary <- entropy_summary %>%
  mutate(flag_zero = SampEn == 0)

entropy_summary %>%
  summarise(
    n_total = n(),
    n_zero = sum(flag_zero, na.rm = TRUE),
    pct_zero = mean(flag_zero, na.rm = TRUE) * 100
  )


###### PHASE 

################################################################################
########################## RUN MSE ENTROPY ON CIRCULAR PHASE (DIFFERENCED)

TAU_VEC <- 1:20
M_VAL   <- 2
R_VAL   <- 0.15

df_in <- PHASE_LONG_CLEAN_RED

# OPTIONAL: if PHASE is in DEGREES, uncomment this conversion:
# df_in <- df_in %>%
#   mutate(PHASE = PHASE * pi / 180)

# Helper: MSE returns data.frame with tau, m, rSD, SampEn (as you observed)
run_mse_one <- function(x, tau = TAU_VEC, m = M_VAL, r = R_VAL) {
  out <- RespirAnalyzer::MSE(x, tau = tau, m = m, r = r, I = length(x))
  out %>%
    as_tibble() %>%
    transmute(
      tau   = tau,
      m     = m,
      rSD   = rSD,
      SampEn = SampEn
    )
}

# Step A: circular first difference within each ID × BAND_PAIR series
df_dphi <- df_in %>%
  group_by(ID, BAND_PAIR) %>%
  arrange(TIME, .by_group = TRUE) %>%
  mutate(
    dphi = atan2(
      sin(PHASE - lag(PHASE)),
      cos(PHASE - lag(PHASE))
    )
  ) %>%
  filter(!is.na(dphi)) %>%
  ungroup()

# Step B: run MSE per series (one row per ID × BAND_PAIR × tau)
phase_entropy_summary <- df_dphi %>%
  group_by(ID, BAND_PAIR) %>%
  summarise(
    BAND = first(BAND),
    PAIR = first(PAIR),
    n_dphi = n(),
    sd_dphi = sd(dphi, na.rm = TRUE),
    mse_tbl = list(run_mse_one(dphi, tau = TAU_VEC, m = M_VAL, r = R_VAL)),
    .groups = "drop"
  ) %>%
  unnest(mse_tbl) %>%
  relocate(ID, BAND_PAIR, BAND, PAIR, tau, SampEn, rSD, m, n_dphi, sd_dphi)

# Flag zeros (same idea as before)
phase_entropy_summary <- phase_entropy_summary %>%
  mutate(flag_zero = SampEn == 0)

phase_entropy_summary %>%
  summarise(
    n_total = n(),
    n_zero = sum(flag_zero, na.rm = TRUE),
    pct_zero = mean(flag_zero, na.rm = TRUE) * 100
  )

################################################################################
######################## ENTROPY PLOT: COHERENCE vs PHASE (FAM only)

# - Side-by-side Coherence vs Phase within each BAND facet
# - Uses mean + 95% CI overlays (colored by BAND)
# - Observed points are colored by Modality (Coherence vs Phase) but NO modality legend
# - Only BAND legend is shown at bottom (dummy legend)
# - Removes x-axis label (no "Modality" label)
#
# Requires objects:
#   entropy_summary        (coherence MSE summary; has: ID, BAND, PAIR, tau, SampEn, ...)
#   phase_entropy_summary  (phase MSE summary; has: ID, BAND, PAIR, tau, SampEn, ...)
################################################################################

library(dplyr)
library(ggplot2)
library(tidyr)
library(patchwork)
library(cowplot)

###############################################################################
# Easy knobs
###############################################################################
BASE_FONT_SIZE  <- 20
TITLE_SIZE      <- 24
STRIP_SIZE      <- 18
AXIS_TITLE_SIZE <- 20
AXIS_TEXT_SIZE  <- 20
LABEL_SIZE      <- 6

###############################################################################
# Tau handling knobs (CHANGE HERE ONLY)
###############################################################################
MODE <- "mean_range"     # "single" or "mean_range"
TAU_CHOICE <- 4          # used if MODE == "single"
TAU_RANGE  <- 1:16       # used if MODE == "mean_range"

###############################################################################
# Title/labels
###############################################################################
OUTCOME_TITLE <- if (MODE == "single") {
  paste0("Sample Entropy (differenced) by Band: Coherence vs Phase (FAM only; tau = ", TAU_CHOICE, ")")
} else {
  paste0("Sample Entropy (differenced) by Band: Coherence vs Phase (FAM only; mean tau = ",
         min(TAU_RANGE), "–", max(TAU_RANGE), ")")
}
Y_LABEL <- "Averaged Entropy (SampEn)"

###############################################################################
# Palettes
###############################################################################
band_palette <- c(
  "FOB" = "#D55E00",
  "HF"  = "#6A3D9A",
  "LF"  = "#0072B2",
  "RSB" = "#009E73"
)

mod_palette <- c(
  "Coherence" = "black",
  "Phase"     = "grey55"
)

###############################################################################
# 1) Collapse tau within each ID x BAND for each modality (FAM only)
###############################################################################
collapse_mse <- function(df, modality_label) {
  if (MODE == "single") {
    df %>%
      filter(PAIR == "FAM", tau == TAU_CHOICE) %>%
      transmute(ID, BAND, Modality = modality_label, VALUE = SampEn)
  } else {
    df %>%
      filter(PAIR == "FAM", tau %in% TAU_RANGE) %>%
      group_by(ID, BAND) %>%
      summarise(VALUE = mean(SampEn, na.rm = TRUE), .groups = "drop") %>%
      mutate(Modality = modality_label) %>%
      select(ID, BAND, Modality, VALUE)
  }
}

plot_data <- bind_rows(
  collapse_mse(entropy_summary, "Coherence"),
  collapse_mse(phase_entropy_summary, "Phase")
) %>%
  mutate(
    BAND = factor(BAND, levels = c("FOB", "RSB", "HF", "LF")),
    Modality = factor(Modality, levels = c("Coherence", "Phase"))
  )

###############################################################################
# 2) Summary stats per BAND x Modality (mean + 95% CI)
###############################################################################
summary_stats <- plot_data %>%
  group_by(BAND, Modality) %>%
  summarise(
    n        = sum(!is.na(VALUE)),
    mean_val = mean(VALUE, na.rm = TRUE),
    sd_val   = sd(VALUE, na.rm = TRUE),
    se_val   = sd_val / sqrt(n),
    tcrit    = ifelse(n > 1, qt(0.975, df = n - 1), NA_real_),
    ci_low   = mean_val - tcrit * se_val,
    ci_high  = mean_val + tcrit * se_val,
    .groups  = "drop"
  ) %>%
  complete(BAND, Modality, fill = list(
    n = 0, mean_val = NA_real_, sd_val = NA_real_,
    se_val = NA_real_, tcrit = NA_real_,
    ci_low = NA_real_, ci_high = NA_real_
  )) %>%
  mutate(
    label_txt = ifelse(
      n > 0,
      sprintf("M = %.2f\n(%.2f, %.2f)", mean_val, ci_low, ci_high),
      ""
    )
  )

###############################################################################
# 3) MAIN PLOT: jittered points + band-colored mean/CI
#    - points colored by Modality but NO modality legend
#    - band mean/CI colored by BAND
#    - combined scale contains both palettes, but legend is suppressed here
###############################################################################
p_main <- ggplot(plot_data, aes(x = Modality, y = VALUE)) +
  geom_point(
    aes(color = Modality),
    position = position_jitter(width = 0.18, height = 0),
    size = 2.0,
    show.legend = FALSE
  ) +
  geom_errorbar(
    data = summary_stats %>% filter(n > 0),
    aes(x = Modality, ymin = ci_low, ymax = ci_high, color = BAND),
    width = 0.20, linewidth = 1.1,
    inherit.aes = FALSE
  ) +
  geom_point(
    data = summary_stats %>% filter(n > 0),
    aes(x = Modality, y = mean_val, color = BAND),
    size = 3.2,
    inherit.aes = FALSE
  ) +
  facet_wrap(~ BAND, nrow = 1, scales = "fixed") +
  scale_color_manual(
    values = c(mod_palette, band_palette),
    breaks = names(band_palette),   # if legend were shown, it would show ONLY bands
    drop = FALSE
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.06, 0.08))) +
  labs(
    #   title = OUTCOME_TITLE,
    x = NULL,              # <-- removes "Modality" axis label
    y = Y_LABEL,
    color = NULL
  ) +
  theme_minimal(base_family = "serif", base_size = BASE_FONT_SIZE) +
  theme(
    plot.title = element_text(size = TITLE_SIZE, hjust = 0.5),
    strip.text = element_text(face = "bold", size = STRIP_SIZE),
    axis.title.y = element_text(size = AXIS_TITLE_SIZE),
    axis.text    = element_text(size = AXIS_TEXT_SIZE),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    legend.position = "none",       # <-- legend handled by dummy legend below
    plot.margin = margin(10, 10, 2, 10)
  )

###############################################################################
# 4) STATS STRIP
###############################################################################
p_stats <- ggplot(
  summary_stats %>% filter(label_txt != ""),
  aes(x = Modality, y = 1, label = label_txt, color = BAND)
) +
  geom_text(size = LABEL_SIZE, family = "serif", lineheight = 0.95, vjust = 1) +
  facet_wrap(~ BAND, nrow = 1) +
  scale_color_manual(values = band_palette, drop = FALSE) +
  scale_y_continuous(limits = c(0.5, 1.1), expand = c(0, 0)) +
  theme_void(base_family = "serif", base_size = BASE_FONT_SIZE) +
  theme(
    strip.text = element_blank(),
    legend.position = "none",
    plot.margin = margin(0, 10, 2, 10)
  )

###############################################################################
# 5) DUMMY LEGEND: band only (bottom)
###############################################################################
p_leg_band <- ggplot(
  summary_stats %>% filter(!is.na(mean_val)),
  aes(x = BAND, y = mean_val, color = BAND)
) +
  geom_point(size = 3) +
  scale_color_manual(values = band_palette, drop = FALSE) +
  labs(color = "Band") +
  theme_void(base_family = "serif", base_size = BASE_FONT_SIZE) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = AXIS_TEXT_SIZE),
    legend.text  = element_text(size = AXIS_TEXT_SIZE)
  )

leg_band <- patchwork::wrap_elements(cowplot::get_legend(p_leg_band))

###############################################################################
# 6) Combine
###############################################################################
p_entropy_final <- (p_main / p_stats / leg_band) +
  plot_layout(heights = c(1, 0.22, 0.14))

p_entropy_final

ggsave("ENTROPY_BY_BANDS_FACET.png",
       p_entropy_final, width = 12, height = 7, dpi = 300, bg = "white")

message("DONE: saved ENTROPY_BY_BANDS_FACET.png")
