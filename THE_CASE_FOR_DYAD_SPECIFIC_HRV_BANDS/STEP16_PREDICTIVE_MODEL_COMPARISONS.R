################################################################################
####################### LIBRARY
library(dplyr)
library(nlme)
library(emmeans)
library(dplyr)
library(broom)
library(stringr)
library(boot)

################################################################################
####################### SET UP ENVIRONMENT

#WORKING DIRECTORY
work_dir <- "C:\\Users\\cjh37695\\Dropbox\\FREQ_OPT_BANDS\\ANALYISIS\\"
setwd(work_dir)

##SUMMARY DF JUST IN CASE 
FOB_SUM <- read.csv("FINAL_FREQUENCY_OPTIMIZED_BAND_SELECTION.csv")

##LOAD IN DATA FOR ANALYSIS 
df <- read.csv("COHERENCE_FULL_DF_WIDE_RAW.csv")


################################################################################
###################### DESCRIPTIVES ON BAND 
names(FOB_SUM)

FOB_SUM <- FOB_SUM %>%
  mutate(
    FOB_hz_lo = pmin(FOB_hz_start, FOB_hz_end, na.rm = TRUE),
    FOB_hz_hi = pmax(FOB_hz_start, FOB_hz_end, na.rm = TRUE),
    
    FOB_band = case_when(
      is.na(FOB_hz_lo) | is.na(FOB_hz_hi) ~ NA_character_,
      
      # fully LF: 0.04 <= f < 0.15
       FOB_hz_hi < 0.15 ~ "LOW",
      
      # fully HF: 0.15 <= f <= 0.40
      FOB_hz_lo >= 0.15  ~ "HIGH",
      
      # straddles LF/HF boundary at 0.15
      FOB_hz_lo < 0.15 & FOB_hz_hi > 0.15 ~ "BOTH",
      
      # optional: outside bands (e.g., below 0.04 or above 0.40)
      TRUE ~ "Out of range"
    )
  )

table(FOB_SUM$FOB_band)
103/177 # 58% LF
6/177 # 3% HF
68/177 # 38% BOTH 


## COMPUTE SAMPLE AVERAGE FOB

FOB_SUM$FOB_MEAN <- rowMeans(
  cbind(
    pmin(FOB_SUM$FOB_hz_start, FOB_SUM$FOB_hz_end),
    pmax(FOB_SUM$FOB_hz_start, FOB_SUM$FOB_hz_end)
  ),
  na.rm = TRUE
)

n  <- sum(!is.na(FOB_SUM$FOB_MEAN))
mean_FOB <- mean(FOB_SUM$FOB_MEAN, na.rm = TRUE)
sd_FOB   <- sd(FOB_SUM$FOB_MEAN, na.rm = TRUE)
se_FOB   <- sd_FOB / sqrt(n)

# 95% CI using t-distribution (preferred)
t_crit <- qt(.975, df = n - 1)

CI_lower <- mean_FOB - t_crit * se_FOB
CI_upper <- mean_FOB + t_crit * se_FOB

mean_FOB
se_FOB
CI_lower
CI_upper

# EFFECT OF FREQUENCY ON BANDWIDTH 
FOB_SUM$FOB_width <- abs(FOB_SUM$FOB_hz_end - FOB_SUM$FOB_hz_start)

M1 <- lm(FOB_width ~ scale(FOB_MEAN, scale = F, center = T), data = FOB_SUM)
summary(M1)
confint(M1)

## DESCRIPTIVES 
names(df)
df_RED <- df %>%
  select(c(ID, GMC_C_AGE, SAMPLE, C_RACE, GMC_INCOME))

FUll_FOB <- left_join(FOB_SUM, df_RED, by = "ID")

FUll_FOB$SAMPLE  <- as.factor(FUll_FOB$SAMPLE)
FUll_FOB$C_RACE  <- as.factor(FUll_FOB$C_RACE)

#SAMPLE 
t.test(FOB_MEAN~SAMPLE, data = FUll_FOB)

#AGE
AGE_cor_test <- cor.test(
  FUll_FOB$GMC_C_AGE,
  FUll_FOB$FOB_MEAN,
  use = "complete.obs"
)

AGE_cor_test

#INCOME 
INC_cor_test <- cor.test(
  FUll_FOB$GMC_INCOME,
  FUll_FOB$FOB_MEAN,
  use = "complete.obs"
)

INC_cor_test

#RACE 
RACE_ANOVA <- aov(FOB_MEAN ~ C_RACE, data = FUll_FOB)
summary(RACE_ANOVA)
################################################################################
###################### RUN SOME MODELS 

## ENSURE CATEGORIES ARE FACTORS 
df$P_SEX   <- as.factor(df$P_SEX)
df$C_SEX   <- as.factor(df$C_SEX)
df$SAMPLE  <- as.factor(df$SAMPLE)
df$C_RACE  <- as.factor(df$C_RACE)



###############################################################################
# FREQUENCY-OPTIMIZED BAND PREDICTIVE VALIDITY
# Hierarchical within-band model ladder + REPEATED 5-fold CV + MSE/RMSE/R2
# NOTE: GMC_* variables are already grand-mean centered
###############################################################################

# -----------------------
# Libraries
# -----------------------
library(dplyr)
library(tibble)

###############################################################################
# 0) SETTINGS: covariates, outcomes, band triplets, model ladder, CV settings
###############################################################################
names(df)

set.seed(125)  # base seed for reproducibility across the whole run

# Ensure categorical covariates are factors (adjust if your coding differs)
df$P_SEX   <- as.factor(df$P_SEX)
df$C_SEX   <- as.factor(df$C_SEX)
df$SAMPLE  <- as.factor(df$SAMPLE)
df$C_RACE  <- as.factor(df$C_RACE)

# toggle covariates on/off
use_covars <- TRUE

covars <- if (use_covars) {
  c("GMC_C_AGE", "C_SEX", "SAMPLE") #, "C_RACE"
} else {
  character(0)
}

# outcomes (simple)
outcomes_simple <- c("ACCp1", "PCONp1", "FCONp1")

# outcomes with baseline controls (time 2 controlling time 1)
outcomes_baseline <- tibble::tribble(
  ~outcome,   ~baseline,
  "INT_P_2",  "GMC_INT_P_1",
  "EXT_P_2",  "GMC_EXT_P_1",
  "TOT_P_2",  "GMC_TOT_P_1",
  "SOM_P_2",  "GMC_SOM_P_1",
  "ANX_P_2",  "GMC_ANX_P_1",
  "WTH_P_2",  "GMC_WTH_P_1",
  "SOC_P_2",  "GMC_SOC_P_1",
  "THO_P_2",  "GMC_THO_P_1",
  "ATT_P_2",  "GMC_ATT_P_1",
  "RBK_P_2",  "GMC_RBK_P_1",
  "AGG_P_2",  "GMC_AGG_P_1"
)

# -----------------------
# Repeated K-fold CV settings (reduces seed-based frailty)
# -----------------------
K_FOLDS <- 5
N_REPEATS <- 100      # typical range: 50–200 for N~80–130
MIN_N <- 10

# -----------------------
# Band triplets (COH and PHS are treated as separate "band sets")
# Each set provides (mu, a, b) for that band/metric combo.
# -----------------------
band_triplets <- list(
  COH_FOB = c(mu = "GMC_COH_MU_FM_FOB", a = "GMC_COH_A_FM_FOB", b = "GMC_COH_B_FM_FOB"),
  PHS_FOB = c(mu = "GMC_PHS_MU_FM_FOB", a = "GMC_PHS_A_FM_FOB", b = "GMC_PHS_B_FM_FOB"),
  
  COH_RSB = c(mu = "GMC_COH_MU_FM_RSB", a = "GMC_COH_A_FM_RSB", b = "GMC_COH_B_FM_RSB"),
  PHS_RSB = c(mu = "GMC_PHS_MU_FM_RSB", a = "GMC_PHS_A_FM_RSB", b = "GMC_PHS_B_FM_RSB"),
  
  COH_HF  = c(mu = "GMC_COH_MU_FM_HF",  a = "GMC_COH_A_FM_HF",  b = "GMC_COH_B_FM_HF"),
  PHS_HF  = c(mu = "GMC_PHS_MU_FM_HF",  a = "GMC_PHS_A_FM_HF",  b = "GMC_PHS_B_FM_HF"),
  
  COH_LF  = c(mu = "GMC_COH_MU_FM_LF",  a = "GMC_COH_A_FM_LF",  b = "GMC_COH_B_FM_LF"),
  PHS_LF  = c(mu = "GMC_PHS_MU_FM_LF",  a = "GMC_PHS_A_FM_LF",  b = "GMC_PHS_B_FM_LF"),
  
  COH_RND_FOB = c(mu = "GMC_COH_MU_RD_FOB", a = "GMC_COH_A_RD_FOB", b = "GMC_COH_B_RD_FOB"),
  PHS_RND_FOB = c(mu = "GMC_PHS_MU_RD_FOB", a = "GMC_PHS_A_RD_FOB", b = "GMC_PHS_B_RD_FOB")
)

# -----------------------
# Model ladder generator within each triplet
# - Main-effects-only models (single, pair, all three)
# - Interaction models for each pair (with their main effects)
# - Optional full 2-way interaction model
# -----------------------
make_model_ladder <- function(trip) {
  a  <- trip[["a"]]
  b  <- trip[["b"]]
  mu <- trip[["mu"]]
  
  list(
    M1_a       = c(a),
    M2_b       = c(b),
    M3_mu      = c(mu),
    
    M4_a_b     = c(a, b),
    M5_a_b_mu  = c(a, b, mu),
    
    # interaction models include main effects (hierarchy principle)
    M6_aXb     = c(a, b, paste0(a, ":", b)),
    M7_aXmu    = c(a, mu, paste0(a, ":", mu)),
    M8_bXmu    = c(b, mu, paste0(b, ":", mu)),
    
    # full pairwise interaction model (all 2-way interactions)
    M9_all2way = c(a, b, mu,
                   paste0(a, ":", b),
                   paste0(a, ":", mu),
                   paste0(b, ":", mu))
  )
}

model_specs <- lapply(band_triplets, make_model_ladder)

###############################################################################
# 1) RESULTS OBJECT
# Stores MEAN performance across repeats; also stores SD to quantify stability.
###############################################################################
results <- data.frame(
  outcome         = character(),
  baseline        = character(),
  band_set        = character(),
  model_id        = character(),
  n_used          = integer(),
  
  mse_cv_mean     = numeric(),
  mse_cv_sd       = numeric(),
  
  rmse_cv_mean    = numeric(),
  rmse_cv_sd      = numeric(),
  
  r2_cv_mean      = numeric(),
  r2_cv_sd        = numeric(),
  
  r2_fullsample   = numeric(),
  stringsAsFactors = FALSE
)

###############################################################################
# 2) REPEATED K-FOLD CV RUNNER (LM)
###############################################################################
run_repeated_kfold_lm <- function(dat, fml, K = 5, repeats = 100) {
  n <- nrow(dat)
  y_name <- all.vars(fml)[1]
  y_obs <- dat[[y_name]]
  
  mse_vec <- numeric(repeats)
  rmse_vec <- numeric(repeats)
  r2_vec <- numeric(repeats)
  
  for (r in seq_len(repeats)) {
    
    folds <- sample(rep(1:K, length.out = n))
    y_hat <- rep(NA_real_, n)
    
    for (k in 1:K) {
      train <- dat[folds != k, , drop = FALSE]
      test  <- dat[folds == k, , drop = FALSE]
      
      m <- lm(fml, data = train)
      y_hat[folds == k] <- predict(m, newdata = test)
    }
    
    mse <- mean((y_obs - y_hat)^2)
    rmse <- sqrt(mse)
    
    sse <- sum((y_obs - y_hat)^2)
    sst <- sum((y_obs - mean(y_obs))^2)
    r2_cv <- 1 - (sse / sst)
    
    mse_vec[r] <- mse
    rmse_vec[r] <- rmse
    r2_vec[r] <- r2_cv
  }
  
  # Full-sample R2 (fit once; descriptive)
  m_full <- lm(fml, data = dat)
  r2_full <- summary(m_full)$r.squared
  
  list(
    mse_mean  = mean(mse_vec),
    mse_sd    = sd(mse_vec),
    
    rmse_mean = mean(rmse_vec),
    rmse_sd   = sd(rmse_vec),
    
    r2_mean   = mean(r2_vec),
    r2_sd     = sd(r2_vec),
    
    r2_full   = r2_full
  )
}

###############################################################################
# 3) LOOP BLOCK A: SIMPLE OUTCOMES (no baseline control)
###############################################################################
for (y in outcomes_simple) {
  
  for (band_set in names(model_specs)) {
    
    for (model_id in names(model_specs[[band_set]])) {
      
      preds <- model_specs[[band_set]][[model_id]]
      rhs <- paste(c(preds, covars), collapse = " + ")
      fml <- as.formula(paste(y, "~", rhs))
      
      model_vars <- unique(c(y, covars, all.vars(fml)[-1]))
      dat <- df[complete.cases(df[, model_vars]), model_vars, drop = FALSE]
      
      n <- nrow(dat)
      if (n < MIN_N) next
      
      met <- run_repeated_kfold_lm(dat, fml, K = K_FOLDS, repeats = N_REPEATS)
      
      results <- rbind(
        results,
        data.frame(
          outcome       = y,
          baseline      = NA_character_,
          band_set      = band_set,
          model_id      = model_id,
          n_used        = n,
          
          mse_cv_mean   = met$mse_mean,
          mse_cv_sd     = met$mse_sd,
          
          rmse_cv_mean  = met$rmse_mean,
          rmse_cv_sd    = met$rmse_sd,
          
          r2_cv_mean    = met$r2_mean,
          r2_cv_sd      = met$r2_sd,
          
          r2_fullsample = met$r2_full,
          stringsAsFactors = FALSE
        )
      )
    }
  }
}

###############################################################################
# 4) LOOP BLOCK B: BASELINE-CONTROL OUTCOMES
###############################################################################
for (row in 1:nrow(outcomes_baseline)) {
  
  y  <- outcomes_baseline$outcome[row]
  b0 <- outcomes_baseline$baseline[row]
  
  for (band_set in names(model_specs)) {
    
    for (model_id in names(model_specs[[band_set]])) {
      
      preds <- model_specs[[band_set]][[model_id]]
      rhs <- paste(c(b0, preds, covars), collapse = " + ")
      fml <- as.formula(paste(y, "~", rhs))
      
      model_vars <- unique(c(y, b0, covars, all.vars(fml)[-1]))
      dat <- df[complete.cases(df[, model_vars]), model_vars, drop = FALSE]
      
      n <- nrow(dat)
      if (n < MIN_N) next
      
      met <- run_repeated_kfold_lm(dat, fml, K = K_FOLDS, repeats = N_REPEATS)
      
      results <- rbind(
        results,
        data.frame(
          outcome       = y,
          baseline      = b0,
          band_set      = band_set,
          model_id      = model_id,
          n_used        = n,
          
          mse_cv_mean   = met$mse_mean,
          mse_cv_sd     = met$mse_sd,
          
          rmse_cv_mean  = met$rmse_mean,
          rmse_cv_sd    = met$rmse_sd,
          
          r2_cv_mean    = met$r2_mean,
          r2_cv_sd      = met$r2_sd,
          
          r2_fullsample = met$r2_full,
          stringsAsFactors = FALSE
        )
      )
    }
  }
}

###############################################################################
# 5) VIEW + SORT
###############################################################################
results_sorted <- results %>%
  arrange(outcome, baseline, mse_cv_mean)   # sort by best (lowest) MSE within outcome/baseline

options(max.print = 10000)
print(results_sorted, row.names = FALSE)

###############################################################################
# 6) FILTERED + COLORED DISPLAY
# Keep ONLY models that are within the MSE equivalence band (red + yellow).
# Then, within that reduced set, re-compute and highlight:
# BLUE   = best CV-R2 within the reduced set
# PURPLE = best full-sample R2 within the reduced set
###############################################################################

# --- ANSI helpers ---
ansi_red    <- function(x) paste0("\033[31m", x, "\033[0m")
ansi_blue   <- function(x) paste0("\033[34m", x, "\033[0m")
ansi_yellow <- function(x) paste0("\033[33m", x, "\033[0m")
ansi_purple <- function(x) paste0("\033[35m", x, "\033[0m")

ansi_bold   <- function(x) paste0("\033[1m",  x, "\033[0m")

red_bold    <- function(x) ansi_red(ansi_bold(x))
blue_bold   <- function(x) ansi_blue(ansi_bold(x))
purple_bold <- function(x) ansi_purple(ansi_bold(x))
yellow      <- function(x) ansi_yellow(x)

fmt <- function(x, digits = 4) formatC(x, format = "f", digits = digits)

strip_ansi <- function(x) gsub("\033\\[[0-9;]*m", "", x)

print_ansi_table <- function(df, sep = "  ") {
  
  df <- as.data.frame(df, check.names = FALSE)
  df[is.na(df)] <- ""
  
  col_widths <- sapply(names(df), function(col) {
    max(
      nchar(strip_ansi(col)),
      max(nchar(strip_ansi(as.character(df[[col]]))), na.rm = TRUE)
    )
  })
  
  pad_vis <- function(x, width) {
    x <- as.character(x)
    vis_len <- nchar(strip_ansi(x))
    paste0(x, strrep(" ", max(0, width - vis_len)))
  }
  
  header <- mapply(pad_vis, names(df), col_widths)
  cat(paste(header, collapse = sep), "\n")
  
  cat(paste(mapply(function(w) strrep("-", w), col_widths), collapse = sep), "\n")
  
  for (i in seq_len(nrow(df))) {
    row <- mapply(pad_vis, df[i, ], col_widths)
    cat(paste(row, collapse = sep), "\n")
  }
}

###############################################################################
# 6) COLORED DISPLAY WITH SEQUENTIAL FILTERING LADDER
#
# Stage 1: CV-MSE
#   RED    = best MSE
#   YELLOW = within 0.5 SD of best MSE
#
# Stage 2: CV-R² (ONLY among Stage-1 models)
#   BLUE   = best CV-R²
#   ORANGE = within 0.5 SD of best CV-R²
#
# Stage 3: Full-sample R² (ONLY among Stage-2 models)
#   PURPLE = best full-sample R²
###############################################################################

# --- ANSI helpers ---
ansi_red    <- function(x) paste0("\033[31m", x, "\033[0m")
ansi_yellow <- function(x) paste0("\033[33m", x, "\033[0m")
ansi_blue   <- function(x) paste0("\033[34m", x, "\033[0m")
ansi_orange <- function(x) paste0("\033[38;5;208m", x, "\033[0m")
ansi_purple <- function(x) paste0("\033[35m", x, "\033[0m")

ansi_bold   <- function(x) paste0("\033[1m", x, "\033[0m")

red_bold    <- function(x) ansi_red(ansi_bold(x))
blue_bold   <- function(x) ansi_blue(ansi_bold(x))
purple_bold <- function(x) ansi_purple(ansi_bold(x))
yellow      <- function(x) ansi_yellow(x)
orange      <- function(x) ansi_orange(x)

fmt <- function(x, digits = 4) formatC(x, format = "f", digits = digits)

###############################################################################
# Sequential filtering ladder (print ONLY highlighted rows)
# + GREEN = bottom 3 (smallest) CV-MSE SD within the Stage-1 (MSE-equivalent) set
###############################################################################

# --- ANSI helpers (add green + orange if not already defined) ---
ansi_green  <- function(x) paste0("\033[32m", x, "\033[0m")
ansi_orange <- function(x) paste0("\033[38;5;208m", x, "\033[0m")  # 256-color orange

green       <- function(x) ansi_green(x)
orange      <- function(x) ansi_orange(x)

green_bold  <- function(x) ansi_green(ansi_bold(x))

results_colored <- results_sorted %>%
  group_by(outcome, baseline) %>%
  mutate(
    
    ###################################
    # STAGE 1: CV-MSE filter
    ###################################
    best_mse = min(mse_cv_mean, na.rm = TRUE),
    
    # SD from best-MSE row
    mse_sd_best = mse_cv_sd[which.min(mse_cv_mean)],
    
    mse_equiv_threshold = best_mse + .25 * mse_sd_best,
    
    mse_is_best  = mse_cv_mean == best_mse,
    mse_is_equiv = mse_cv_mean <= mse_equiv_threshold,
    
    ###################################
    # STAGE 1b: Highlight smallest 3 MSE SDs (within Stage-1 set)
    ###################################
    mse_sd_rank_stage1 = ifelse(
      mse_is_equiv,
      dplyr::min_rank(mse_cv_sd),   # 1 = smallest; ties share rank
      NA_integer_
    ),
    
    mse_sd_is_best3 = !is.na(mse_sd_rank_stage1) & mse_sd_rank_stage1 <= 3,
    
    ###################################
    # STAGE 2: CV-R² filter (only within Stage-1)
    ###################################
    best_r2_stage1 = ifelse(
      any(mse_is_equiv, na.rm = TRUE),
      max(r2_cv_mean[mse_is_equiv], na.rm = TRUE),
      NA_real_
    ),
    
    idx_best_r2_stage1 = ifelse(
      any(mse_is_equiv, na.rm = TRUE),
      which.max(ifelse(mse_is_equiv, r2_cv_mean, -Inf)),
      NA_integer_
    ),
    
    # SD from best-R² row (so orange actually shows up)
    r2_sd_best = ifelse(
      is.na(idx_best_r2_stage1),
      NA_real_,
      r2_cv_sd[idx_best_r2_stage1]
    ),
    
    r2_equiv_threshold = ifelse(
      is.na(best_r2_stage1) | is.na(r2_sd_best),
      NA_real_,
      best_r2_stage1 - .25 * r2_sd_best
    ),
    
    r2_is_best  = mse_is_equiv & !is.na(best_r2_stage1) & (r2_cv_mean == best_r2_stage1),
    r2_is_equiv = mse_is_equiv & !is.na(r2_equiv_threshold) & (r2_cv_mean >= r2_equiv_threshold),
    
    ###################################
    # STAGE 3: Full-sample R² (only within Stage-2)
    ###################################
    best_r2full_stage2 = ifelse(
      any(r2_is_equiv, na.rm = TRUE),
      max(r2_fullsample[r2_is_equiv], na.rm = TRUE),
      NA_real_
    ),
    
    r2full_is_best =
      !is.na(best_r2full_stage2) &
      r2_is_equiv &
      (r2_fullsample == best_r2full_stage2)
    
  ) %>%
  ungroup() %>%
  mutate(
    
    ###################################
    # Display formatting
    ###################################
    mse_disp =
      ifelse(mse_is_best,
             red_bold(fmt(mse_cv_mean, 4)),
             ifelse(mse_is_equiv,
                    yellow(fmt(mse_cv_mean, 4)),
                    fmt(mse_cv_mean, 4))),
    
    rmse_disp = fmt(rmse_cv_mean, 4),
    
    r2_disp =
      ifelse(r2_is_best,
             blue_bold(fmt(r2_cv_mean, 4)),
             ifelse(r2_is_equiv,
                    orange(fmt(r2_cv_mean, 4)),
                    fmt(r2_cv_mean, 4))),
    
    r2full_disp =
      ifelse(r2full_is_best,
             purple_bold(fmt(r2_fullsample, 4)),
             fmt(r2_fullsample, 4)),
    
    # ✅ GREEN highlight for the smallest 3 MSE SDs within Stage-1 set
    mse_sd_disp =
      ifelse(mse_sd_is_best3,
             green_bold(fmt(mse_cv_sd, 4)),
             fmt(mse_cv_sd, 4)),
    
    rmse_sd_disp = fmt(rmse_cv_sd, 4),
    r2_sd_disp   = fmt(r2_cv_sd, 4)
    
  ) %>%
  # ✅ KEEP ONLY ROWS THAT ARE "HIGHLIGHTED" IN ANY WAY (including green SD rows)
  filter(mse_is_equiv | r2_is_equiv | r2full_is_best | mse_sd_is_best3) %>%
  select(
    outcome, baseline, band_set, model_id, n_used,
    mse_disp, mse_sd_disp,
    rmse_disp, rmse_sd_disp,
    r2_disp, r2_sd_disp,
    r2full_disp
  )

options(max.print = 10000000)
print_ansi_table(results_colored)

## SAVE RESULTS 
RES_TAB <- results_colored %>%
  mutate(across(everything(), strip_ansi)) %>%
  as.data.frame()

write.csv(RES_TAB,
          "RAW_SUMMARY_BAND_CONTRASTS.csv",
          row.names = FALSE,
          na = "")

################################################################################
################################ REGRESSION MODELS 
library(glmmTMB)
library(parameters)
options(scipen = 999)
names(df)


## PSYCHOLOGICAL CONTROL 
PCON <- lm(PCONp1 ~ GMC_C_AGE + C_SEX + SAMPLE +
             GMC_COH_A_FM_FOB + GMC_PHS_MU_FM_FOB, data = df)

summary(PCON)
model_parameters(
  PCON,
  standardize = "refit",   # TRUE standardization (recommended)
  ci = 0.95
)

PCON_RND <- glmmTMB(
  PCONp1 ~ GMC_C_AGE + C_SEX +
    GMC_COH_A_FM_FOB + GMC_PHS_MU_FM_FOB +
    (1 | SAMPLE),
  data = df,
  family = gaussian()
)

summary(PCON_RND)
model_parameters(
  PCON_RND,
  standardize = "refit",   # TRUE standardization (recommended)
  ci = 0.95
)

## Anxiety 
Anx <- lm(ANX_P_2 ~ GMC_ANX_P_1 + GMC_C_AGE + C_SEX + SAMPLE + 
               GMC_PHS_A_FM_RSB + 
               GMC_PHS_B_FM_RSB +
               (GMC_PHS_A_FM_RSB*GMC_PHS_B_FM_RSB) 
               GMC_COH_MU_FM_RSB, data = df)

summary(Anx)
model_parameters(
  Anx,
  standardize = "refit",   # TRUE standardization (recommended)
  ci = 0.95
)

Anx_RND <- glmmTMB(
  ANX_P_2 ~ GMC_ANX_P_1 + GMC_C_AGE + C_SEX +
    GMC_PHS_A_FM_RSB + 
    GMC_PHS_B_FM_RSB +
    (GMC_PHS_A_FM_RSB*GMC_PHS_B_FM_RSB)  + 
    (1 | SAMPLE),  data = df,
  family = gaussian()
)


summary(Anx_RND)
model_parameters(
  Anx_RND,
  standardize = "refit",   # TRUE standardization (recommended)
  ci = 0.95
)
names(df)


## Aggression
AGGRES <- lm(AGG_P_2 ~ GMC_AGG_P_1 + GMC_C_AGE + C_SEX + SAMPLE + 
             GMC_COH_A_FM_FOB + 
             GMC_COH_MU_FM_FOB +
             (GMC_COH_A_FM_FOB*GMC_COH_MU_FM_FOB) + 
               GMC_PHS_A_FM_FOB +
               GMC_PHS_A_FM_HF +
               GMC_PHS_B_FM_HF + 
               GMC_PHS_MU_FM_HF +
               GMC_COH_B_FM_HF + 
               GMC_PHS_B_FM_LF +
               GMC_PHS_B_FM_RSB +
               GMC_COH_B_FM_RSB + 
               GMC_COH_MU_FM_RSB, data = df)

summary(AGGRES)
model_parameters(
  INTERN,
  standardize = "refit",   # TRUE standardization (recommended)
  ci = 0.95
)

AGGRES_RND <- glmmTMB(
  AGG_P_2 ~ GMC_AGG_P_1 + GMC_C_AGE + C_SEX + 
    GMC_COH_A_FM_FOB + 
    GMC_COH_MU_FM_FOB +
    (GMC_COH_A_FM_FOB*GMC_COH_MU_FM_FOB) + 
    GMC_PHS_A_FM_FOB +
    GMC_PHS_A_FM_HF +
    GMC_PHS_B_FM_HF + 
    GMC_PHS_MU_FM_HF +
    GMC_COH_B_FM_HF + 
    GMC_PHS_B_FM_LF +
    GMC_PHS_B_FM_RSB +
    GMC_COH_B_FM_RSB + 
    GMC_COH_MU_FM_RSB + 
    (1 | SAMPLE),  data = df,
  family = gaussian()
)


summary(AGGRES_RND)
model_parameters(
  AGGRES_RND,
  standardize = "refit",   # TRUE standardization (recommended)
  ci = 0.95
)
names(df)

## Attention
ATTPROB <- lm(ATT_P_2 ~  GMC_ATT_P_1 + GMC_C_AGE + C_SEX + SAMPLE + 
               GMC_PHS_A_FM_FOB + 
                GMC_PHS_MU_FM_FOB +  
                (GMC_PHS_A_FM_FOB*GMC_PHS_MU_FM_FOB), data = df)

summary(ATTPROB)
model_parameters(
  EXTERN,
  standardize = "refit",   # TRUE standardization (recommended)
  ci = 0.95
)

ATTPROB_RND <- glmmTMB(
  ATT_P_2 ~  GMC_ATT_P_1 + GMC_C_AGE + C_SEX + 
    GMC_PHS_A_FM_FOB + 
    GMC_PHS_MU_FM_FOB +  
    (GMC_PHS_A_FM_FOB*GMC_PHS_MU_FM_FOB) + 
    (1 | SAMPLE),
  data = df,
  family = gaussian()
)

summary(ATTPROB_RND)
model_parameters(
  ATTPROB_RND,
  standardize = "refit",   # TRUE standardization (recommended)
  ci = 0.95
)

names(df)
## RULE BREAKING  
RULEBREAK <- lm(RBK_P_2 ~ GMC_RBK_P_1 + GMC_C_AGE + C_SEX + SAMPLE + 
            GMC_PHS_A_FM_HF, data = df)

summary(RULEBREAK)
model_parameters(
  RULEBREAK,
  standardize = "refit",   # TRUE standardization (recommended)
  ci = 0.95
)

RULEBREAK_RND <- glmmTMB(
  RBK_P_2 ~ GMC_RBK_P_1 + GMC_C_AGE + C_SEX + 
    GMC_PHS_A_FM_HF  + 
    (1 | SAMPLE),  data = df,
  family = gaussian()
)


summary(RULEBREAK_RND)
model_parameters(
  RULEBREAK_RND,
  standardize = "refit",   # TRUE standardization (recommended)
  ci = 0.95
)
names(df)

## SOCIAL PROBLEMS 
SOCP <- lm(SOC_P_2 ~ GMC_SOC_P_1 + GMC_C_AGE + C_SEX + SAMPLE +
             GMC_PHS_MU_FM_FOB, data = df)

summary(SOCP)
model_parameters(
  SOCP,
  standardize = "refit",   # TRUE standardization (recommended)
  ci = 0.95
)

SOCP_RND <- glmmTMB(
  SOC_P_2 ~ GMC_SOC_P_1 + GMC_C_AGE + C_SEX + 
    GMC_PHS_MU_FM_FOB + 
    (1 | SAMPLE),
  data = df,
  family = gaussian()
)

summary(SOCP_RND)
model_parameters(
  SOCP_RND,
  standardize = "refit",   # TRUE standardization (recommended)
  ci = 0.95
)
names(df)

hist(df$GMC_PHS_MU_FM_FOB)


# WITHDRAWN 
WITHDRAW <- lm(WTH_P_2 ~ GMC_WTH_P_1 + GMC_C_AGE + C_SEX + SAMPLE +
             GMC_PHS_MU_FM_LF +
             GMC_COH_B_FM_HF, data = df)

summary(WITHDRAW)
model_parameters(
  WITHDRAW,
  standardize = "refit",   # TRUE standardization (recommended)
  ci = 0.95
)

WITHDRAW_RND <- glmmTMB(
  WTH_P_2 ~ GMC_WTH_P_1 + GMC_C_AGE + C_SEX + 
    GMC_PHS_MU_FM_LF +
    GMC_COH_B_FM_HF + 
    (1 | SAMPLE),
  data = df,
  family = gaussian()
)

summary(WITHDRAW_RND)
model_parameters(
  WITHDRAW_RND,
  standardize = "refit",   # TRUE standardization (recommended)
  ci = 0.95
)
names(df)

hist(df$GMC_PHS_MU_FM_FOB)


# Thought Problems 
WITHDRAW <- lm(THO_P_2 ~ GMC_THO_P_1 + GMC_C_AGE + C_SEX + SAMPLE + 
                 GMC_PHS_MU_FM_FOB +
                 GMC_COH_MU_FM_FOB + 
                 GMC_COH_B_FM_LF + 
                 GMC_COH_MU_FM_RSB, data = df)

summary(WITHDRAW)
model_parameters(
  WITHDRAW,
  standardize = "refit",   # TRUE standardization (recommended)
  ci = 0.95
)

WITHDRAW_RND <- glmmTMB(
  THO_P_2 ~ GMC_THO_P_1 + GMC_C_AGE + C_SEX + 
    GMC_PHS_MU_FM_FOB +
    GMC_COH_MU_FM_FOB + 
    GMC_COH_B_FM_LF + 
    GMC_COH_MU_FM_RSB +
    (1 | SAMPLE),
  data = df,
  family = gaussian()
)

summary(WITHDRAW_RND)
model_parameters(
  WITHDRAW_RND,
  standardize = "refit",   # TRUE standardization (recommended)
  ci = 0.95
)
names(df)


# TOTAL PROBLEMS 

TOTAL_P <- lm(GMC_TOT_P_1 + GMC_C_AGE + C_SEX + SAMPLE + 
                GMC_PHS_MU_FM_FOB +
                GMC_COH_B_FM_LF + 
                GMC_COH_B_FM_HF, data = df)

summary(TOTAL_P)
model_parameters(
  TOTAL_P,
  standardize = "refit",   # TRUE standardization (recommended)
  ci = 0.95
)

TOTAL_P_RND <- glmmTMB(
  TOT_P_2 ~ GMC_TOT_P_1 + GMC_C_AGE + C_SEX + 
    GMC_PHS_MU_FM_FOB +
    GMC_COH_B_FM_LF + 
    GMC_COH_B_FM_HF + 
    (1 | SAMPLE),
  data = df,
  family = gaussian()
)

summary(TOTAL_P_RND)
model_parameters(
  TOTAL_P_RND,
  standardize = "refit",   # TRUE standardization (recommended)
  ci = 0.95
)
names(df)
hist(df$GMC_PHS_MU_FM_FOB)


sim <- simulateResiduals(SOCP_RND)
plot(sim)
