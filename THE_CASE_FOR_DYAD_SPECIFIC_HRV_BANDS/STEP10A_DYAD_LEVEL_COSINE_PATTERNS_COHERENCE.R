###############################################################################
# DYAD-SPECIFIC COSINE FITS (COHERENCE) — WITH BOUNDS VALIDATION & QC PLOTS
# UPDATED: constrain amplitude a to be NON-NEGATIVE (a >= 0)
#
# Key changes:
# - A_MIN set to 0 (instead of -0.5)
# - a_multiplier replaced with positive-only multipliers (includes 0)
# - a_grid uses abs(sdy) * multipliers and is clipped to [A_MIN, A_MAX]
# - is_valid_params updated to use mu ± a (no abs()) since a >= 0
# - (optional) adds COH_C_CANON = c %% b for interpretability (saved in summary only)
###############################################################################

library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)

###############################################################################
# 0) PATHS + LOAD

work_dir <- "/home/cjh37695/WAVELET_WORKING/TEST_RESP"
setwd(work_dir)

coh_file <- "SAMPLE_TRUNC_COHERENCE_DS_SERIES.csv"
COH_WIDE <- read.csv(coh_file, check.names = FALSE, stringsAsFactors = FALSE)

# Drop accidental rowname column
COH_WIDE <- COH_WIDE[, -1]

## REMOVE DYADS FAILING QC CHECKS
BAD_IDS <- c(1035)  # VISUAL INSPECTION
COH_WIDE$FAIL_QC <- ifelse(COH_WIDE$ID %in% BAD_IDS, 1, 0)
COH_WIDE <- subset(COH_WIDE, FAIL_QC == 0)

# Drop unnamed/blank columns safely
nms <- names(COH_WIDE)
bad <- is.na(nms) | trimws(nms) == ""
if (any(bad)) {
  message("Dropping ", sum(bad), " columns with NA/blank names.")
  COH_WIDE <- COH_WIDE[, !bad, drop = FALSE]
}

###############################################################################
# 1) SETTINGS + BOUNDS DEFINITION

bands  <- c("FOB","RSB","HF","LF")
b_grid <- seq(20, 240, by = 5)  # period grid in seconds

# Multi-start grids (Option B: big grid; sensitive-start safe)
# UPDATED: positive-only multipliers (include 0 so flat solutions are reachable)
a_multiplier_pos <- c(0, 0.1, 0.25, 0.5, 1, 2)

# c starts: dense across observed window (0..tmax)
c_grid_len <- 9L

# nls controls
nls_maxiter <- 200
nls_tol     <- 1e-6

# PARAMETER BOUNDS (coherence-specific)
MU_MIN <- 0
MU_MAX <- 1

# UPDATED: constrain amplitude to be non-negative
A_MIN  <- 0
A_MAX  <- 0.5

B_MIN  <- 20
B_MAX  <- 600  # twice task length to allow very slow patterns

message("=== PARAMETER BOUNDS ===")
message("mu: [", MU_MIN, ", ", MU_MAX, "]")
message("a:  [", A_MIN, ", ", A_MAX, "] (a >= 0)")
message("b:  [", B_MIN, ", ", B_MAX, "] seconds")
message("c:  [0, max_time] (per dyad)")

###############################################################################
# 2) HELPERS

# Wide (one row ID×BAND) -> Long (one row per timepoint)
row_to_long <- function(one_row_df) {
  time_cols <- names(one_row_df)[str_detect(names(one_row_df), "^t[_\\.]")]
  if (length(time_cols) == 0) stop("No time columns found matching ^t_ / ^t.")
  
  one_row_df %>%
    dplyr::select(ID, BAND, period_start, period_end, all_of(time_cols)) %>%
    tidyr::pivot_longer(
      cols      = all_of(time_cols),
      names_to  = "t_name",
      values_to = "y"
    ) %>%
    dplyr::mutate(
      t = suppressWarnings(as.numeric(sub("^t[_\\.]", "", t_name)))
    ) %>%
    dplyr::filter(!is.na(t), !is.na(y))  # drop missingness
}

compute_r2 <- function(y, yhat) {
  ss_res <- sum((y - yhat)^2, na.rm = TRUE)
  ss_tot <- sum((y - mean(y, na.rm = TRUE))^2, na.rm = TRUE)
  if (isTRUE(all.equal(ss_tot, 0))) return(NA_real_)
  1 - (ss_res / ss_tot)
}

# UPDATED: Validate if parameters are within bounds (a >= 0, so no abs())
is_valid_params <- function(mu, a, b, c, c_max) {
  if (!all(is.finite(c(mu, a, b, c)))) return(FALSE)
  if (mu < MU_MIN || mu > MU_MAX) return(FALSE)
  if (a < A_MIN  || a  > A_MAX)  return(FALSE)
  if (b < B_MIN  || b  > B_MAX)  return(FALSE)
  if (c < 0      || c  > c_max)  return(FALSE)
  # Ensure mu ± a stays in [0,1]
  if (mu + a > 1) return(FALSE)
  if (mu - a < 0) return(FALSE)
  return(TRUE)
}

# Extract full coefficient table + residual sigma + fit indices from an nls/nlsLM
extract_nls_full <- function(mod, dat) {
  sm <- tryCatch(summary(mod), error = function(e) NULL)
  if (is.null(sm)) return(NULL)
  
  ct <- sm$coefficients
  if (is.null(ct)) return(NULL)
  
  sigma <- tryCatch(sm$sigma, error = function(e) NA_real_)
  
  yhat <- tryCatch(predict(mod, newdata = dat), error = function(e) rep(NA_real_, nrow(dat)))
  resid_vec <- dat$y - yhat
  
  SSE <- sum(resid_vec^2, na.rm = TRUE)
  ss_tot <- sum((dat$y - mean(dat$y, na.rm = TRUE))^2, na.rm = TRUE)
  R2 <- ifelse(isTRUE(all.equal(ss_tot, 0)), NA_real_, 1 - (SSE / ss_tot))
  
  list(
    ct = ct,
    sigma = sigma,
    SSE = SSE,
    R2 = R2,
    AIC = tryCatch(AIC(mod), error = function(e) NA_real_),
    BIC = tryCatch(BIC(mod), error = function(e) NA_real_),
    logLik = tryCatch(as.numeric(logLik(mod)), error = function(e) NA_real_)
  )
}

# Multi-start dyad×band cosine fit WITH BOUNDS VALIDATION
# Returns best IN-BOUNDS fit (or NULL if none found)
fit_cosine_multistart <- function(dat, b_grid, c_grid_len = 9L) {
  
  if (nrow(dat) < 10) {
    return(list(converged = FALSE, best = NULL, reason = "too_few_obs"))
  }
  
  mu0 <- mean(dat$y, na.rm = TRUE)
  sdy <- sd(dat$y, na.rm = TRUE)
  if (!is.finite(sdy) || sdy == 0) sdy <- 0.01
  
  # UPDATED: positive-only a grid, clipped into [A_MIN, A_MAX]
  a_grid <- pmin(A_MAX, pmax(A_MIN, abs(sdy) * a_multiplier_pos))
  
  tmax <- max(dat$t, na.rm = TRUE)
  if (!is.finite(tmax) || tmax <= 0) tmax <- max(dat$t, na.rm = TRUE)
  
  c_grid <- unique(round(seq(0, tmax, length.out = c_grid_len), 6))
  
  start_grid <- expand.grid(
    mu = mu0,
    a  = a_grid,
    b  = b_grid,
    c  = c_grid,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  
  use_nlsLM <- requireNamespace("minpack.lm", quietly = TRUE)
  
  valid_fits <- list()
  
  # UPDATED bounds: a lower is 0
  lower <- c(mu = MU_MIN, a = A_MIN, b = B_MIN, c = 0)
  upper <- c(mu = MU_MAX, a = A_MAX, b = B_MAX, c = tmax)
  
  fml <- y ~ mu + a * cos(((2*pi)/b) * (t - c))
  
  for (i in seq_len(nrow(start_grid))) {
    st <- as.list(start_grid[i, ])
    
    mod <- tryCatch({
      if (use_nlsLM) {
        minpack.lm::nlsLM(
          formula = fml,
          data    = dat,
          start   = st,
          lower   = lower,
          upper   = upper,
          control = minpack.lm::nls.lm.control(
            maxiter = nls_maxiter,
            ftol    = nls_tol,
            ptol    = nls_tol
          )
        )
      } else {
        # Fallback: nls doesn't support bounds, so we'll validate after
        nls(
          formula = fml,
          data    = dat,
          start   = st,
          control = nls.control(maxiter = nls_maxiter, tol = nls_tol, warnOnly = TRUE)
        )
      }
    }, error = function(e) NULL)
    
    if (!is.null(mod)) {
      coefs <- tryCatch(coef(mod), error = function(e) NULL)
      if (is.null(coefs)) next
      
      mu <- unname(coefs["mu"])
      a  <- unname(coefs["a"])
      b  <- unname(coefs["b"])
      c  <- unname(coefs["c"])
      
      if (is_valid_params(mu, a, b, c, tmax)) {
        aic <- tryCatch(AIC(mod), error = function(e) NA_real_)
        if (is.finite(aic)) {
          valid_fits[[length(valid_fits) + 1]] <- list(
            mod = mod,
            mu = mu,
            a  = a,
            b  = b,
            c  = c,
            AIC = aic
          )
        }
      }
    }
  }
  
  if (length(valid_fits) == 0) {
    return(list(converged = FALSE, best = NULL, reason = "no_valid_in_bounds_fit"))
  }
  
  best_idx <- which.min(sapply(valid_fits, `[[`, "AIC"))
  best_fit_obj <- valid_fits[[best_idx]]
  best_fit <- best_fit_obj$mod
  
  yhat <- tryCatch(predict(best_fit, newdata = dat),
                   error = function(e) rep(NA_real_, nrow(dat)))
  SSE  <- sum((dat$y - yhat)^2, na.rm = TRUE)
  
  # Optional: canonicalize c for interpretability
  c_canon <- if (is.finite(best_fit_obj$b) && best_fit_obj$b > 0) (best_fit_obj$c %% best_fit_obj$b) else NA_real_
  
  list(
    converged = TRUE,
    best = best_fit,
    mu = best_fit_obj$mu,
    a  = best_fit_obj$a,
    b  = best_fit_obj$b,
    c  = best_fit_obj$c,
    c_canon = c_canon,
    R2 = compute_r2(dat$y, yhat),
    AIC = best_fit_obj$AIC,
    BIC = tryCatch(BIC(best_fit), error = function(e) NA_real_),
    SSE = SSE
  )
}

###############################################################################
# 3) CREATE PLOT DIRECTORY

plot_dir <- "COHERENCE_COSINE_FIT_PLOTS"
if (!dir.exists(plot_dir)) {
  dir.create(plot_dir)
  message("Created directory: ", plot_dir)
}

###############################################################################
# 4) MAIN LOOP: FIT + BUILD TWO OUTPUT TABLES + GENERATE QC PLOTS

summary_rows_core <- list()
summary_rows_full <- list()
k <- 0L
j <- 0L

for (band in bands) {
  
  message("==== BAND: ", band, " ====")
  
  df_band <- COH_WIDE %>%
    dplyr::filter(toupper(as.character(BAND)) == band)
  
  if (nrow(df_band) == 0) {
    warning("No rows found for band: ", band)
    next
  }
  
  for (r in seq_len(nrow(df_band))) {
    
    one_row <- df_band[r, , drop = FALSE]
    id      <- one_row$ID
    
    dat_long <- row_to_long(one_row)
    
    fit <- fit_cosine_multistart(dat_long, b_grid = b_grid, c_grid_len = c_grid_len)
    
    # -------------------------
    # A) CORE SUMMARY ROW
    # -------------------------
    k <- k + 1L
    summary_rows_core[[k]] <- data.frame(
      ID = id,
      BAND = band,
      period_start = one_row$period_start,
      period_end   = one_row$period_end,
      
      COH_MU = ifelse(isTRUE(fit$converged), fit$mu, NA_real_),
      COH_A  = ifelse(isTRUE(fit$converged), fit$a,  NA_real_),
      COH_B  = ifelse(isTRUE(fit$converged), fit$b,  NA_real_),
      COH_C  = ifelse(isTRUE(fit$converged), fit$c,  NA_real_),
      
      # OPTIONAL but helpful
      COH_C_CANON = ifelse(isTRUE(fit$converged), fit$c_canon, NA_real_),
      
      R2  = ifelse(isTRUE(fit$converged), fit$R2,  NA_real_),
      AIC = ifelse(isTRUE(fit$converged), fit$AIC, NA_real_),
      BIC = ifelse(isTRUE(fit$converged), fit$BIC, NA_real_),
      SSE = ifelse(isTRUE(fit$converged), fit$SSE, NA_real_),
      
      converged = isTRUE(fit$converged),
      stringsAsFactors = FALSE
    )
    
    # -------------------------
    # B) FULL OUTPUT ROW
    # -------------------------
    j <- j + 1L
    
    if (!isTRUE(fit$converged) || is.null(fit$best)) {
      
      summary_rows_full[[j]] <- data.frame(
        ID = id,
        BAND = band,
        period_start = one_row$period_start,
        period_end   = one_row$period_end,
        
        mu_est = NA_real_, mu_se = NA_real_, mu_t = NA_real_, mu_p = NA_real_,
        a_est  = NA_real_, a_se  = NA_real_, a_t  = NA_real_, a_p  = NA_real_,
        b_est  = NA_real_, b_se  = NA_real_, b_t  = NA_real_, b_p  = NA_real_,
        c_est  = NA_real_, c_se  = NA_real_, c_t  = NA_real_, c_p  = NA_real_,
        
        sigma = NA_real_,
        SSE   = NA_real_,
        R2    = NA_real_,
        AIC   = NA_real_,
        BIC   = NA_real_,
        logLik= NA_real_,
        
        converged = FALSE,
        stringsAsFactors = FALSE
      )
      
    } else {
      
      full <- extract_nls_full(fit$best, dat_long)
      
      if (is.null(full)) {
        
        summary_rows_full[[j]] <- data.frame(
          ID = id,
          BAND = band,
          period_start = one_row$period_start,
          period_end   = one_row$period_end,
          
          mu_est = fit$mu, mu_se = NA_real_, mu_t = NA_real_, mu_p = NA_real_,
          a_est  = fit$a,  a_se  = NA_real_, a_t  = NA_real_, a_p  = NA_real_,
          b_est  = fit$b,  b_se  = NA_real_, b_t  = NA_real_, b_p  = NA_real_,
          c_est  = fit$c,  c_se  = NA_real_, c_t  = NA_real_, c_p  = NA_real_,
          
          sigma = NA_real_,
          SSE   = fit$SSE,
          R2    = fit$R2,
          AIC   = fit$AIC,
          BIC   = fit$BIC,
          logLik= NA_real_,
          
          converged = TRUE,
          stringsAsFactors = FALSE
        )
        
      } else {
        
        ct <- full$ct
        
        get_row <- function(par) {
          if (!par %in% rownames(ct)) return(c(NA, NA, NA, NA))
          c(ct[par,"Estimate"], ct[par,"Std. Error"], ct[par,"t value"], ct[par,"Pr(>|t|)"])
        }
        
        mu_row <- get_row("mu")
        a_row  <- get_row("a")
        b_row  <- get_row("b")
        c_row  <- get_row("c")
        
        summary_rows_full[[j]] <- data.frame(
          ID = id,
          BAND = band,
          period_start = one_row$period_start,
          period_end   = one_row$period_end,
          
          mu_est = mu_row[1], mu_se = mu_row[2], mu_t = mu_row[3], mu_p = mu_row[4],
          a_est  = a_row[1],  a_se  = a_row[2],  a_t  = a_row[3],  a_p  = a_row[4],
          b_est  = b_row[1],  b_se  = b_row[2],  b_t  = b_row[3],  b_p  = b_row[4],
          c_est  = c_row[1],  c_se  = c_row[2],  c_t  = c_row[3],  c_p  = c_row[4],
          
          sigma = full$sigma,
          SSE   = full$SSE,
          R2    = full$R2,
          AIC   = full$AIC,
          BIC   = full$BIC,
          logLik= full$logLik,
          
          converged = TRUE,
          stringsAsFactors = FALSE
        )
      }
    }
    
    # -------------------------
    # C) QC PLOT: Observed vs Fitted
    # -------------------------
    if (isTRUE(fit$converged) && !is.null(fit$best)) {
      t_obs <- dat_long$t
      y_fitted <- tryCatch(
        predict(fit$best, newdata = data.frame(t = t_obs)),
        error = function(e) fit$mu + fit$a * cos((2*pi/fit$b) * (t_obs - fit$c))
      )
      y_fitted <- pmax(0, pmin(1, y_fitted))  # clip to [0,1]
      
      plot_df <- data.frame(
        t = t_obs,
        y_obs = dat_long$y,
        y_fitted = y_fitted
      )
      
      p <- ggplot(plot_df, aes(x = t)) +
        geom_point(aes(y = y_obs), color = "grey50", alpha = 0.6, size = 0.8) +
        geom_line(aes(y = y_fitted), color = "#D55E00", linewidth = 1.2) +
        labs(
          title = paste0("ID: ", id, " | Band: ", band),
          subtitle = paste0(
            "mu = ", round(fit$mu, 4),
            ", a = ",  round(fit$a, 4),
            ", b = ",  round(fit$b, 2),
            ", c = ",  round(fit$c, 2),
            " (canon ", round(fit$c_canon, 2), ")",
            " | R² = ", round(fit$R2, 3)
          ),
          x = "Time (seconds)",
          y = "Coherence"
        ) +
        ylim(0, 1) +
        theme_minimal() +
        theme(
          plot.title = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, size = 9)
        )
      
      filename <- paste0("COHERENCE_COSINE_FIT_", id, "_", band, ".png")
      filepath <- file.path(plot_dir, filename)
      ggsave(filepath, p, width = 8, height = 5, dpi = 150, bg = "white")
    }
    
    if (k %% 50 == 0) message("...completed fits: ", k)
  }
}

###############################################################################
# 5) BIND + SAVE

DYAD_COHERENCE_COSINE_SUMMARY <- dplyr::bind_rows(summary_rows_core)
DYAD_COHERENCE_COSINE_FULL_OUTPUT <- dplyr::bind_rows(summary_rows_full)

write.csv(DYAD_COHERENCE_COSINE_SUMMARY,
          "DYAD_COHERENCE_COSINE_SUMMARY.csv",
          row.names = FALSE)

write.csv(DYAD_COHERENCE_COSINE_FULL_OUTPUT,
          "DYAD_COHERENCE_COSINE_FULL_OUTPUT.csv",
          row.names = FALSE)

message(
  "\nDONE: saved DYAD_COHERENCE_COSINE_SUMMARY.csv and DYAD_COHERENCE_COSINE_FULL_OUTPUT.csv\n",
  "  Summary rows: ", nrow(DYAD_COHERENCE_COSINE_SUMMARY), "\n",
  "  Full rows:    ", nrow(DYAD_COHERENCE_COSINE_FULL_OUTPUT), "\n",
  "  Converged:    ", sum(DYAD_COHERENCE_COSINE_SUMMARY$converged, na.rm = TRUE), "\n",
  "  QC plots saved to: ", plot_dir, "/"
)
