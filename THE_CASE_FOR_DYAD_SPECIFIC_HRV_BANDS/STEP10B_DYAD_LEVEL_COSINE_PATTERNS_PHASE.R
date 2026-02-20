###############################################################################
# DYAD-SPECIFIC COSINE FITS (PHASE; CIRCULAR) — WITH MU + POSITIVE AMPLITUDE
# - Reads SAMPLE_TRUNC_PHASE_DS_SERIES.csv (wide: one row per ID×BAND)
# - Applies same QC exclusions
# - Drops unnamed columns
# - For each BAND × ID: pivots to long, drops NAs, wraps theta to [-pi, pi]
#
# CIRCULAR COSINE MODEL (WITH MU; 0 still meaningful)
#   theta(t) = wrap_pi( mu + a * cos(((2*pi)/b) * (t - c)) )
#
# Fitting approach:
# - Fit in circular-safe way by minimizing circular SSE:
#     SSE_circ = sum( wrap_pi(theta - theta_hat)^2 )
# - Uses multi-start grid over (mu,a,b,c)
# - Enforces realistic bounds during fitting (L-BFGS-B)
# - Selects best in-bounds by pseudo-AIC
# - Generates QC plots: observed data vs fitted curve for each dyad (ORIGINAL STYLE)
#
# Outputs:
#   1) DYAD_PHASE_CIRCULAR_COSINE_SUMMARY.csv
#        ID, BAND, period_start, period_end,
#        PHS_MU, PHS_A, PHS_B, PHS_C, PHS_C_CANON,
#        R2_circ, AIC_circ, BIC_circ, SSE_circ, converged
#   2) DYAD_PHASE_CIRCULAR_COSINE_FULL_OUTPUT.csv
#        ID, BAND, period_start, period_end,
#        mu/a/b/c (est, approx_se, approx_z, approx_p),
#        sigma_circ, SSE_circ, R2_circ, AIC_circ, BIC_circ, converged
#
# Notes:
# - Pseudo AIC/BIC assumes Gaussian circular residuals; interpret cautiously.
# - Wald p-values from Hessian are approximate.
# - a is constrained to be NON-NEGATIVE to remove sign/phase redundancy.
# - c has periodic redundancy; PHS_C_CANON = c %% b is saved for comparability.
###############################################################################

library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)

###############################################################################
# 0) PATHS + LOAD

work_dir <- "/home/cjh37695/WAVELET_WORKING/TEST_RESP"
setwd(work_dir)

phs_file <- "SAMPLE_TRUNC_PHASE_DS_SERIES.csv"
PHS_WIDE <- read.csv(phs_file, check.names = FALSE, stringsAsFactors = FALSE)

# Drop accidental rowname column
PHS_WIDE <- PHS_WIDE[, -1]

## REMOVE DYADS FAILING QC CHECKS (same list)
BAD_IDS <- c(1035)  # VISUAL INSPECTION
PHS_WIDE$FAIL_QC <- ifelse(PHS_WIDE$ID %in% BAD_IDS, 1, 0)
PHS_WIDE <- subset(PHS_WIDE, FAIL_QC == 0)

# Drop unnamed/blank columns safely
nms <- names(PHS_WIDE)
bad <- is.na(nms) | trimws(nms) == ""
if (any(bad)) {
  message("Dropping ", sum(bad), " columns with NA/blank names.")
  PHS_WIDE <- PHS_WIDE[, !bad, drop = FALSE]
}

###############################################################################
# 1) SETTINGS + BOUNDS DEFINITION

bands  <- c("FOB","RSB","HF","LF")
b_grid <- seq(20, 240, by = 5)  # period grid in seconds

# Multi-start grids
# a constrained >= 0, so use positive multipliers (include 0)
a_multiplier_pos <- c(0, 0.1, 0.25, 0.5, 1, 2)

# mu starts: circular mean +/- offsets (wrapped)
mu_offsets <- c(0, pi/8, -pi/8, pi/4, -pi/4)

# c starts: dense across observed window (0..tmax)
c_grid_len <- 9L

# PARAMETER BOUNDS (phase-specific)
MU_MIN <- -pi
MU_MAX <-  pi
A_MIN  <- 0
A_MAX  <- pi
B_MIN  <- 20
B_MAX  <- 600  # keep as you had; tighten if desired

message("=== PARAMETER BOUNDS ===")
message("mu: [", round(MU_MIN, 3), ", ", round(MU_MAX, 3), "] radians")
message("a:  [", round(A_MIN, 3), ", ", round(A_MAX, 3), "] radians (a >= 0)")
message("b:  [", B_MIN, ", ", B_MAX, "] seconds")
message("c:  [0, max_time] (per dyad)")

###############################################################################
# 2) HELPERS

wrap_pi <- function(x) ((x + pi) %% (2*pi)) - pi

# Wide (one row ID×BAND) -> Long (one row per timepoint)
row_to_long_phase <- function(one_row_df) {
  time_cols <- names(one_row_df)[str_detect(names(one_row_df), "^t[_\\.]")]
  if (length(time_cols) == 0) stop("No time columns found matching ^t_ / ^t.")
  
  one_row_df %>%
    dplyr::select(ID, BAND, period_start, period_end, all_of(time_cols)) %>%
    tidyr::pivot_longer(
      cols      = all_of(time_cols),
      names_to  = "t_name",
      values_to = "theta"
    ) %>%
    dplyr::mutate(
      t     = suppressWarnings(as.numeric(sub("^t[_\\.]", "", t_name))),
      theta = wrap_pi(as.numeric(theta))
    ) %>%
    dplyr::filter(!is.na(t), !is.na(theta))  # drop missingness
}

# Circular mean
circ_mean <- function(theta) atan2(mean(sin(theta), na.rm = TRUE),
                                   mean(cos(theta), na.rm = TRUE))

# Circular cosine prediction (WITH mu)
theta_hat <- function(t, mu, a, b, c) {
  wrap_pi(mu + a * cos(((2*pi)/b) * (t - c)))
}

# Circular residuals and SSE
circ_resid <- function(theta, theta_pred) wrap_pi(theta - theta_pred)
sse_circ <- function(theta, theta_pred) sum(circ_resid(theta, theta_pred)^2, na.rm = TRUE)

# "Circular R2" using SSE on wrapped residuals vs SSE around circular mean
circ_r2 <- function(theta, theta_pred) {
  SSE  <- sse_circ(theta, theta_pred)
  mu0  <- circ_mean(theta)
  SSE0 <- sse_circ(theta, rep(mu0, length(theta)))
  if (!is.finite(SSE0) || SSE0 == 0) return(NA_real_)
  1 - SSE / SSE0
}

# Pseudo AIC/BIC from SSE under Gaussian assumption on circular residuals
# n = number of obs, k = number parameters (mu,a,b,c = 4)
pseudo_aic_bic <- function(SSE, n, k = 4) {
  if (!is.finite(SSE) || SSE <= 0 || n <= 0) return(c(AIC = NA_real_, BIC = NA_real_))
  ll  <- -0.5 * n * (log(2*pi) + 1 + log(SSE/n))  # sigma^2 = SSE/n
  AIC <- -2 * ll + 2 * k
  BIC <- -2 * ll + log(n) * k
  c(AIC = AIC, BIC = BIC)
}

# Validate if parameters are within bounds
is_valid_params_phase <- function(mu, a, b, c, c_max) {
  if (!all(is.finite(c(mu, a, b, c)))) return(FALSE)
  if (mu < MU_MIN || mu > MU_MAX) return(FALSE)
  if (a  < A_MIN  || a  > A_MAX)  return(FALSE)
  if (b  < B_MIN  || b  > B_MAX)  return(FALSE)
  if (c  < 0      || c  > c_max)  return(FALSE)
  return(TRUE)
}

###############################################################################
# 3) DYAD×BAND FITTER (multi-start grid + optim refinement WITH BOUNDS VALIDATION)

fit_phase_circular_multistart <- function(dat, b_grid, c_grid_len = 9L) {
  
  if (nrow(dat) < 10) {
    return(list(converged = FALSE, reason = "too_few_obs"))
  }
  
  tmax <- max(dat$t, na.rm = TRUE)
  if (!is.finite(tmax) || tmax <= 0) {
    return(list(converged = FALSE, reason = "bad_time"))
  }
  
  mu0 <- circ_mean(dat$theta)
  
  # amplitude starts from SD(theta), constrained to [0, pi]
  sdy <- sd(dat$theta, na.rm = TRUE)
  if (!is.finite(sdy) || sdy == 0) sdy <- 0.5
  a_grid <- pmin(A_MAX, pmax(A_MIN, abs(sdy) * a_multiplier_pos))
  
  mu_grid <- unique(wrap_pi(mu0 + mu_offsets))
  c_grid  <- unique(round(seq(0, tmax, length.out = c_grid_len), 6))
  
  start_grid <- expand.grid(
    mu = mu_grid,
    a  = a_grid,
    b  = b_grid,
    c  = c_grid,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  
  valid_fits <- list()
  
  obj_fun <- function(par, t, theta) {
    mu <- par[1]; a <- par[2]; b <- par[3]; c <- par[4]
    if (!is.finite(b) || b <= 0) return(Inf)
    if (!is.finite(c) || c < 0 || c > max(t, na.rm = TRUE)) return(Inf)
    if (!is.finite(mu) || mu < MU_MIN || mu > MU_MAX) return(Inf)
    if (!is.finite(a)  || a  < A_MIN  || a  > A_MAX)  return(Inf)
    pred <- theta_hat(t, mu, a, b, c)
    sse_circ(theta, pred)
  }
  
  lower <- c(mu = MU_MIN, a = A_MIN, b = B_MIN, c = 0)
  upper <- c(mu = MU_MAX, a = A_MAX, b = B_MAX, c = tmax)
  
  for (i in seq_len(nrow(start_grid))) {
    
    st <- c(
      mu = start_grid$mu[i],
      a  = start_grid$a[i],
      b  = start_grid$b[i],
      c  = start_grid$c[i]
    )
    
    opt <- tryCatch(
      optim(
        par     = st,
        fn      = obj_fun,
        t       = dat$t,
        theta   = dat$theta,
        method  = "L-BFGS-B",
        lower   = lower,
        upper   = upper,
        control = list(maxit = 200)
      ),
      error = function(e) NULL
    )
    
    if (!is.null(opt) && is.finite(opt$value)) {
      
      mu_hat <- unname(opt$par[1])
      a_hat  <- unname(opt$par[2])
      b_hat  <- unname(opt$par[3])
      c_hat  <- unname(opt$par[4])
      
      if (is_valid_params_phase(mu_hat, a_hat, b_hat, c_hat, tmax)) {
        
        pred <- theta_hat(dat$t, mu_hat, a_hat, b_hat, c_hat)
        SSE  <- sse_circ(dat$theta, pred)
        R2   <- circ_r2(dat$theta, pred)
        n    <- nrow(dat)
        ab   <- pseudo_aic_bic(SSE, n, k = 4)
        
        valid_fits[[length(valid_fits) + 1]] <- list(
          opt = opt,
          mu = mu_hat,
          a  = a_hat,
          b  = b_hat,
          c  = c_hat,
          SSE = SSE,
          R2  = R2,
          AIC = ab["AIC"],
          BIC = ab["BIC"]
        )
      }
    }
  }
  
  if (length(valid_fits) == 0) {
    return(list(converged = FALSE, reason = "no_valid_in_bounds_fit"))
  }
  
  best_idx <- which.min(sapply(valid_fits, `[[`, "AIC"))
  best_fit_obj <- valid_fits[[best_idx]]
  best_opt <- best_fit_obj$opt
  
  mu_hat <- best_fit_obj$mu
  a_hat  <- best_fit_obj$a
  b_hat  <- best_fit_obj$b
  c_hat  <- best_fit_obj$c
  
  pred <- theta_hat(dat$t, mu_hat, a_hat, b_hat, c_hat)
  SSE  <- best_fit_obj$SSE
  R2   <- best_fit_obj$R2
  
  n  <- nrow(dat)
  ab <- pseudo_aic_bic(SSE, n, k = 4)
  
  c_canon <- if (is.finite(b_hat) && b_hat > 0) (c_hat %% b_hat) else NA_real_
  
  # Approx SEs / z / p from Hessian at optimum (optimHess)
  se <- c(mu = NA_real_, a = NA_real_, b = NA_real_, c = NA_real_)
  z  <- c(mu = NA_real_, a = NA_real_, b = NA_real_, c = NA_real_)
  p  <- c(mu = NA_real_, a = NA_real_, b = NA_real_, c = NA_real_)
  
  H <- tryCatch(
    optimHess(
      par   = best_opt$par,
      fn    = obj_fun,
      t     = dat$t,
      theta = dat$theta
    ),
    error = function(e) NULL
  )
  
  if (!is.null(H)) {
    V <- tryCatch(solve(H), error = function(e) NULL)
    if (!is.null(V)) {
      se_vals <- sqrt(pmax(diag(V), 0))
      if (length(se_vals) == 4) {
        names(se_vals) <- c("mu","a","b","c")
        se <- se_vals
        z_vals <- best_opt$par / se
        p_vals <- 2 * (1 - pnorm(abs(z_vals)))
        names(z_vals) <- names(p_vals) <- c("mu","a","b","c")
        z <- z_vals
        p <- p_vals
      }
    }
  }
  
  list(
    converged = TRUE,
    mu    = mu_hat,
    a     = a_hat,
    b     = b_hat,
    c     = c_hat,
    c_canon = c_canon,
    SSE   = SSE,
    R2    = R2,
    AIC   = ab["AIC"],
    BIC   = ab["BIC"],
    se    = se,
    z     = z,
    p     = p,
    sigma = sqrt(SSE / n),
    opt   = best_opt
  )
}

###############################################################################
# 4) CREATE PLOT DIRECTORY

plot_dir <- "PHASE_CIRCULAR_COSINE_FIT_PLOTS"
if (!dir.exists(plot_dir)) {
  dir.create(plot_dir)
  message("Created directory: ", plot_dir)
}

###############################################################################
# 5) MAIN LOOP: FIT + BUILD TWO OUTPUT TABLES + GENERATE QC PLOTS

summary_rows_core <- list()
summary_rows_full <- list()
k <- 0L
j <- 0L

for (band in bands) {
  
  message("==== PHASE BAND: ", band, " ====")
  
  df_band <- PHS_WIDE %>%
    dplyr::filter(toupper(as.character(BAND)) == band)
  
  if (nrow(df_band) == 0) {
    warning("No rows found for band: ", band)
    next
  }
  
  for (r in seq_len(nrow(df_band))) {
    
    one_row <- df_band[r, , drop = FALSE]
    id      <- one_row$ID
    
    dat_long <- row_to_long_phase(one_row)
    
    fit <- fit_phase_circular_multistart(dat_long, b_grid = b_grid, c_grid_len = c_grid_len)
    
    # -------------------------
    # A) CORE SUMMARY ROW
    # -------------------------
    k <- k + 1L
    summary_rows_core[[k]] <- data.frame(
      ID = id,
      BAND = band,
      period_start = one_row$period_start,
      period_end   = one_row$period_end,
      
      PHS_MU = ifelse(isTRUE(fit$converged), fit$mu, NA_real_),
      PHS_A  = ifelse(isTRUE(fit$converged), fit$a,  NA_real_),
      PHS_B  = ifelse(isTRUE(fit$converged), fit$b,  NA_real_),
      PHS_C  = ifelse(isTRUE(fit$converged), fit$c,  NA_real_),
      PHS_C_CANON = ifelse(isTRUE(fit$converged), fit$c_canon, NA_real_),
      
      R2_circ  = ifelse(isTRUE(fit$converged), fit$R2,  NA_real_),
      AIC_circ = ifelse(isTRUE(fit$converged), fit$AIC, NA_real_),
      BIC_circ = ifelse(isTRUE(fit$converged), fit$BIC, NA_real_),
      SSE_circ = ifelse(isTRUE(fit$converged), fit$SSE, NA_real_),
      
      converged = isTRUE(fit$converged),
      stringsAsFactors = FALSE
    )
    
    # -------------------------
    # B) FULL OUTPUT ROW
    # -------------------------
    j <- j + 1L
    
    if (!isTRUE(fit$converged)) {
      
      summary_rows_full[[j]] <- data.frame(
        ID = id,
        BAND = band,
        period_start = one_row$period_start,
        period_end   = one_row$period_end,
        
        mu_est = NA_real_, mu_se = NA_real_, mu_z = NA_real_, mu_p = NA_real_,
        a_est  = NA_real_, a_se  = NA_real_, a_z  = NA_real_, a_p  = NA_real_,
        b_est  = NA_real_, b_se  = NA_real_, b_z  = NA_real_, b_p  = NA_real_,
        c_est  = NA_real_, c_se  = NA_real_, c_z  = NA_real_, c_p  = NA_real_,
        
        sigma_circ = NA_real_,
        SSE_circ   = NA_real_,
        R2_circ    = NA_real_,
        AIC_circ   = NA_real_,
        BIC_circ   = NA_real_,
        
        converged = FALSE,
        stringsAsFactors = FALSE
      )
      
    } else {
      
      summary_rows_full[[j]] <- data.frame(
        ID = id,
        BAND = band,
        period_start = one_row$period_start,
        period_end   = one_row$period_end,
        
        mu_est = fit$mu,           mu_se = fit$se["mu"], mu_z = fit$z["mu"], mu_p = fit$p["mu"],
        a_est  = fit$a,            a_se  = fit$se["a"],  a_z  = fit$z["a"],  a_p  = fit$p["a"],
        b_est  = fit$b,            b_se  = fit$se["b"],  b_z  = fit$z["b"],  b_p  = fit$p["b"],
        c_est  = fit$c,            c_se  = fit$se["c"],  c_z  = fit$z["c"],  c_p  = fit$p["c"],
        
        sigma_circ = fit$sigma,
        SSE_circ   = fit$SSE,
        R2_circ    = fit$R2,
        AIC_circ   = fit$AIC,
        BIC_circ   = fit$BIC,
        
        converged = TRUE,
        stringsAsFactors = FALSE
      )
    }
    
    # -------------------------
    # C) QC PLOT: Observed vs Fitted (ORIGINAL PLOTTING STYLE)
    # -------------------------
    if (isTRUE(fit$converged)) {
      
      t_obs <- dat_long$t
      theta_fitted <- theta_hat(t_obs, fit$mu, fit$a, fit$b, fit$c)
      
      plot_df <- data.frame(
        t = t_obs,
        theta_obs = dat_long$theta,
        theta_fitted = theta_fitted
      )
      
      p <- ggplot(plot_df, aes(x = t)) +
        geom_point(aes(y = theta_obs), color = "grey50", alpha = 0.6, size = 0.8) +
        geom_line(aes(y = theta_fitted), color = "#0072B2", linewidth = 1.2) +
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
          y = "Phase (radians)"
        ) +
        ylim(-pi, pi) +
        theme_minimal() +
        theme(
          plot.title = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, size = 9)
        )
      
      filename <- paste0("PHASE_CIRCULAR_COSINE_FIT_", id, "_", band, ".png")
      filepath <- file.path(plot_dir, filename)
      ggsave(filepath, p, width = 8, height = 5, dpi = 150, bg = "white")
    }
    
    if (k %% 50 == 0) message("...completed phase fits: ", k)
  }
}

###############################################################################
# 6) BIND + SAVE

DYAD_PHASE_CIRCULAR_COSINE_SUMMARY <- dplyr::bind_rows(summary_rows_core)
DYAD_PHASE_CIRCULAR_COSINE_FULL_OUTPUT <- dplyr::bind_rows(summary_rows_full)

write.csv(DYAD_PHASE_CIRCULAR_COSINE_SUMMARY,
          "DYAD_PHASE_CIRCULAR_COSINE_SUMMARY.csv",
          row.names = FALSE)

write.csv(DYAD_PHASE_CIRCULAR_COSINE_FULL_OUTPUT,
          "DYAD_PHASE_CIRCULAR_COSINE_FULL_OUTPUT.csv",
          row.names = FALSE)

message(
  "\nDONE: saved DYAD_PHASE_CIRCULAR_COSINE_SUMMARY.csv and DYAD_PHASE_CIRCULAR_COSINE_FULL_OUTPUT.csv\n",
  "  Summary rows: ", nrow(DYAD_PHASE_CIRCULAR_COSINE_SUMMARY), "\n",
  "  Full rows:    ", nrow(DYAD_PHASE_CIRCULAR_COSINE_FULL_OUTPUT), "\n",
  "  Converged:    ", sum(DYAD_PHASE_CIRCULAR_COSINE_SUMMARY$converged, na.rm = TRUE), "\n",
  "  QC plots saved to: ", plot_dir, "/"
)
