#Library

library(dplyr)
library(psych)
library(stringr)
library(purrr)
library(MplusAutomation)

########## SET WORKING DIRECTORY, OUTPUT DIRECTORY, AND FILE NAME ########## 

setwd("/home/cjh37695/ABCD_PROJECTS/")

OUT_DIR <- file.path(getwd(), "SUICIDALITY_SOCIAL_FACTORS")

# create name with analytic step 
FILENAME <- "ABCD_SSP_"

# get today's date appended
DATE_STAMP <- format(Sys.Date(), "%m.%d.%y")

########## LOAD IN BASE_DATA df ########## 

df <- read.csv(
  file.path(OUT_DIR, "ABCD_SSP_BASE_DATA_05.21.26.csv")
)

codebook <- read.csv(
  file.path(OUT_DIR, "ABCD_SSP_CODEBOOK_05.21.26.csv")
)

################################################################################
##################### COMPUTE WAVE-SPECIFIC PAST/PRESENT SUICIDE INDICATORS

# This "total" approach is used because we are interested broadly in if these symptoms are 
# present for a person by that wave or not. This is adopted because there are 
# idiosyncrasies for people between waves that make us weary to adopt dynamic 
# wave-specific scores as accurate indices of waxing and waning suicide risk. 

# helper function: Create a broad dummy from two component variables.
# Each broad dummy = 1 if either present-year or past-year item is endorsed.
# If both component items are missing, broad dummy = NA; otherwise = 0.

past_pres_SuRrk_dummy <- function(data, new_var, var1, var2) {
  
  requested_vars <- c(var1, var2)
  available_vars <- requested_vars[requested_vars %in% names(data)]
  missing_vars   <- setdiff(requested_vars, available_vars)
  
  if (length(available_vars) == 0) {
    warning("Neither ", var1, " nor ", var2, " found. ", new_var, " not created.")
    return(data)
  }
  
  vals <- data[available_vars]
  
  n_full_na <- sum(rowSums(!is.na(vals)) == 0)
  n_usable  <- sum(rowSums(!is.na(vals)) > 0)
  
  data[[new_var]] <- ifelse(
    rowSums(vals == 1, na.rm = TRUE) > 0,
    1,
    ifelse(
      rowSums(!is.na(vals)) == 0,
      NA,
      0
    )
  )
  
  cat("\nCreated:", new_var, "\n")
  cat("  Found variables:", paste(available_vars, collapse = ", "), "\n")
  
  if (length(missing_vars) > 0) {
    cat("  Missing variables:", paste(missing_vars, collapse = ", "), "\n")
  } else {
    cat("  Missing variables: none\n")
  }
  
  cat("  Rows fully NA across available input variable(s):", n_full_na, "\n")
  cat("  Rows with at least one non-missing input value:", n_usable, "\n")
  
  return(data)
}

###### applied to passive suicidal idealizations 
#GCDS version: SuIdPas[pr|pt]_[W], ABCD version: mh_y_ksads__suic__pass__[pres|past]_dx

df <- past_pres_SuRrk_dummy(df, "SuIdPasT_1", "SuIdPaspr_1", "SuIdPaspt_1")
df <- past_pres_SuRrk_dummy(df, "SuIdPasT_3", "SuIdPaspr_3", "SuIdPaspt_3")
df <- past_pres_SuRrk_dummy(df, "SuIdPasT_5", "SuIdPaspr_5", "SuIdPaspt_5")
df <- past_pres_SuRrk_dummy(df, "SuIdPasT_7", "SuIdPaspr_7", "SuIdPaspt_7")
df <- past_pres_SuRrk_dummy(df, "SuIdPasT_9", "SuIdPaspr_9", "SuIdPaspt_9")
df <- past_pres_SuRrk_dummy(df, "SuIdPasT_11", "SuIdPaspr_11", "SuIdPaspt_11")
df <- past_pres_SuRrk_dummy(df, "SuIdPasT_13", "SuIdPaspr_13", "SuIdPaspt_13")

###### applied to suicide attempts 
#GCDS version: SuAtt[pr|pt]_[W], ABCD version:  mh_y_ksads__suic__atmpt__[pres|past]_dx 

df <- past_pres_SuRrk_dummy(df, "SuAttT_1", "SuAttpr_1", "SuAttpt_1")
df <- past_pres_SuRrk_dummy(df, "SuAttT_3", "SuAttpr_3", "SuAttpt_3")
df <- past_pres_SuRrk_dummy(df, "SuAttT_5", "SuAttpr_5", "SuAttpt_5")
df <- past_pres_SuRrk_dummy(df, "SuAttT_7", "SuAttpr_7", "SuAttpt_7")
df <- past_pres_SuRrk_dummy(df, "SuAttT_9", "SuAttpr_9", "SuAttpt_9")
df <- past_pres_SuRrk_dummy(df, "SuAttT_11", "SuAttpr_11", "SuAttpt_11")
df <- past_pres_SuRrk_dummy(df, "SuAttT_13", "SuAttpr_13", "SuAttpt_13")

###### applied to suicidal ideations: active methods
#GCDS version: SuIdMTD[pr|pt]_[W], ABCD version:  mh_y_ksads__suic__actv__mthd__[pres|past]_dx 

df <- past_pres_SuRrk_dummy(df, "SuIdMtdT_1", "SuIdMTDpr_1", "SuIdMTDpt_1")
df <- past_pres_SuRrk_dummy(df, "SuIdMtdT_3", "SuIdMTDpr_3", "SuIdMTDpt_3")
df <- past_pres_SuRrk_dummy(df, "SuIdMtdT_5", "SuIdMTDpr_5", "SuIdMTDpt_5")
df <- past_pres_SuRrk_dummy(df, "SuIdMtdT_7", "SuIdMTDpr_7", "SuIdMTDpt_7")
df <- past_pres_SuRrk_dummy(df, "SuIdMtdT_9", "SuIdMTDpr_9", "SuIdMTDpt_9")
df <- past_pres_SuRrk_dummy(df, "SuIdMtdT_11", "SuIdMTDpr_11", "SuIdMTDpt_11")
df <- past_pres_SuRrk_dummy(df, "SuIdMtdT_13", "SuIdMTDpr_13", "SuIdMTDpt_13")


###### applied to suicidal ideations: active nonspecific
#GCDS version: SuIdNonS[pr|pt]_[W], ABCD version: mh_y_ksads__suic__actv__[pres|past]_dx 

df <- past_pres_SuRrk_dummy(df, "SuIdNonST_1", "SuIdNonSpr_1", "SuIdNonSpt_1")
df <- past_pres_SuRrk_dummy(df, "SuIdNonST_3", "SuIdNonSpr_3", "SuIdNonSpt_3")
df <- past_pres_SuRrk_dummy(df, "SuIdNonST_5", "SuIdNonSpr_5", "SuIdNonSpt_5")
df <- past_pres_SuRrk_dummy(df, "SuIdNonST_7", "SuIdNonSpr_7", "SuIdNonSpt_7")
df <- past_pres_SuRrk_dummy(df, "SuIdNonST_9", "SuIdNonSpr_9", "SuIdNonSpt_9")
df <- past_pres_SuRrk_dummy(df, "SuIdNonST_11", "SuIdNonSpr_11", "SuIdNonSpt_11")
df <- past_pres_SuRrk_dummy(df, "SuIdNonST_13", "SuIdNonSpr_13", "SuIdNonSpt_13")


##----------- remove the wave-specific past- and present-dignosis variables 

df <- df %>%
  dplyr::select(
    -dplyr::matches(
      "^(SuAttpt|SuAttpr|SuIdNonSpt|SuIdMtdpr|SuIdMtdpt|SuIdNonSpr|SuIdActpt|SuIdActpr|SuIdPaspt|SuIdPaspr)_\\d+$"
    )
  )

#check it 
names(df)

################################################################################
####### COMPUTE EVENT-HISTORY SUICIDE-RISK VARIABLES FOR SURVIVAL ANALYSIS

# Recode wave-specific suicide outcomes into first-onset event-history indicators.
# Participants remain in the risk set until their first observed event. 
# later waves are set to NA after onset.
# This format is required for discrete-time survival models estimating 
# wave-specific first-event hazard.
## People who experienced the event "drop out"

################################################################################
##################### COMPUTE EVENT-HISTORY VARIABLES FOR SURVIVAL ANALYSIS

# Recode wave-specific event indicators into first-onset event-history indicators.
# Participants remain in the risk set until first observed onset; later waves are set to NA.
# This supports discrete-time survival models estimating wave-specific first-event hazard.

survival_event_coding <- function(data, base_var, waves, prefix = "SV_") {
  
  # Construct input and output variable names
  input_vars  <- paste0(base_var, "_", waves)
  output_vars <- paste0(prefix, base_var, "_", waves)
  
  # Keep only waves that are actually present in the data
  available <- input_vars %in% names(data)
  
  if (!any(available)) {
    warning("No variables found for ", base_var, ". No survival variables created.")
    return(data)
  }
  
  input_vars  <- input_vars[available]
  output_vars <- output_vars[available]
  waves_found <- waves[available]
  waves_miss  <- waves[!available]
  
  # Copy original event indicators into survival-coded variables
  data[output_vars] <- data[input_vars]
  
  # Set later waves to NA after first observed event
  for (i in seq_along(output_vars)) {
    
    if (i == 1) next
    
    earlier_vars <- output_vars[1:(i - 1)]
    
    prior_event <- rowSums(data[earlier_vars] == 1, na.rm = TRUE) > 0
    
    data[[output_vars[i]]][prior_event] <- NA
  }
  
  # Print diagnostic summary
  cat("\nCreated survival-coded variables for:", base_var, "\n")
  cat("  Waves found:", paste0("_", waves_found, collapse = ", "), "\n")
  
  if (length(waves_miss) > 0) {
    cat("  Waves missing:", paste0("_", waves_miss, collapse = ", "), "\n")
  } else {
    cat("  Waves missing: none\n")
  }
  
  cat("  New variables:", paste(output_vars, collapse = ", "), "\n")
  
  cat("  Post-coding event counts:\n")
  print(lapply(data[output_vars], table, useNA = "ifany"))
  
  return(data)
}

######  APPLIED  

# indicate the wave suffixes to be included 
surv_waves <- c(1, 3, 5, 7, 9, 11, 13)

# Passive suicidal ideation
df <- survival_event_coding(
  data = df,
  base_var = "SuIdPasT",
  waves = surv_waves
)

# Suicide attempts
df <- survival_event_coding(
  data = df,
  base_var = "SuAttT",
  waves = surv_waves
)

# Active suicidal ideation with method
df <- survival_event_coding(
  data = df,
  base_var = "SuIdMtdT",
  waves = surv_waves
)

# Active nonspecific suicidal ideation
df <- survival_event_coding(
  data = df,
  base_var = "SuIdNonST",
  waves = surv_waves
)


################################################################################
##################### SUMMARIZE SUICIDE-RISK PREVALENCE AND FIRST-ONSET EVENTS

# This table summarizes each suicide-risk indicator by wave:
# 1) Original wave-specific prevalence: number and percent endorsing the original variable.
# 2) First-onset prevalence: number and percent endorsing the SV_ event-history variable.
# 3) Wave-specific chi-square tests by Y_SEX for both original and first-onset indicators.

summarize_survival_indicator <- function(data, base_var, waves, sex_var = "Y_SEX") {
  
  out <- lapply(waves, function(w) {
    
    orig_var <- paste0(base_var, "_", w)
    sv_var   <- paste0("SV_", base_var, "_", w)
    
    # Skip wave if neither original nor SV variable exists
    if (!orig_var %in% names(data) & !sv_var %in% names(data)) {
      warning("Neither ", orig_var, " nor ", sv_var, " found. Skipping wave ", w, ".")
      return(NULL)
    }
    
    #### Original wave-specific prevalence
    
    if (orig_var %in% names(data)) {
      
      orig_n_obs <- sum(!is.na(data[[orig_var]]))
      orig_n_1   <- sum(data[[orig_var]] == 1, na.rm = TRUE)
      orig_pct   <- ifelse(orig_n_obs > 0, 100 * orig_n_1 / orig_n_obs, NA_real_)
      
      orig_chisq_p <- tryCatch({
        tab <- table(data[[orig_var]], data[[sex_var]], useNA = "no")
        if (all(dim(tab) >= c(2, 2))) {
          suppressWarnings(chisq.test(tab)$p.value)
        } else {
          NA_real_
        }
      }, error = function(e) NA_real_)
      
    } else {
      
      orig_n_obs <- NA_integer_
      orig_n_1   <- NA_integer_
      orig_pct   <- NA_real_
      orig_chisq_p <- NA_real_
    }
    
    
    #### First-onset / event-history prevalence
    
    if (sv_var %in% names(data)) {
      
      sv_n_obs <- sum(!is.na(data[[sv_var]]))
      sv_n_1   <- sum(data[[sv_var]] == 1, na.rm = TRUE)
      sv_pct   <- ifelse(sv_n_obs > 0, 100 * sv_n_1 / sv_n_obs, NA_real_)
      
      sv_chisq_p <- tryCatch({
        tab <- table(data[[sv_var]], data[[sex_var]], useNA = "no")
        if (all(dim(tab) >= c(2, 2))) {
          suppressWarnings(chisq.test(tab)$p.value)
        } else {
          NA_real_
        }
      }, error = function(e) NA_real_)
      
    } else {
      
      sv_n_obs <- NA_integer_
      sv_n_1   <- NA_integer_
      sv_pct   <- NA_real_
      sv_chisq_p <- NA_real_
    }
    
    
    #### Return one row
    
    data.frame(
      INDICATOR = base_var,
      WAVE = paste0("_", w),
      
      ORIGINAL_VAR = ifelse(orig_var %in% names(data), orig_var, NA),
      ORIGINAL_N_OBSERVED = orig_n_obs,
      ORIGINAL_N_EVENT = orig_n_1,
      ORIGINAL_PCT = round(orig_pct, 2),
      ORIGINAL_CHISQ_BY_SEX_P = round(orig_chisq_p, 4),
      
      SV_VAR = ifelse(sv_var %in% names(data), sv_var, NA),
      SV_RISKSET_N_OBSERVED = sv_n_obs,
      SV_NEW_ONSET_N = sv_n_1,
      SV_NEW_ONSET_PCT = round(sv_pct, 2),
      SV_CHISQ_BY_SEX_P = round(sv_chisq_p, 4),
      
      stringsAsFactors = FALSE
    )
  })
  
  dplyr::bind_rows(out)
}


################################################################################
##################### APPLY SUMMARY FUNCTION

surv_waves <- c(1, 3, 5, 7, 9, 11, 13) ## the final waves have a suscpicious missing pattern - don't use

SUICIDE_RISK_SUMMARY <- dplyr::bind_rows(
  
  summarize_survival_indicator(
    data = df,
    base_var = "SuIdPasT",
    waves = surv_waves
  ),
  
  summarize_survival_indicator(
    data = df,
    base_var = "SuAttT",
    waves = surv_waves
  ),
  
  summarize_survival_indicator(
    data = df,
    base_var = "SuIdMtdT",
    waves = surv_waves
  ),
  
  summarize_survival_indicator(
    data = df,
    base_var = "SuIdNonST",
    waves = surv_waves
  )
)

SUICIDE_RISK_SUMMARY


##################### CUMULATIVE FIRST-ONSET SUMMARY

CUMULATIVE_ONSET_SUMMARY <- SUICIDE_RISK_SUMMARY %>%
  dplyr::group_by(INDICATOR) %>%
  dplyr::summarise(
    BASELINE_N = dplyr::first(ORIGINAL_N_OBSERVED),
    CUMULATIVE_NEW_ONSET_N = sum(SV_NEW_ONSET_N, na.rm = TRUE),
    CUMULATIVE_NEW_ONSET_PCT_BASELINE = round(
      100 * CUMULATIVE_NEW_ONSET_N / 11775, 2
    ),
    .groups = "drop"
  )

CUMULATIVE_ONSET_SUMMARY


################################################################################
##################### REMOVE ORIGINAL WAVE-SPECIFIC SUICIDE-RISK VARIABLES

# After creating first-onset event-history variables and verifying summary output,
# remove the original repeated wave-specific suicide-risk indicators.
# The SV_ variables are retained as the analysis-ready DTS/event-history outcomes.

df_RED <- df %>%
  dplyr::select(
    -dplyr::matches("^(SuIdPasT|SuAttT|SuIdMtdT|SuIdNonST)_\\d+$")
  )


################################################################################
###################### SCALE COMPUTATION AND LONGITUDINAL CHANGE VISUALIZATION

make_scale_means <- function(data, item_sets, new_names = names(item_sets),
                             min_completion_prop = 0.80, digits = 3,
                             min_items_for_omega = 4,
                             min_complete_for_omega = 200) {

  if (!requireNamespace("psych", quietly = TRUE)) {
    stop("Package 'psych' required. Install with install.packages('psych')")
  }

  if (is.null(new_names)) {
    stop("Please name item_sets or provide new_names.")
  }

  if (length(item_sets) != length(new_names)) {
    stop("item_sets and new_names must have the same length.")
  }

  df <- data
  desc_list <- vector("list", length(item_sets))

  for (i in seq_along(item_sets)) {

    vars <- item_sets[[i]]
    new_var <- new_names[i]

    missing_vars <- setdiff(vars, names(df))
    if (length(missing_vars) > 0) {
      stop(
        paste0("Scale '", new_var, "' contains variables not found in data: ",
               paste(missing_vars, collapse = ", "))
      )
    }

    scale_df <- df[, vars, drop = FALSE]
    n_items <- length(vars)

    # Completion / missingness for row scoring
    n_missing_items <- rowSums(is.na(scale_df))
    n_observed_items <- n_items - n_missing_items
    prop_observed <- n_observed_items / n_items

    all_missing <- n_missing_items == n_items
    fails_completion_rule <- prop_observed < min_completion_prop & !all_missing

    # Compute scale mean
    scale_mean <- rowMeans(scale_df, na.rm = TRUE)
    scale_mean[all_missing] <- NA
    scale_mean[fails_completion_rule] <- NA
    scale_mean[is.nan(scale_mean)] <- NA

    df[[new_var]] <- scale_mean

    # Descriptive Ns
    N_valid <- sum(!is.na(scale_mean))
    N_missing_total <- sum(is.na(scale_mean))
    N_missing_all_items <- sum(all_missing)
    N_missing_completion_rule <- sum(fails_completion_rule)

    # Alpha
    alpha_val <- tryCatch({
      suppressWarnings(
        psych::alpha(scale_df, warnings = FALSE, check.keys = FALSE)$total$raw_alpha
      )
    }, error = function(e) NA)

    # Omega
    complete_cases <- stats::complete.cases(scale_df)
    n_complete <- sum(complete_cases)
    omega_val <- NA
    omega_note <- ""

    if (n_items < min_items_for_omega) {
      omega_note <- paste0("omega not computed: fewer than ", min_items_for_omega, " items")
    } else if (n_complete < min_complete_for_omega) {
      omega_note <- paste0("omega not computed: fewer than ", min_complete_for_omega,
                           " complete cases (n=", n_complete, ")")
    } else {
      omega_attempt <- tryCatch({
        suppressWarnings(
          psych::omega(scale_df[complete_cases, , drop = FALSE],
                       nfactors = 1,
                       plot = FALSE,
                       warnings = FALSE)
        )
      }, error = function(e) NULL)

      if (is.null(omega_attempt)) {
        omega_note <- "omega failed"
        omega_val <- NA
      } else {
        omega_val <- omega_attempt$omega.tot
        if (is.null(omega_val) || length(omega_val) == 0 || !is.finite(omega_val)) {
          omega_val <- NA
          omega_note <- "omega non-finite / unstable"
        } else {
          omega_note <- "ok"
        }
      }
    }

    # Reliability flag
    rel_vals <- c(alpha_val, omega_val)
    rel_vals <- rel_vals[!is.na(rel_vals)]

    rel_flag <- if (length(rel_vals) == 0) {
      ""
    } else if (all(rel_vals < 0.70)) {
      "X"
    } else {
      ""
    }

    # Safe descriptives
    safe_mean <- if (all(is.na(scale_mean))) NA else round(mean(scale_mean, na.rm = TRUE), digits)
    safe_sd   <- if (all(is.na(scale_mean))) NA else round(sd(scale_mean, na.rm = TRUE), digits)
    safe_min  <- if (all(is.na(scale_mean))) NA else round(min(scale_mean, na.rm = TRUE), digits)
    safe_max  <- if (all(is.na(scale_mean))) NA else round(max(scale_mean, na.rm = TRUE), digits)

    desc_list[[i]] <- data.frame(
      variable = new_var,
      mean = safe_mean,
      sd = safe_sd,
      min = safe_min,
      max = safe_max,
      n_items = n_items,
      alpha = ifelse(is.na(alpha_val), NA, round(alpha_val, digits)),
      omega = ifelse(is.na(omega_val), NA, round(omega_val, digits)),
      rel_flag = rel_flag,
      omega_note = omega_note,
      N_valid = N_valid,
      N_missing_total = N_missing_total,
      N_missing_all_items = N_missing_all_items,
      N_missing_completion_rule = N_missing_completion_rule,
      n_complete_for_omega = n_complete,
      stringsAsFactors = FALSE
    )
  }

  descriptives_table <- do.call(rbind, desc_list)

  return(list(
    data = df,
    descriptives = descriptives_table
  ))
}


## compute scale means 

names(df_RED)

item_sets <- list(
## adversity predictors
  mFconY_1 = c("Fcon1Y_1","Fcon2rY_1","Fcon3Y_1","Fcon4rY_1","Fcon5Y_1","Fcon6Y_1",
               "Fcon7rY_1", "Fcon8Y_1", "Fcon9Y_1"),
  mFconY_3 = c("Fcon1Y_3","Fcon2rY_3","Fcon3Y_3","Fcon4rY_3","Fcon5Y_3","Fcon6Y_3",
               "Fcon7rY_3", "Fcon8Y_3", "Fcon9Y_3"),
  mFconY_5 = c("Fcon1Y_5","Fcon2rY_5","Fcon3Y_5","Fcon4rY_5","Fcon5Y_5","Fcon6Y_5",
               "Fcon7rY_5", "Fcon8Y_5", "Fcon9Y_5"),
  mFconY_7 = c("Fcon1Y_7","Fcon2rY_7","Fcon3Y_7","Fcon4rY_7","Fcon5Y_7","Fcon6Y_7",
               "Fcon7rY_7", "Fcon8Y_7", "Fcon9Y_7"),
  mFconY_9 = c("Fcon1Y_9","Fcon2rY_9","Fcon3Y_9","Fcon4rY_9","Fcon5Y_9","Fcon6Y_9",
               "Fcon7rY_9", "Fcon8Y_9", "Fcon9Y_9"),
  mFconY_11 = c("Fcon1Y_11","Fcon2rY_11","Fcon3Y_11","Fcon4rY_11","Fcon5Y_11","Fcon6Y_11",
                "Fcon7rY_11", "Fcon8Y_11", "Fcon9Y_11"),
  mFconY_13 = c("Fcon1Y_13","Fcon2rY_13","Fcon3Y_13","Fcon4rY_13","Fcon5Y_13","Fcon6Y_13",
                "Fcon7rY_13", "Fcon8Y_13", "Fcon9Y_13"),
  OvtVic_5 = c("PEQ1_OV_5", "PEQ2_OV_5", "PEQ3_OV_5"),
  OvtVic_7 = c("PEQ1_OV_7", "PEQ2_OV_7", "PEQ3_OV_7"),
  OvtVic_9 = c("PEQ1_OV_9", "PEQ2_OV_9", "PEQ3_OV_9"),
  OvtVic_11 = c("PEQ1_OV_11", "PEQ2_OV_11", "PEQ3_OV_11"),
  OvtVic_13 = c("PEQ1_OV_13", "PEQ2_OV_13", "PEQ3_OV_13"),
  RelVic_5 = c("PEQ1_RLV_5", "PEQ2_RLV_5", "PEQ3_RLV_5"),
  RelVic_7 = c("PEQ1_RLV_7", "PEQ2_RLV_7", "PEQ3_RLV_7"),
  RelVic_9 = c("PEQ1_RLV_9", "PEQ2_RLV_9", "PEQ3_RLV_9"),
  RelVic_11 = c("PEQ1_RLV_11", "PEQ2_RLV_11", "PEQ3_RLV_11"),
  RelVic_13 = c("PEQ1_RLV_13", "PEQ2_RLV_13", "PEQ3_RLV_13"),
  RepVic_5 = c("PEQ1_RPV_5", "PEQ2_RPV_5", "PEQ3_RPV_5"),
  RepVic_7 = c("PEQ1_RPV_7", "PEQ2_RPV_7", "PEQ3_RPV_7"),
  RepVic_9 = c("PEQ1_RPV_9", "PEQ2_RPV_9", "PEQ3_RPV_9"),
  RepVic_11 = c("PEQ1_RPV_11", "PEQ2_RPV_11", "PEQ3_RPV_11"),
  RepVic_13 = c("PEQ1_RPV_13", "PEQ2_RPV_13", "PEQ3_RPV_13"),
  #moderators 
  ## drop item 5 of parental monitoring for bad loading - Cite Scharf et al 2026
  mPmonY_1 = c("Y1pMon_1", "Y2pMon_1", "Y3pMon_1", "Y4pMon_1"),#, "Y5pMon_1"), 
  mPmonY_3 = c("Y1pMon_3", "Y2pMon_3", "Y3pMon_3", "Y4pMon_3"),#, "Y5pMon_3"),
  mPmonY_5 = c("Y1pMon_5", "Y2pMon_5", "Y3pMon_5", "Y4pMon_5"),# "Y5pMon_5"),
  mPmonY_7 = c("Y1pMon_7", "Y2pMon_7", "Y3pMon_7", "Y4pMon_7"),# "Y5pMon_7"),
  mPmonY_9 = c("Y1pMon_9", "Y2pMon_9", "Y3pMon_9", "Y4pMon_9"),#, "Y5pMon_9"),
  mPmonY_11 = c("Y1pMon_11", "Y2pMon_11", "Y3pMon_11", "Y4pMon_11"), # Y5pMon_11),
  mPmonY_13 = c("Y1pMon_13", "Y2pMon_13", "Y3pMon_13", "Y4pMon_13") # Y5pMon_13),
)


## RUN 
out <- make_scale_means(df_RED, item_sets)

# GET DESCRIPTIVES 
out$descriptives

df_RED_ADD <- as.data.frame(out$data)

DESCRIPS <- as.data.frame(out$descriptives)

#TODO: Get all 5 items for peer network health and check psychometrics - using ABCD Sum for now

## Cyberbullying 
# the items do nto simply scale together. Intead, we are going to use
# CBB3_FRQ (mh_y_cb_001a__01__01): How often have you been Cyberbullied in the last year? We will input 
# 0's (None) from question CBB2_REC (mh_y_cb_001a__01) asking if they had been cyberbullied in the last year. 

################################################################################
##################### COMPUTE WAVE-SPECIFIC CYBERBULLYING FREQUENCY SCORE

# CBB2_REC indicates whether cyberbullying occurred in the past 12 months.
# CBB3_FRQ indicates frequency among those with past-year cyberbullying.
# New score = 0 if CBB2_REC == 0; otherwise uses CBB3_FRQ.

df_RED_ADD <- df_RED_ADD %>%
  dplyr::mutate(
    CyBlFq_5  = dplyr::case_when(
      CBB2_REC_5 == 0 ~ 0,
      CBB2_REC_5 == 1 ~ CBB3_FRQ_5,
      TRUE ~ NA_real_
    ),
    
    CyBlFq_7  = dplyr::case_when(
      CBB2_REC_7 == 0 ~ 0,
      CBB2_REC_7 == 1 ~ CBB3_FRQ_7,
      TRUE ~ NA_real_
    ),
    
    CyBlFq_9  = dplyr::case_when(
      CBB2_REC_9 == 0 ~ 0,
      CBB2_REC_9 == 1 ~ CBB3_FRQ_9,
      TRUE ~ NA_real_
    ),
    
    CyBlFq_11 = dplyr::case_when(
      CBB2_REC_11 == 0 ~ 0,
      CBB2_REC_11 == 1 ~ CBB3_FRQ_11,
      TRUE ~ NA_real_
    ),
    
    CyBlFq_13 = dplyr::case_when(
      CBB2_REC_13 == 0 ~ 0,
      CBB2_REC_13 == 1 ~ CBB3_FRQ_13,
      TRUE ~ NA_real_
    )
  )


## REMOVE VARIABLES 
names(df_RED_ADD)

df_RED_ADD_RED <- df_RED_ADD %>%
  select(-c(
    "CBB1_PRS_5",      "CBB1_PRS_7",      "CBB1_PRS_9",      "CBB1_PRS_11",     "CBB1_PRS_13",     "CBB2_REC_5" ,     "CBB2_REC_7"  ,    "CBB2_REC_9" ,     "CBB2_REC_11",    
    "CBB2_REC_13",     "CBB3_FRQ_5",      "CBB3_FRQ_7",      "CBB3_FRQ_9"  ,    "CBB3_FRQ_11" ,    "CBB3_FRQ_13" ,    "CBB4_PWR_5"  ,    "CBB4_PWR_7" ,     "CBB4_PWR_9",     
    "CBB4_PWR_11",     "CBB4_PWR_13",     "CExtSTP_1", 
    "Fcon1Y_1",        "Fcon1Y_3",        "Fcon1Y_5",        "Fcon1Y_7"   ,     "Fcon1Y_9"  ,      "Fcon1Y_11"  ,     "Fcon1Y_13" ,     
    "Fcon2rY_1",       "Fcon2rY_3",       "Fcon2rY_5",       "Fcon2rY_7"  ,     "Fcon2rY_9" ,      "Fcon2rY_11" ,     "Fcon2rY_13" ,     "Fcon3Y_1"   ,     "Fcon3Y_3",       
    "Fcon3Y_5",        "Fcon3Y_7",        "Fcon3Y_9" ,       "Fcon3Y_11"  ,     "Fcon3Y_13" ,      "Fcon4rY_1"  ,     "Fcon4rY_3" ,      "Fcon4rY_5"  ,     "Fcon4rY_7",      
     "Fcon4rY_9",       "Fcon4rY_11",      "Fcon4rY_13",      "Fcon5Y_1"  ,      "Fcon5Y_3" ,       "Fcon5Y_5"  ,      "Fcon5Y_7" ,       "Fcon5Y_9"  ,      "Fcon5Y_11",      
     "Fcon5Y_13",       "Fcon6Y_1",        "Fcon6Y_3",        "Fcon6Y_5"  ,      "Fcon6Y_7" ,       "Fcon6Y_9"  ,      "Fcon6Y_11",       "Fcon6Y_13" ,      "Fcon7rY_1",      
     "Fcon7rY_3",       "Fcon7rY_5",       "Fcon7rY_7",       "Fcon7rY_9" ,      "Fcon7rY_11",      "Fcon7rY_13",      "Fcon8Y_1" ,       "Fcon8Y_3"   ,     "Fcon8Y_5" ,      
    "Fcon8Y_7",        "Fcon8Y_9",        "Fcon8Y_11" ,      "Fcon8Y_13"  ,     "Fcon9Y_1",        "Fcon9Y_3"   ,     "Fcon9Y_5"  ,      "Fcon9Y_7"   ,     "Fcon9Y_9",       
    "Fcon9Y_11",       "Fcon9Y_13",       "FconMy_1"  ,      "FconMy_3"   ,     "FconMy_5"  ,      "FconMy_7"   ,     "FconMy_9"  ,      "FconMy_11"  ,     "FconMy_13",      
    "PEQ1_OV_5",       "PEQ1_OV_7",       "PEQ1_OV_9" ,      "PEQ1_OV_11" ,     "PEQ1_OV_13"  ,    "PEQ1_RLV_5" ,     "PEQ1_RLV_7" ,     "PEQ1_RLV_9" ,     "PEQ1_RLV_11",    
    "PEQ1_RLV_13",     "PEQ1_RPV_5" ,     "PEQ1_RPV_7",      "PEQ1_RPV_9" ,     "PEQ1_RPV_11" ,    "PEQ1_RPV_13",     "PEQ2_OV_5" ,      "PEQ2_OV_7"   ,    "PEQ2_OV_9" ,     
    "PEQ2_OV_11",      "PEQ2_OV_13" ,     "PEQ2_RLV_5",      "PEQ2_RLV_7" ,     "PEQ2_RLV_9",      "PEQ2_RLV_11" ,    "PEQ2_RLV_13",     "PEQ2_RPV_5" ,     "PEQ2_RPV_7",     
    "PEQ2_RPV_9",      "PEQ2_RPV_11",     "PEQ2_RPV_13",     "PEQ3_OV_5"  ,     "PEQ3_OV_7" ,      "PEQ3_OV_9"  ,     "PEQ3_OV_11",      "PEQ3_OV_13" ,     "PEQ3_RLV_5",     
    "PEQ3_RLV_7",      "PEQ3_RLV_9",      "PEQ3_RLV_11",     "PEQ3_RLV_13",     "PEQ3_RPV_5" ,     "PEQ3_RPV_7"  ,    "PEQ3_RPV_9",      "PEQ3_RPV_11" ,    "PEQ3_RPV_13",    
    "PEQs_OV_5",       "PEQs_OV_7",       "PEQs_OV_9" ,      "PEQs_OV_11",      "PEQs_OV_13"  ,    "PEQs_RLV_5" ,     "PEQs_RLV_7" ,     "PEQs_RLV_9"  ,    "PEQs_RLV_11",    
    "PEQs_RLV_13",     "PEQs_RPV_5",      "PEQs_RPV_7",      "PEQs_RPV_9",      "PEQs_RPV_11"  ,   "PEQs_RPV_13" ,    "PNH1y_5"   ,      "PNH1y_7"     ,    "PNH1y_9"   ,     
    "PNH1y_13",        "PNH2y_5",         "PNH2y_7"   ,      "PNH2y_9"  ,       "PNH2y_13"  ,      "PNH3_1y_5"  ,     "PNH3_1y_7" ,      "PNH3_1y_9"  ,     "PNH3_1y_13",     
    "PNH3y_5",         "PNH3y_7",         "PNH3y_9"   ,      "PNH3y_13",
    "Y1pMon_1",        "Y1pMon_3" ,       "Y1pMon_5"  ,      "Y1pMon_7" ,      
    "Y1pMon_9",        "Y1pMon_11",       "Y1pMon_13" ,      "Y2pMon_1" ,       "Y2pMon_3"   ,     "Y2pMon_5" ,       "Y2pMon_7"  ,      "Y2pMon_9"   ,     "Y2pMon_11" ,     
    "Y2pMon_13",       "Y3pMon_1",        "Y3pMon_3"  ,      "Y3pMon_5" ,       "Y3pMon_7"   ,     "Y3pMon_9" ,       "Y3pMon_11"  ,     "Y3pMon_13"  ,     "Y4pMon_1"  ,     
    "Y4pMon_3",        "Y4pMon_5",        "Y4pMon_7"  ,      "Y4pMon_9"  ,      "Y4pMon_11"  ,     "Y4pMon_13",       "Y5pMon_1" ,       "Y5pMon_3"   ,     "Y5pMon_5"  ,     
    "Y5pMon_7",        "Y5pMon_9",        "Y5pMon_11" ,      "Y5pMon_13",       "YMpMon_1"   ,     "YMpMon_3"  ,      "YMpMon_5"  ,      "YMpMon_7"   ,     "YMpMon_9"  ,     
    "YMpMon_11",       "YMpMon_13"
    ))



## centering function
CenterForge <- function(data, vars,
                        types = c("GMC", "LMC", "PMC"),
                        save_means = TRUE,
                        set_name = NULL,
                        check_numeric = TRUE,
                        poly = c("none", "quadratic", "cubic")) {
  
  if (!is.data.frame(data)) {
    stop("data must be a data frame.")
  }
  
  if (missing(vars) || length(vars) == 0) {
    stop("You must provide at least one variable name in 'vars'.")
  }
  
  if (!all(vars %in% names(data))) {
    missing_vars <- vars[!vars %in% names(data)]
    stop("These variables are not in the data: ",
         paste(missing_vars, collapse = ", "))
  }
  
  types <- unique(toupper(types))
  allowed_types <- c("GMC", "LMC", "PMC")
  
  if (!all(types %in% allowed_types)) {
    bad_types <- types[!types %in% allowed_types]
    stop("Invalid type(s): ", paste(bad_types, collapse = ", "),
         ". Allowed types are: ", paste(allowed_types, collapse = ", "))
  }
  
  poly <- match.arg(poly)
  
  df_sub <- data[, vars, drop = FALSE]
  
  if (check_numeric) {
    non_numeric <- vars[!vapply(df_sub, is.numeric, logical(1))]
    if (length(non_numeric) > 0) {
      stop("These variables are not numeric: ",
           paste(non_numeric, collapse = ", "))
    }
  }
  
  out <- data
  
  if (is.null(set_name)) {
    set_name <- "SET"
  }
  
  # helper to add polynomial terms
  add_poly_terms <- function(df, centered_names, poly_type) {
    if (poly_type %in% c("quadratic", "cubic")) {
      for (nm in centered_names) {
        df[[paste0("Q", nm)]] <- df[[nm]]^2
      }
    }
    
    if (poly_type == "cubic") {
      for (nm in centered_names) {
        df[[paste0("C", nm)]] <- df[[nm]]^3
      }
    }
    
    df
  }
  
  # GMC
  if ("GMC" %in% types) {
    var_means <- vapply(df_sub, mean, numeric(1), na.rm = TRUE)
    
    if (save_means) {
      for (j in seq_along(var_means)) {
        out[[paste0("MEAN_GMC_", vars[j])]] <- var_means[j]
      }
    }
    
    gmc_df <- as.data.frame(
      Map(function(x, m) x - m, df_sub, var_means)
    )
    gmc_names <- paste0("GMC_", vars)
    names(gmc_df) <- gmc_names
    out <- cbind(out, gmc_df)
    
    out <- add_poly_terms(out, gmc_names, poly)
  }
  
  # LMC
  if ("LMC" %in% types) {
    pooled_mean <- mean(as.matrix(df_sub), na.rm = TRUE)
    
    if (save_means) {
      out[[paste0("MEAN_LMC_", set_name)]] <- pooled_mean
    }
    
    lmc_df <- as.data.frame(
      lapply(df_sub, function(x) x - pooled_mean)
    )
    lmc_names <- paste0("LMC_", vars)
    names(lmc_df) <- lmc_names
    out <- cbind(out, lmc_df)
    
    out <- add_poly_terms(out, lmc_names, poly)
  }
  
  # PMC
  if ("PMC" %in% types) {
    person_means <- rowMeans(df_sub, na.rm = TRUE)
    person_means[is.nan(person_means)] <- NA
    
    if (save_means) {
      out[[paste0("MEAN_PMC_", set_name)]] <- person_means
    }
    
    pmc_df <- as.data.frame(
      lapply(df_sub, function(x) x - person_means)
    )
    pmc_names <- paste0("PMC_", vars)
    names(pmc_df) <- pmc_names
    out <- cbind(out, pmc_df)
    
    out <- add_poly_terms(out, pmc_names, poly)
  }
  
  return(out)
}


#### APPLIED VARIABLE CENTERING 
names(df_RED_ADD_RED)
#covariates
age = c("Y_AGE_1", "Y_AGE_3", "Y_AGE_5", "Y_AGE_7", "Y_AGE_9",  "Y_AGE_11", "Y_AGE_13" ) # pooled longitudinal grand mean centered
demo = c("HiParEdu", "INCOME6L", "Y_SEX") # grand mean centered 
Scrn = c("WKD_SMU_1", "WKD_SMU_3", "WKD_SMU_9") #pooled longitudinal grand mean centered
int = c("CIntSTP_1", "CIntSTP_3", "CIntSTP_5", "CIntSTP_7", "CIntSTP_9", "CIntSTP_11", "CIntSTP_13")

df_RED_ADD_RED_CENT <- CenterForge(
  data = df_RED_ADD_RED,
  vars = age,
  types = c("LMC"),
  save_means = FALSE,
  set_name = "age",
  poly = "none"
)

df_RED_ADD_RED_CENT <- CenterForge(
  data = df_RED_ADD_RED_CENT,
  vars = demo,
  types = c("GMC"),
  save_means = FALSE,
  set_name = "demo",
  poly = "none"
)

df_RED_ADD_RED_CENT <- CenterForge(
  data = df_RED_ADD_RED_CENT,
  vars = Scrn,
  types = c("LMC"),
  save_means = FALSE,
  set_name = "Scrn",
  poly = "none"
)

df_RED_ADD_RED_CENT <- CenterForge(
  data = df_RED_ADD_RED_CENT,
  vars = int,
  types = c("LMC"),
  save_means = FALSE,
  set_name = "int",
  poly = "none"
)

#predictors
Fcon = c("mFconY_1", "mFconY_3", "mFconY_5", "mFconY_7", "mFconY_9", "mFconY_11", "mFconY_13") #PMC with GMC
OvVc = c("OvtVic_5", "OvtVic_7", "OvtVic_9", "OvtVic_11", "OvtVic_13") #PMC with GMC
RLVC = c("RelVic_5",  "RelVic_7", "RelVic_9", "RelVic_11", "RelVic_13") #PMC with GMC
RPVC = c("RepVic_5", "RepVic_7", "RepVic_9", "RepVic_11", "RepVic_13")  #PMC with GMC
CYBL = c("CyBlFq_5", "CyBlFq_7", "CyBlFq_9", "CyBlFq_11", "CyBlFq_13") #PMC with GMC

df_RED_ADD_RED_CENT <- CenterForge(
  data = df_RED_ADD_RED_CENT,
  vars = Fcon,
  types = c("PMC"),
  save_means = TRUE,
  set_name = "Fcon",
  poly = "none"
)

df_RED_ADD_RED_CENT <- CenterForge(
  data = df_RED_ADD_RED_CENT,
  vars = OvVc,
  types = c("PMC"),
  save_means = TRUE,
  set_name = "OvVc",
  poly = "none"
)

df_RED_ADD_RED_CENT <- CenterForge(
  data = df_RED_ADD_RED_CENT,
  vars = RLVC,
  types = c("PMC"),
  save_means = TRUE,
  set_name = "RLVC",
  poly = "none"
)

df_RED_ADD_RED_CENT <- CenterForge(
  data = df_RED_ADD_RED_CENT,
  vars = RPVC,
  types = c("PMC"),
  save_means = TRUE,
  set_name = "RPVC",
  poly = "none"
)

df_RED_ADD_RED_CENT <- CenterForge(
  data = df_RED_ADD_RED_CENT,
  vars = CYBL,
  types = c("PMC"),
  save_means = TRUE,
  set_name = "CYBL",
  poly = "none"
)

#moderators
Pmon = c("mPmonY_1", "mPmonY_3", "mPmonY_5", "mPmonY_7", "mPmonY_9", "mPmonY_11", "mPmonY_13") #PMC with GMC
PNH = c("PNHyM_5", "PNHyM_7", "PNHyM_9", "PNHyM_13") #PMC with GMC

df_RED_ADD_RED_CENT <- CenterForge(
  data = df_RED_ADD_RED_CENT,
  vars = Pmon,
  types = c("PMC"),
  save_means = TRUE,
  set_name = "Pmon",
  poly = "none"
)

df_RED_ADD_RED_CENT <- CenterForge(
  data = df_RED_ADD_RED_CENT,
  vars = PNH,
  types = c("PMC"),
  save_means = TRUE,
  set_name = "PNH",
  poly = "none"
)


################################################################################
#################### reduce data 

df_RED_ADD_RED_CENT_RED <- df_RED_ADD_RED_CENT %>%
  select(-c("mPmonY_1", "mPmonY_3", "mPmonY_5", "mPmonY_7", "mPmonY_9", "mPmonY_11", "mPmonY_13",
            "PNHyM_5", "PNHyM_7", "PNHyM_9", "PNHyM_13",
            "mFconY_1", "mFconY_3", "mFconY_5", "mFconY_7", "mFconY_9", "mFconY_11", "mFconY_13", 
            "OvtVic_5", "OvtVic_7", "OvtVic_9", "OvtVic_11", "OvtVic_13", 
            "RelVic_5",  "RelVic_7", "RelVic_9", "RelVic_11", "RelVic_13",
            "RepVic_5", "RepVic_7", "RepVic_9", "RepVic_11", "RepVic_13",
            "CyBlFq_5", "CyBlFq_7", "CyBlFq_9", "CyBlFq_11", "CyBlFq_13",
            "HiParEdu", "WKD_SMU_1", "WKD_SMU_3", "WKD_SMU_9"))

################################################################################
#################### SAVE DATA  
names(df_RED_ADD_RED_CENT_RED)

# save dataframe
write.csv(
  df_RED_ADD_RED_CENT_RED, 
  file.path(OUT_DIR, paste0(FILENAME, "_BASE_DATA_Vprep_", DATE_STAMP, ".csv")), 
  row.names=FALSE, na="")

# Save descriptive statstics
write.csv(
  DESCRIPS,
  file.path(OUT_DIR, paste0(FILENAME, "_Var_Descriptives_", DATE_STAMP, ".csv")),
  row.names = FALSE,
  na = ""
)

## Mplus .dat and input file 
prepareMplusData(df_RED_ADD_RED_CENT_RED,file.path(OUT_DIR,paste0(FILENAME, "_BASE_DATA_Vprep_", DATE_STAMP, ".dat")), inpfile =T)


################################################################################
#################### WAVE MISSINGNESS 

library(dplyr)

waves <- c("_1", "_3", "_5", "_7", "_9", "_11", "_13")

wave_qc <- lapply(waves, function(w) {
  
  # Identify variables ending in the wave suffix
  vars_w <- grep(paste0(w, "$"), names(df_RED_ADD_RED_CENT_RED), value = TRUE)
  
  cat("\n=============================\n")
  cat("Wave suffix:", w, "\n")
  cat("Variables identified:\n")
  print(vars_w)
  cat("Number of variables:", length(vars_w), "\n")
  
  if (length(vars_w) == 0) {
    return(data.frame(
      wave = w,
      n_vars = 0,
      n_all_na = NA,
      n_non_na = NA,
      total_n = nrow(df_RED_ADD_RED_CENT_RED)
    ))
  }
  
  # Rows that are ALL NA for this wave
  all_na <- df_RED_ADD_RED_CENT_RED %>%
    select(all_of(vars_w)) %>%
    apply(1, function(x) all(is.na(x)))
  
  data.frame(
    wave = w,
    n_vars = length(vars_w),
    n_all_na = sum(all_na),
    n_non_na = sum(!all_na),
    total_n = nrow(df_RED_ADD_RED_CENT_RED)
  )
})

wave_qc <- bind_rows(wave_qc)
wave_qc





