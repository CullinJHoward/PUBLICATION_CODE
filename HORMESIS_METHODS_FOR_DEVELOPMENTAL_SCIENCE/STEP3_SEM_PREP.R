#Library

library(dplyr)
library(psych)
library(stringr)
library(purrr)
library(MplusAutomation)
library(tidyverse)

# Set working directory 

setwd("/home/cjh37695/ABCD_PROJECTS/HORMESIS_METHODS/")

# LOAD DF 

df <- read.csv("ABCD_HORM_METH_PREP_4.23.26.csv")

options(max.print = 10000)
names(df)

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


## USE IT 

item_sets <- list(
  mLkCogEn = c("ErlHs1PR","ErlHs2PR","ErlHs3PR","ErlHs4PR","ErlHs5PR",
               "ErlHs6PR","ErlHs7PR","ErlHs8PR","ErlHs9PR","ErlHs10P_11"),
  # mFconP_1 = c("Fcon1P_1","Fcon2rP_1","Fcon3P_1","Fcon4rP_1","Fcon5P_1","Fcon6P_1",
  #              "Fcon8P_1" ),
  # mFconP_3 = c("Fcon1P_3","Fcon2rP_3","Fcon3P_3","Fcon4rP_3","Fcon5P_3","Fcon6P_3",
  #              "Fcon8P_3" ),
  # mFconP_5 = c("Fcon1P_5","Fcon2rP_5","Fcon3P_5","Fcon4rP_5","Fcon5P_5","Fcon6P_5",
  #              "Fcon8P_5"),
  # mFconP_7 = c("Fcon1P_7","Fcon2rP_7","Fcon3P_7","Fcon4rP_7","Fcon5P_7","Fcon6P_7",
  #              "Fcon8P_7"),
  # mFconP_9 = c("Fcon1P_9","Fcon2rP_9","Fcon3P_9","Fcon4rP_9","Fcon5P_9","Fcon6P_9",
  #              "Fcon8P_9"),
  # mFconP_11 = c("Fcon1P_11","Fcon2rP_11","Fcon3P_11","Fcon4rP_11","Fcon5P_11","Fcon6P_11",
  #              "Fcon8P_11"),
  # mFconY_1 = c("Fcon1Y_1","Fcon2rY_1","Fcon3Y_1","Fcon4rY_1","Fcon5Y_1","Fcon6Y_1",
  #              "Fcon8Y_1" ),
  # mFconY_3 = c("Fcon1Y_3","Fcon2rY_3","Fcon3Y_3","Fcon4rY_3","Fcon5Y_3","Fcon6Y_3",
  #              "Fcon8Y_3" ),
  # mFconY_5 = c("Fcon1Y_5","Fcon2rY_5","Fcon3Y_5","Fcon4rY_5","Fcon5Y_5","Fcon6Y_5",
  #              "Fcon8Y_5" ),
  # mFconY_7 = c("Fcon1Y_7","Fcon2rY_7","Fcon3Y_7","Fcon4rY_7","Fcon5Y_7","Fcon6Y_7",
  #              "Fcon8Y_7" ),
  # mFconY_9 = c("Fcon1Y_9","Fcon2rY_9","Fcon3Y_9","Fcon4rY_9","Fcon5Y_9","Fcon6Y_9",
  #              "Fcon8Y_9" ),
  # mFconY_11 = c("Fcon1Y_11","Fcon2rY_11","Fcon3Y_11","Fcon4rY_11","Fcon5Y_11","Fcon6Y_11",
  #               "Fcon8Y_11" ),
  mNBdang_1 = c("NBHsaf1P_1R", "NBHsaf2P_1R", "NBHsaf3P_1R"),
  mNBdang_3 = c("NBHsaf1P_3R", "NBHsaf2P_3R", "NBHsaf3P_3R"),
  mNBdang_5 = c("NBHsaf1P_5R", "NBHsaf2P_5R", "NBHsaf3P_5R"),
  mNBdang_7 = c("NBHsaf1P_7R", "NBHsaf2P_7R", "NBHsaf3P_7R"),
  mNBdang_9 = c("NBHsaf1P_9R", "NBHsaf2P_9R", "NBHsaf3P_9R"),
  mNBdang_11 = c("NBHsaf1P_11R", "NBHsaf2P_11R", "NBHsaf3P_11R"),
  mPedsup_7 = c("PNGedu1Y_7", "PNGedu2Y_7", "PNGedu3Y_7"),
  mPsupr_7 = c("PNGsup1Y_7", "PNGsup2Y_7", "PNGsup3Y_7", "PNGsup4Y_7", "PNGsup5Y_7"),
  mSchAdv_1 = c("SCHenv1Y_1R", "SCHenv2Y_1R", "SCHenv3Y_1R", "SCHenv4Y_1R", "SCHenv5Y_1R", "SCHenv6Y_1R"),
  mSchAdv_3 = c("SCHenv1Y_3R", "SCHenv2Y_3R", "SCHenv3Y_3R", "SCHenv4Y_3R", "SCHenv5Y_3R", "SCHenv6Y_3R"),
  mSchAdv_5 = c("SCHenv1Y_5R", "SCHenv2Y_5R", "SCHenv3Y_5R", "SCHenv4Y_5R", "SCHenv5Y_5R", "SCHenv6Y_5R"),
  mSchAdv_7 = c("SCHenv1Y_7R", "SCHenv2Y_7R", "SCHenv3Y_7R", "SCHenv4Y_7R", "SCHenv5Y_7R", "SCHenv6Y_7R"),
  mSchAdv_9 = c("SCHenv1Y_9R", "SCHenv2Y_9R", "SCHenv3Y_9R", "SCHenv4Y_9R", "SCHenv5Y_9R", "SCHenv6Y_9R"),
  mNgUrg_1 = c("UPPSnurg1Y_1", "UPPSnurg2Y_1", "UPPSnurg3Y_1","UPPSnurg4Y_1"),
  mNgUrg_5 = c("UPPSnurg1Y_5", "UPPSnurg2Y_5", "UPPSnurg3Y_5","UPPSnurg4Y_5"),
  mNgUrg_9 = c("UPPSnurg1Y_9", "UPPSnurg2Y_9", "UPPSnurg3Y_9","UPPSnurg4Y_9"),
  mPers_1 = c("UPPSpers1Y_1", "UPPSpers2Y_1", "UPPSpers3Y_1","UPPSpers4Y_1"),
  mPers_5 = c("UPPSpers1Y_5", "UPPSpers2Y_5", "UPPSpers3Y_5","UPPSpers4Y_5"),
  mPers_9 = c("UPPSpers1Y_9", "UPPSpers2Y_9", "UPPSpers3Y_9","UPPSpers4Y_9"),
  mPlan_1 = c("UPPSplan1Y_1", "UPPSplan2Y_1", "UPPSplan3Y_1","UPPSplan4Y_1"),
  mPlan_5 = c("UPPSplan1Y_5", "UPPSplan2Y_5", "UPPSplan3Y_5","UPPSplan4Y_5"),
  mPlan_9 = c("UPPSplan1Y_9", "UPPSplan2Y_9", "UPPSplan3Y_9","UPPSplan4Y_9"),
  mPsUrg_1 = c("UPPSpurg1Y_1", "UPPSpurg2Y_1", "UPPSpurg3Y_1","UPPSpurg4Y_1"),
  mPsUrg_5 = c("UPPSpurg1Y_5", "UPPSpurg2Y_5", "UPPSpurg3Y_5","UPPSpurg4Y_5"),
  mPsUrg_9 = c("UPPSpurg1Y_9", "UPPSpurg2Y_9", "UPPSpurg3Y_9","UPPSpurg4Y_9"),
  mSenSe_1 = c("UPPSsens1Y_1", "UPPSsens2Y_1", "UPPSsens3Y_1","UPPSsens4Y_1"),
  mSenSe_5 = c("UPPSsens1Y_5", "UPPSsens2Y_5", "UPPSsens3Y_5","UPPSsens4Y_5"),
  mSenSe_9 = c("UPPSsens1Y_9", "UPPSsens2Y_9", "UPPSsens3Y_9","UPPSsens4Y_9"),
  mBISy_1 = c("yBBbis1Y_1", "yBBbis2Y_1", "yBBbis3Y_1","yBBbis4Y_1", "yBBbis5Y_1",
              "yBBbis6Y_1", "yBBbis7Y_1"),
  mBISy_5 = c("yBBbis1Y_5", "yBBbis2Y_5", "yBBbis3Y_5","yBBbis4Y_5", "yBBbis5Y_5",
              "yBBbis6Y_5", "yBBbis7Y_5"),
  mBISy_9 = c("yBBbis1Y_9", "yBBbis2Y_9", "yBBbis3Y_9","yBBbis4Y_9",
              "yBBbis6Y_9", "yBBbis7Y_9"),
  mDRVy_1 = c("yBBdr1Y_1", "yBBdr2Y_1", "yBBdr3Y_1","yBBdr4Y_1"),
  mDRVy_5 = c("yBBdr1Y_5", "yBBdr2Y_5", "yBBdr3Y_5","yBBdr4Y_5"),
  mDRVy_9 = c("yBBdr1Y_9", "yBBdr2Y_9", "yBBdr3Y_9","yBBdr4Y_9"),
  mFnSky_1 = c("yBBfs1Y_1", "yBBfs2Y_1", "yBBfs3Y_1","yBBfs4Y_1"),
  mFnSky_5 = c("yBBfs1Y_5", "yBBfs2Y_5", "yBBfs3Y_5","yBBfs4Y_5"),
  mFnSky_9 = c("yBBfs1Y_9", "yBBfs2Y_9", "yBBfs3Y_9","yBBfs4Y_9"),
  mRwRspy_1 = c("yBBrr1Y_1", "yBBrr2Y_1", "yBBrr3Y_1","yBBrr4Y_1", "yBBrr5Y_1"),
  mRwRspy_5 = c("yBBrr1Y_5", "yBBrr2Y_5", "yBBrr3Y_5","yBBrr4Y_5", "yBBrr5Y_5"),
  mRwRspy_9 = c("yBBrr1Y_9", "yBBrr2Y_9", "yBBrr3Y_9","yBBrr4Y_9", "yBBrr5Y_9"),
  mAttune_7 = c("yDERatt1P_7", "yDERatt2P_7", "yDERatt3P_7","yDERatt4P_7", "yDERatt5P_7","yDERatt6P_7"),
  mAttune_9 = c("yDERatt1P_9", "yDERatt2P_9", "yDERatt3P_9","yDERatt4P_9", "yDERatt5P_9","yDERatt6P_9"),
  mAttune_11 = c("yDERatt1P_11", "yDERatt2P_11", "yDERatt3P_11","yDERatt4P_11", "yDERatt5P_11","yDERatt6P_11"),
  mDist_7 = c("yDERgol1P_7", "yDERgol2P_7", "yDERgol3P_7","yDERgol4P_7"),
  mDist_9 = c("yDERgol1P_9", "yDERgol2P_9", "yDERgol3P_9","yDERgol4P_9"),
  mDist_11 = c("yDERgol1P_11", "yDERgol2P_11", "yDERgol3P_11","yDERgol4P_11"),
  mNgSec_7 = c("yDERnac1P_7", "yDERnac2P_7", "yDERnac3P_7","yDERnac4P_7",
              "yDERnac5P_7", "yDERnac6P_7", "yDERnac7P_7"),
  mNgSec_9 = c("yDERnac1P_9", "yDERnac2P_9", "yDERnac3P_9","yDERnac4P_9",
               "yDERnac5P_9", "yDERnac6P_9", "yDERnac7P_9"),
  mNgSec_11 = c("yDERnac1P_11", "yDERnac2P_11", "yDERnac3P_11","yDERnac4P_11",
               "yDERnac5P_11", "yDERnac6P_11", "yDERnac7P_11"),
  mCata_7 = c("yDERstr1P_7", "yDERstr2P_7", "yDERstr3P_7","yDERstr4P_7",
               "yDERstr5P_7", "yDERstr6P_7", "yDERstr7P_7", "yDERstr8P_7",
              "yDERstr9P_7", "yDERstr10P_7", "yDERstr11P_7", "yDERstr12P_7"),
  mCata_9 = c("yDERstr1P_9", "yDERstr2P_9", "yDERstr3P_9","yDERstr4P_9",
              "yDERstr5P_9", "yDERstr6P_9", "yDERstr7P_9", "yDERstr8P_9",
              "yDERstr9P_9", "yDERstr10P_9", "yDERstr11P_9", "yDERstr12P_9"),
  mCata_11 = c("yDERstr1P_11", "yDERstr2P_11", "yDERstr3P_11","yDERstr4P_11",
              "yDERstr5P_11", "yDERstr6P_11", "yDERstr7P_11", "yDERstr8P_11",
              "yDERstr9P_11", "yDERstr10P_11", "yDERstr11P_11", "yDERstr12P_11"),
  mEmoSup_7 = c("yEmSup1Y_7", "yEmSup2Y_7", "yEmSup3Y_7"),
  mEmoSup_9 = c("yEmSup1Y_9", "yEmSup2Y_9", "yEmSup3Y_9"),
  mEmoSup_11 = c("yEmSup1Y_11", "yEmSup2Y_11", "yEmSup3Y_11"),
  mReApp_7 = c("yReAp1Y_7", "yReAp2Y_7", "yReAp3Y_7"),
  mReApp_9 = c("yReAp1Y_9", "yReAp2Y_9", "yReAp3Y_9"),
  mReApp_11 = c("yReAp1Y_11", "yReAp2Y_11", "yReAp3Y_11"),
  ## drop item 5 of parental monitoring for bad loading - Cite Scharf et al 2026
  mPmon_1 = c("P1pMon_1", "P2pMon_1", "P3pMon_1", "P4pMon_1"),#, "P5pMon_1"), 
  mPmon_3 = c("P1pMon_3", "P2pMon_3", "P3pMon_3", "P4pMon_3"),#, "P5pMon_3"),
  mPmon_5 = c("P1pMon_5", "P2pMon_5", "P3pMon_5", "P4pMon_5"),# "P5pMon_5"),
  mPmon_7 = c("P1pMon_7", "P2pMon_7", "P3pMon_7", "P4pMon_7"),# "P5pMon_7"),
  mPmon_9 = c("P1pMon_9", "P2pMon_9", "P3pMon_9", "P4pMon_9"),#, "P5pMon_9"),
  mPmon_11 = c("P1pMon_11", "P2pMon_11", "P3pMon_11", "P4pMon_11", "P5pMon_11"),
  mPNH_5 = c("PNH1y_5", "PNH2_1y_5", "PNH2y_5", "PNH3_1y_5", "PNH3y_5"),
  mPNH_7 = c("PNH1y_7", "PNH2_1y_7", "PNH2y_7", "PNH3_1y_7", "PNH3y_7"),
  mPNH_9 = c("PNH1y_9", "PNH2_1y_9", "PNH2y_9", "PNH3_1y_9", "PNH3y_9")
)


## CEHCK FOR MISSINGNESS 
missing_vars_by_scale <- lapply(names(item_sets), function(scale_name) {
  vars <- item_sets[[scale_name]]
  missing_vars <- setdiff(vars, names(df))
  if (length(missing_vars) > 0) {
    data.frame(
      scale = scale_name,
      missing_var = missing_vars
    )
  }
})

missing_vars_by_scale <- do.call(rbind, missing_vars_by_scale)
missing_vars_by_scale

## RUN 
out <- make_scale_means(df, item_sets)

# GET DESCRIPTIVES 
out$descriptives

df_NEWVARS <- as.data.frame(out$data)

DESCRIPS <- as.data.frame(out$descriptives)

## REMOVE VARIABLES 
df_NEWVARS_RED <- df_NEWVARS %>%
  select(-c(
"PBBbis5P_1",    "PBBbis5P_5",    "PBBbis5P_9",
"yAff1P_5",      "yAff2P_5",      "yAff3P_5",      "yAff4P_5",     
"yAff5P_5",      "yAff6P_5",      "yAggr1P_5",     "yAggr2P_5",     "yAggr3P_5",     "yAggr4P_5",     "yAggr5P_5",     "yAggr6P_5",    
"yAggr7P_5",     "yAttn1P_5",     "yAttn2P_5",     "yAttn3P_5",     "yAttn4P_5",     "yAttn5P_5",     "yAttn6P_5",
"yDepm1P_5",     "yDepm2P_5",     "yDepm3P_5",     "yDepm4P_5",     "yDepm5P_5",     "yEfCn1P_5",     "yEfCn2P_5",    
"yEfCn3P_5",     "yEfCn4P_5",     "yEfCn5P_5",     "yEfCn6P_5",     "yEfCn7P_5",
"yFear1P_5",     "yFear2P_5",    
"yFear3P_5",     "yFear4P_5",     "yFear5P_5",     "yFear6P_5",     "yFrus1P_5",     "yFrus2P_5",     "yFrus3P_5",     "yFrus4P_5",    
"yFrus5P_5",     "yFrus6P_5",     "yPinh1P_5",     "yPinh2P_5",     "yPinh3P_5",     "yPinh4P_5",     "yPinh5P_5",
"yShy1P_5",      "yShy2P_5",      "yShy3P_5",      "yShy4P_5",      "yShy5P_5",      "yUrg1P_5",      "yUrg2P_5",      "yUrg3P_5",     
"yUrg4P_5",      "yUrg5P_5",      "yUrg6P_5",      "yUrg7P_5",      "yUrg8P_5",      "yUrg9P_5",
"ErlHs1PR",      "ErlHs2PR", "ErlHs3PR",      "ErlHs4PR",      "ErlHs5PR",    "ErlHs6PR", "ErlHs7PR", "ErlHs8PR", "ErlHs9PR", "ErlHs10P_11",
"NBHsaf1P_1",   "NBHsaf2P_1",   "NBHsaf3P_1",   "NBHsaf1P_3",  
"NBHsaf2P_3",   "NBHsaf3P_3",   "NBHsaf1P_5",   "NBHsaf2P_5",   "NBHsaf3P_5",   "NBHsaf1P_7",   "NBHsaf2P_7",   "NBHsaf3P_7",  
"NBHsaf1P_9",   "NBHsaf2P_9",   "NBHsaf3P_9",   "NBHsaf1P_11",  "NBHsaf2P_11",  "NBHsaf3P_11",  "SCHenv1Y_1",   "SCHenv2Y_1",  
"SCHenv3Y_1",   "SCHenv4Y_1",   "SCHenv5Y_1",   "SCHenv6Y_1",   "SCHenv1Y_3",   "SCHenv2Y_3",   "SCHenv3Y_3",   "SCHenv4Y_3",  
"SCHenv5Y_3",   "SCHenv6Y_3",   "SCHenv1Y_5",   "SCHenv2Y_5",   "SCHenv3Y_5",   "SCHenv4Y_5",   "SCHenv5Y_5",   "SCHenv6Y_5",  
"SCHenv1Y_7",   "SCHenv2Y_7",   "SCHenv3Y_7",   "SCHenv4Y_7",   "SCHenv5Y_7",   "SCHenv6Y_7",   "SCHenv1Y_9",   "SCHenv2Y_9",  
"SCHenv3Y_9",   "SCHenv4Y_9",   "SCHenv5Y_9",   "SCHenv6Y_9",         
"P1pMon_1", "P2pMon_1", "P3pMon_1", "P4pMon_1", "P5pMon_1",
"P1pMon_3", "P2pMon_3", "P3pMon_3", "P4pMon_3", "P5pMon_3",
"P1pMon_5", "P2pMon_5", "P3pMon_5", "P4pMon_5", "P5pMon_5",
"P1pMon_7", "P2pMon_7", "P3pMon_7", "P4pMon_7", "P5pMon_7",
"P1pMon_9", "P2pMon_9", "P3pMon_9", "P4pMon_9", "P5pMon_9",
"P1pMon_11", "P2pMon_11", "P3pMon_11", "P4pMon_11", "P5pMon_11",
"PNH1y_5", "PNH2_1y_5", "PNH2y_5", "PNH3_1y_5", "PNH3y_5",
"PNH1y_7", "PNH2_1y_7", "PNH2y_7", "PNH3_1y_7", "PNH3y_7",
"PNH1y_9", "PNH2_1y_9", "PNH2y_9", "PNH3_1y_9", "PNH3y_9"))


################################################################################
############### AVERAGE IMAGING VARIABLES 

library(dplyr)
library(stringr)
library(tibble)
library(purrr)

#--------------------------------------------
# 1) Collect imaging variable names
#--------------------------------------------
img_vars <- names(df_NEWVARS_RED)[str_detect(names(df_NEWVARS_RED), "^(dti_FA|rs)")]

#--------------------------------------------
# 2) Parse variable names
#    DTI pattern: dti_FA_CNCG_L_1
#    RS pattern : rsDMN_Lamy_1
#--------------------------------------------
parse_imaging_name <- function(x) {
  
  if (str_detect(x, "^dti_FA_")) {
    m <- str_match(x, "^(dti_FA)_([A-Za-z0-9]+)_([LR])_(\\d+)$")
    
    if (!is.na(m[1, 1])) {
      return(tibble(
        var      = x,
        modality = m[1, 2],
        region   = m[1, 3],
        hemi     = m[1, 4],
        wave     = m[1, 5],
        type     = "DTI"
      ))
    }
  }
  
  if (str_detect(x, "^rs")) {
    m <- str_match(x, "^(rs[A-Za-z0-9]+)_([LR])([A-Za-z0-9]+)_(\\d+)$")
    
    if (!is.na(m[1, 1])) {
      return(tibble(
        var      = x,
        modality = m[1, 2],
        region   = m[1, 4],
        hemi     = m[1, 3],
        wave     = m[1, 5],
        type     = "RS"
      ))
    }
  }
  
  # not a lateralized imaging variable
  tibble(
    var      = x,
    modality = NA_character_,
    region   = NA_character_,
    hemi     = NA_character_,
    wave     = NA_character_,
    type     = NA_character_
  )
}

parsed_vars <- bind_rows(lapply(img_vars, parse_imaging_name))

#--------------------------------------------
# 3) Keep only successfully parsed lateralized vars
#--------------------------------------------
lat_vars <- parsed_vars %>%
  filter(!is.na(modality), hemi %in% c("L", "R"), !is.na(wave))

#--------------------------------------------
# 4) Find L/R pairs within modality-region-wave
#--------------------------------------------
pair_map <- lat_vars %>%
  select(var, modality, region, hemi, wave, type) %>%
  tidyr::pivot_wider(
    names_from  = hemi,
    values_from = var
  ) %>%
  mutate(
    new_mean_var = paste0(modality, "_", region, "_M_", wave),
    note = case_when(
      !is.na(L) & !is.na(R) ~ "complete L/R pair",
      !is.na(L) &  is.na(R) ~ "only L variable present; no pair created",
      is.na(L) & !is.na(R) ~ "only R variable present; no pair created",
      TRUE                  ~ "no usable lateralized variables"
    )
  )

# Keep only complete pairs for averaging/correlation creation
complete_pairs <- pair_map %>%
  filter(!is.na(L) & !is.na(R))

#--------------------------------------------
# 5) Create person-level mean vars and QC stats
#    Rule:
#    - if both missing -> NA
#    - if one present -> use that value
#    - if both present -> mean(L, R)
#--------------------------------------------
qc_results <- vector("list", nrow(complete_pairs))

for (i in seq_len(nrow(complete_pairs))) {
  
  Lvar <- complete_pairs$L[i]
  Rvar <- complete_pairs$R[i]
  Mvar <- complete_pairs$new_mean_var[i]
  
  # create person-level mean with your missing-data rule
  df_NEWVARS_RED[[Mvar]] <- rowMeans(df_NEWVARS_RED[, c(Lvar, Rvar)], na.rm = TRUE)
  df_NEWVARS_RED[[Mvar]][is.nan(df_NEWVARS_RED[[Mvar]])] <- NA_real_
  
  # sample-level correlation across people for that wave
  pair_cor <- suppressWarnings(
    cor(df[[Lvar]], df_NEWVARS_RED[[Rvar]], use = "pairwise.complete.obs", method = "pearson")
  )
  
  qc_results[[i]] <- tibble(
    modality      = complete_pairs$modality[i],
    region        = complete_pairs$region[i],
    wave          = complete_pairs$wave[i],
    var_L         = Lvar,
    var_R         = Rvar,
    correlation   = unname(pair_cor),
    new_mean_var  = Mvar,
    mean_M        = mean(df_NEWVARS_RED[[Mvar]], na.rm = TRUE),
    min_M         = min(df_NEWVARS_RED[[Mvar]], na.rm = TRUE),
    max_M         = max(df_NEWVARS_RED[[Mvar]], na.rm = TRUE),
    note          = complete_pairs$note[i]
  )
}

qc_table_complete <- bind_rows(qc_results)

#--------------------------------------------
# 6) Notes for incomplete structures
#--------------------------------------------
qc_table_incomplete <- pair_map %>%
  filter(is.na(L) | is.na(R)) %>%
  transmute(
    modality     = modality,
    region       = region,
    wave         = wave,
    var_L        = L,
    var_R        = R,
    correlation  = NA_real_,
    new_mean_var = new_mean_var,
    mean_M       = NA_real_,
    min_M        = NA_real_,
    max_M        = NA_real_,
    note         = note
  )

#--------------------------------------------
# 7) Final QC table
#--------------------------------------------
qc_table_lr_means <- bind_rows(qc_table_complete, qc_table_incomplete) %>%
  arrange(modality, region, wave)

# View results
qc_table_lr_means

names(df_NEWVARS_RED)

################################################################################
#################### Impute motion data 

names(df_NEWVARS_RED)

motion_vars_rs <- c("rsMnMOT_1", "rsMnMOT_5", "rsMnMOT_9")
motion_vars_dti <- c("dtiMnMOT_1", "dtiMnMOT_5", "dtiMnMOT_9")

# Save wave-specific sample means before imputing
wave_means_rs <- sapply(df_NEWVARS_RED[motion_vars_rs], function(x) mean(x, na.rm = TRUE))
wave_means_dti <- sapply(df_NEWVARS_RED[motion_vars_dti], function(x) mean(x, na.rm = TRUE))

# Function to impute one person's motion values
impute_motion_row <- function(x, wave_means) {
  # x is a numeric vector of length 3: _1, _5, _9
  observed <- !is.na(x)
  
  if (sum(observed) >= 1) {
    # If at least one value exists, use the person's mean of available scores
    person_mean <- mean(x[observed], na.rm = TRUE)
    x[!observed] <- person_mean
  } else {
    # If all are missing, use sample mean for each wave
    x <- wave_means
  }
  
  return(x)
}

# Apply row-wise
df_NEWVARS_RED[motion_vars_rs] <- t(apply(df_NEWVARS_RED[motion_vars_rs], 1, impute_motion_row, wave_means = wave_means_rs))
df_NEWVARS_RED[motion_vars_dti] <- t(apply(df_NEWVARS_RED[motion_vars_dti], 1, impute_motion_row, wave_means = wave_means_dti))

################################################################################
############### CLEAN UP THE DATASET 
names(df_NEWVARS_RED)

# ADD A NUMERIC ID 

df_NEWVARS_RED$NUMID <- row.names(df_NEWVARS_RED)

#REDUCE DF 
df_NEWVARS_RED_RED <- df_NEWVARS_RED %>%
  select(c("subID", "SiteID", "FamilyID", "ppensity", "SCANman",  "SCANmod", "scanID", "NUMID",
           "INCOME6L",  "Y_HISP", "Y_SEX", "HiParEdu_1",  "MarWrkSt_1", "HiParEdu_1R",
           "PubDvFm_1",       "PubDvFm_3",        "PubDvFm_5",       "PubDvFm_7",       "PubDvFm_9",       "PubDvFm_11",
           "PubDvMl_1",        "PubDvMl_3",       "PubDvMl_5",       "PubDvMl_7",       "PubDvMl_9",       "PubDvMl_11",
           "YHghtm_1",        "YHghtm_3",        "YHghtm_5",        "YHghtm_7",        "YHghtm_9",        "YHghtm_11",      
           "YWghtm_1",        "YWghtm_3",        "YWghtm_5",        "YWghtm_7",        "YWghtm_9",        "YWghtm_11",       
           "Y_AGE_1",         "Y_AGE_3",        "Y_AGE_5",         "Y_AGE_7",         "Y_AGE_9",         "Y_AGE_11",  
           "PRBlPov_1",       "PRCrowd_1",       "PRLowEd_1",       "PRUnEmp_1",     "SocMob_1R",      "HiParEdu_1R",
            "CAggSTP_1",       "CAggSTP_3",       "CAggSTP_5",       "CAggSTP_7",       "CAggSTP_9",       "CAggSTP_11",     
           "CAttSTP_1",       "CAttSTP_3",       "CAttSTP_5",       "CAttSTP_7",       "CAttSTP_9",       "CAttSTP_11",      
           "CAxDpTP_1",       "CAxDpTP_3",      "CAxDpTP_5",       "CAxDpTP_7",      "CAxDpTP_9",       "CAxDpTP_11",      
           "CExtSTP_1",       "CExtSTP_3",       "CExtSTP_5",       "CExtSTP_7",       "CExtSTP_9",       "CExtSTP_11",      
           "CIntSTP_1",       "CIntSTP_3",       "CIntSTP_5",       "CIntSTP_7",       "CIntSTP_9",       "CIntSTP_11",     
           "CRlBrMP_1",       "CRlBrMP_3",       "CRlBrMP_5",       "CRlBrMP_7",       "CRlBrMP_9",       "CRlBrMP_11",      
           "CSocSTP_1",       "CSocSTP_3",       "CSocSTP_5",       "CSocSTP_7",       "CSocSTP_9",       "CSocSTP_11",      
           "CSomSTP_1",       "CSomSTP_3",       "CSomSTP_5",       "CSomSTP_7",       "CSomSTP_9",       "CSomSTP_11",      
           "CThtSTP_1",       "CThtSTP_3",       "CThtSTP_5",       "CThtSTP_7",       "CThtSTP_9",       "CThtSTP_11",     
           "CWtDpTP_1",       "CWtDpTP_3",       "CWtDpTP_5",       "CWtDpTP_7",       "CWtDpTP_9",       "CWtDpTP_11",      
           "DysBpM_5",        "DysBpM_7",       "DysBpM_9",        "DysBpM_11", 
           "SysBpM_5",        "SysBpM_7",        "SysBpM_9",        "SysBpM_11",
           "mLkCogEn",      
           "mNBdang_1",       "mNBdang_3",      "mNBdang_5",        "mNBdang_7",      "mNBdang_9",         "mNBdang_11",
           "mSchAdv_1",       "mSchAdv_3",       "mSchAdv_5",       "mSchAdv_7",       "mSchAdv_9",
           "FAMCONY_1",       "FAMCONY_3",       "FAMCONY_5",       "FAMCONY_7",       "FAMCONY_9",       "FAMCONY_11",
           "mNgUrg_1",        "mNgUrg_5",        "mNgUrg_9",   
           "mPsUrg_1",        "mPsUrg_5",        "mPsUrg_9",  
           "mPers_1",         "mPers_5",         "mPers_9",         
           "mPlan_1",         "mPlan_5",         "mPlan_9",              
           "mSenSe_1",        "mSenSe_5",        "mSenSe_9",        
           "mAttune_7",       "mAttune_9",       "mAttune_11",      
           "mDist_7",         "mDist_9",         "mDist_11",        
           "mEmoSup_7",       "mEmoSup_9",       "mEmoSup_11",      
           "mReApp_7",        "mReApp_9",        "mReApp_11",
           "dts_FA_ALL_1",    "dts_FA_ALL_5",    "dts_FA_ALL_9",
           "dti_FA_CNCG_M_1", "dti_FA_CNCG_M_5", "dti_FA_CNCG_M_9", 
           "dti_FA_ILF_M_1",  "dti_FA_ILF_M_5",  "dti_FA_ILF_M_9",  
           "dti_FA_SLF_M_1",  "dti_FA_SLF_M_5",  "dti_FA_SLF_M_9", 
           "dti_FA_UNC_M_1",  "dti_FA_UNC_M_5",  "dti_FA_UNC_M_9",  
           "dti_FA_iFSF_M_1", "dti_FA_iFSF_M_5", "dti_FA_iFSF_M_9", 
           "dti_FA_pSLF_M_1", "dti_FA_pSLF_M_5", "dti_FA_pSLF_M_9", 
           "dti_FA_tSLF_M_1", "dti_FA_tSLF_M_5", "dti_FA_tSLF_M_9", 
           "dtiMnMOT_1",      "dtiMnMOT_5",      "dtiMnMOT_9",
           "mPmon_1",  "mPmon_3",  "mPmon_5", "mPmon_7", "mPmon_9", "mPmon_11",       
           "mPNH_5",   "mPNH_7",   "mPNH_9"))
           
           

################################################################################
############### POMS, CENTER, AND WINSORIZE STUFF

POMS <- function(data, cols,
                 min_val = NULL, max_val = NULL,
                 scale_to = 100,
                 center = FALSE, center_at = NULL,
                 wins = FALSE, trim = 0.01,
                 plot = FALSE, bins = 30) {
  
  # Required packages:
  # dplyr, ggplot2, rlang, tidyselect, psych
  
  cols_quo  <- rlang::enquo(cols)
  col_names <- names(tidyselect::eval_select(cols_quo, data))
  
  if (length(col_names) == 0) {
    stop("No columns were selected.")
  }
  
  if (!scale_to %in% c(1, 10, 100)) {
    stop("scale_to must be one of: 1, 10, or 100.")
  }
  
  resolve_bound <- function(bound, var, i, bound_name) {
    if (is.null(bound)) return(NULL)
    
    if (length(bound) == 1) {
      return(bound)
    }
    
    if (!is.null(names(bound)) && var %in% names(bound)) {
      return(bound[[var]])
    }
    
    if (length(bound) == length(col_names)) {
      return(bound[[i]])
    }
    
    stop(bound_name, " must be NULL, a single value, a named vector, or a vector the same length as the selected columns.")
  }
  
  prefix <- dplyr::case_when(
    wins  & center  ~ "WPC_",
    wins  & !center ~ "WP_",
    !wins & center  ~ "PC_",
    TRUE            ~ "P_"
  )
  
  message("=== POMS START ===")
  message("vars=", paste(col_names, collapse = ", "),
          "; WIN=", ifelse(wins, "Y", "N"),
          if (wins) paste0(" (trim=", trim, ")") else "",
          "; SCALE=0-", scale_to,
          "; CENTER=", ifelse(center, "Y", "N"),
          if (center && is.null(center_at)) " (at mean)" else "",
          if (center && !is.null(center_at)) paste0(" (at ", center_at, ")") else "",
          "; PLOT=", ifelse(plot, "Y", "N"))
  
  out_data <- data
  
  for (i in seq_along(col_names)) {
    var    <- col_names[i]
    x_orig <- data[[var]]
    
    if (!is.numeric(x_orig)) {
      stop("Column '", var, "' is not numeric. POMS requires numeric variables.")
    }
    
    message("---- [", var, "] ----")
    
    # Step 1: winsorize first if requested
    if (wins) {
      x_work <- psych::winsor(x_orig, trim = trim)
      
      was_winsorized <- !is.na(x_orig) & !is.na(x_work) &
        !dplyr::near(x_orig, x_work)
      
      orig_min <- min(x_orig, na.rm = TRUE)
      orig_max <- max(x_orig, na.rm = TRUE)
      work_min <- min(x_work, na.rm = TRUE)
      work_max <- max(x_work, na.rm = TRUE)
      n_wins   <- sum(was_winsorized, na.rm = TRUE)
      
      message("[", var, "] WIN=Y; trim=", trim,
              "; raw min/max=", orig_min, "/", orig_max,
              "; win min/max=", work_min, "/", work_max,
              "; n winsorized=", n_wins)
    } else {
      x_work <- x_orig
      was_winsorized <- rep(FALSE, length(x_orig))
      n_wins <- 0
      
      work_min <- min(x_work, na.rm = TRUE)
      work_max <- max(x_work, na.rm = TRUE)
      
      message("[", var, "] WIN=N; raw=min/max used unless user-defined; obs min/max=",
              work_min, "/", work_max)
    }
    
    # Step 2: range decision
    user_min <- resolve_bound(min_val, var, i, "min_val")
    user_max <- resolve_bound(max_val, var, i, "max_val")
    
    if (!is.null(user_min) && !is.null(user_max)) {
      this_min <- user_min
      this_max <- user_max
      message("[", var, "] RANGE=user; min/max=", this_min, "/", this_max)
    } else if (is.null(user_min) && is.null(user_max)) {
      this_min <- min(x_work, na.rm = TRUE)
      this_max <- max(x_work, na.rm = TRUE)
      
      if (wins) {
        message("[", var, "] RANGE=observed post-win; min/max=", this_min, "/", this_max)
      } else {
        message("[", var, "] RANGE=observed raw; min/max=", this_min, "/", this_max)
      }
    } else {
      stop("[", var, "] Both min_val and max_val must be supplied together, or both left NULL.")
    }
    
    if (isTRUE(all.equal(this_max, this_min))) {
      stop("[", var, "] min and max are equal; cannot compute POMS.")
    }
    
    # Step 3: POMS transform
    x_poms <- scale_to * (x_work - this_min) / (this_max - this_min)
    message("[", var, "] POMS=Y; formula=", scale_to, "*(x-min)/(max-min)")
    
    # Step 4: centering
    if (center) {
      if (is.null(center_at)) {
        center_value <- mean(x_poms, na.rm = TRUE)
        message("[", var, "] CENTER=Y; at mean of transformed values=", round(center_value, 4))
      } else {
        center_value <- center_at
        message("[", var, "] CENTER=Y; at user-specified transformed value=", center_value)
      }
      
      x_final <- x_poms - center_value
    } else {
      x_final <- x_poms
      message("[", var, "] CENTER=N")
    }
    
    new_name <- paste0(prefix, var)
    out_data[[new_name]] <- x_final
    message("[", var, "] OUTPUT=", new_name)
    
    # Step 5: diagnostic plot
    if (plot) {
      
      df_orig <- data.frame(
        value = x_orig,
        panel = "Original",
        fill_group = factor("Not winsorized",
                            levels = c("Not winsorized", "Winsorized"))
      )
      
      df_trans <- data.frame(
        value = x_final,
        panel = "Transformed",
        fill_group = factor(
          ifelse(was_winsorized, "Winsorized", "Not winsorized"),
          levels = c("Not winsorized", "Winsorized")
        )
      )
      
      plot_df <- rbind(df_orig, df_trans)
      plot_df$panel <- factor(plot_df$panel, levels = c("Original", "Transformed"))
      
      p <- ggplot2::ggplot() +
        ggplot2::geom_histogram(
          data = subset(plot_df, panel == "Original"),
          mapping = ggplot2::aes(x = value, fill = fill_group),
          bins = bins,
          na.rm = TRUE,
          show.legend = TRUE
        ) +
        ggplot2::geom_histogram(
          data = subset(plot_df, panel == "Transformed"),
          mapping = ggplot2::aes(x = value, fill = fill_group),
          bins = bins,
          position = "stack",
          na.rm = TRUE,
          show.legend = TRUE
        ) +
        ggplot2::facet_grid(panel ~ ., scales = "free") +
        ggplot2::scale_fill_manual(
          values = c("Not winsorized" = "gray70",
                     "Winsorized" = "red"),
          drop = FALSE,
          limits = c("Not winsorized", "Winsorized")
        ) +
        ggplot2::labs(
          title = paste0("POMS Transformation (0-", scale_to, "): ", var),
          x = NULL,
          y = "Count",
          fill = NULL
        ) +
        ggplot2::theme_bw() +
        ggplot2::theme(
          legend.position = "bottom"
        )
      
      print(p)
      message("[", var, "] PLOT=printed; # winsorized=", n_wins)
    } else {
      message("[", var, "] PLOT=N")
    }
  }
  
  message("=== POMS END ===")
  return(out_data)
}


# APPLIED 
names(df_NEWVARS_RED_RED)

dev.off() # was crashing out with too many plots. can remove if it runs

df_NEWVARS_RED_RED_POMS <- df_NEWVARS_RED_RED %>%
  POMS(c(mNBdang_1, mNBdang_3, mNBdang_5, mNBdang_7, mNBdang_9, mNBdang_11), min_val = 1, max_val = 5, wins = TRUE, trim = .01,  scale_to = 10, plot = TRUE) %>%
  POMS(c(mSchAdv_1, mSchAdv_3,mSchAdv_5, mSchAdv_7, mSchAdv_9), min_val = 1, max_val = 4, wins = TRUE, trim = .01, scale_to = 10, plot = TRUE)%>%
  POMS(c(FAMCONY_1, FAMCONY_3, FAMCONY_5, FAMCONY_7, FAMCONY_9, FAMCONY_11), min_val = 0, max_val = 9, wins = TRUE, trim = .01, scale_to = 10, plot = TRUE)%>%
  POMS(c(mLkCogEn), min_val = 0, max_val = 4.8, wins = TRUE, trim = .01, scale_to = 10, plot = TRUE) %>%
  POMS(c(dtiMnMOT_1, dtiMnMOT_5, dtiMnMOT_9),  wins = TRUE, trim = .01, scale_to = 10, plot = TRUE)%>%
  POMS(c(dti_FA_CNCG_M_1, dti_FA_CNCG_M_5, dti_FA_CNCG_M_9, 
         dti_FA_ILF_M_1,  dti_FA_ILF_M_5,  dti_FA_ILF_M_9,  
         dti_FA_SLF_M_1,  dti_FA_SLF_M_5,  dti_FA_SLF_M_9, 
         dti_FA_UNC_M_1,  dti_FA_UNC_M_5,  dti_FA_UNC_M_9,  
         dti_FA_iFSF_M_1, dti_FA_iFSF_M_5, dti_FA_iFSF_M_9, 
         dti_FA_pSLF_M_1, dti_FA_pSLF_M_5, dti_FA_pSLF_M_9, 
         dti_FA_tSLF_M_1, dti_FA_tSLF_M_5, dti_FA_tSLF_M_9, 
         dts_FA_ALL_1,    dts_FA_ALL_5,    dts_FA_ALL_9), min_val = 0, max_val = 1, wins = TRUE, trim = .01, scale_to = 10, plot = TRUE)%>%
  POMS(c(PubDvFm_1, PubDvFm_3, PubDvFm_5, PubDvFm_7, PubDvFm_9, PubDvFm_11, PubDvMl_1, PubDvMl_3, PubDvMl_5, PubDvMl_7, PubDvMl_9, PubDvMl_11), min_val = 1, max_val = 4, wins = FALSE, scale_to = 10, plot = FALSE) %>%
  POMS(c(mPmon_1, mPmon_3, mPmon_5, mPmon_7, mPmon_9,mPmon_11), min_val = 1, max_val = 5, wins = TRUE, trim = .01,  scale_to = 10, plot = TRUE) %>%
  POMS(c(mPNH_5, mPNH_7, mPNH_9), min_val = 1, max_val = 5, wins = TRUE, trim = .01,  scale_to = 10, plot = TRUE)

# CREATE A SINGLE PUBERTY METRIC 
df_NEWVARS_RED_RED_POMS <- df_NEWVARS_RED_RED_POMS %>%
  mutate(PUBlev_1 = case_when(
    Y_SEX == 1 ~ P_PubDvMl_1,
    Y_SEX == 2 ~ P_PubDvFm_1,
    TRUE ~ NA_real_
  ))

df_NEWVARS_RED_RED_POMS <- df_NEWVARS_RED_RED_POMS %>%
  mutate(PUBlev_3 = case_when(
    Y_SEX == 1 ~ P_PubDvMl_3,
    Y_SEX == 2 ~ P_PubDvFm_3,
    TRUE ~ NA_real_
  ))

df_NEWVARS_RED_RED_POMS <- df_NEWVARS_RED_RED_POMS %>%
  mutate(PUBlev_5 = case_when(
    Y_SEX == 1 ~ P_PubDvMl_5,
    Y_SEX == 2 ~ P_PubDvFm_5,
    TRUE ~ NA_real_
  ))

df_NEWVARS_RED_RED_POMS <- df_NEWVARS_RED_RED_POMS %>%
  mutate(PUBlev_7 = case_when(
    Y_SEX == 1 ~ P_PubDvMl_7,
    Y_SEX == 2 ~ P_PubDvFm_7,
    TRUE ~ NA_real_
  ))

df_NEWVARS_RED_RED_POMS <- df_NEWVARS_RED_RED_POMS %>%
  mutate(PUBlev_9 = case_when(
    Y_SEX == 1 ~ P_PubDvMl_9,
    Y_SEX == 2 ~ P_PubDvFm_9,
    TRUE ~ NA_real_
  ))

df_NEWVARS_RED_RED_POMS <- df_NEWVARS_RED_RED_POMS %>%
  mutate(PUBlev_11 = case_when(
    Y_SEX == 1 ~ P_PubDvMl_11,
    Y_SEX == 2 ~ P_PubDvFm_11,
    TRUE ~ NA_real_
  ))

################################################################################
################## REDUCE DATA 

names(df_NEWVARS_RED_RED_POMS)

## SPLIT UP THE DATA 

df_NEWVARS_RED_RED_POMS_RED <- df_NEWVARS_RED_RED_POMS %>%
  select(-c( "mNBdang_1",       "mNBdang_3",      
             "mNBdang_5",       "mNBdang_7",       "mNBdang_9",       "mNBdang_11",    "mSchAdv_1",       "mSchAdv_3",       "mSchAdv_5",       "mSchAdv_7",      
             "mSchAdv_9",       "FAMCONY_1",       "FAMCONY_3",       "FAMCONY_5",       "FAMCONY_7",       "FAMCONY_9",       "FAMCONY_11",      "mLkCogEn",       
             "dtiMnMOT_1",      "dtiMnMOT_5",      "dtiMnMOT_9",      "dti_FA_CNCG_M_1", "dti_FA_CNCG_M_5", "dti_FA_CNCG_M_9", "dti_FA_ILF_M_1",  "dti_FA_ILF_M_5", 
             "dti_FA_ILF_M_9",  "dti_FA_SLF_M_1",  "dti_FA_SLF_M_5",  "dti_FA_SLF_M_9",  "dti_FA_UNC_M_1",  "dti_FA_UNC_M_5",  "dti_FA_UNC_M_9",  "dti_FA_iFSF_M_1",
             "dti_FA_iFSF_M_5", "dti_FA_iFSF_M_9", "dti_FA_pSLF_M_1", "dti_FA_pSLF_M_5", "dti_FA_pSLF_M_9", "dti_FA_tSLF_M_1", "dti_FA_tSLF_M_5", "dti_FA_tSLF_M_9",
             "dts_FA_ALL_1",    "dts_FA_ALL_5",    "dts_FA_ALL_9",    "mPmon_1",        
             "mPmon_3",         "mPmon_5",         "mPmon_7",         "mPmon_9",         "mPmon_11",        "mPNH_5",          "mPNH_7",          "mPNH_9",  
             "PubDvFm_1",     "PubDvFm_3",      "PubDvFm_5",     "PubDvFm_7",        "PubDvFm_9",          "PubDvFm_11",        
             "PubDvMl_1",          "PubDvMl_3",         "PubDvMl_5",          "PubDvMl_7",          "PubDvMl_9",         
             "PubDvMl_11",
             "P_PubDvFm_1",        "P_PubDvFm_3",        "P_PubDvFm_5",        "P_PubDvFm_7",        "P_PubDvFm_9",        "P_PubDvFm_11",       "P_PubDvMl_1",       
             "P_PubDvMl_3",        "P_PubDvMl_5",        "P_PubDvMl_7",        "P_PubDvMl_9",        "P_PubDvMl_11"
             ))



######################## SAVE 
names(df_NEWVARS_RED_RED_POMS_RED)


write.csv(df_NEWVARS_RED_RED_POMS_RED,"ABCD_HORM_METH_SEM_MEASUREMEMT_4.24.26.csv", row.names = F)
prepareMplusData(df_NEWVARS_RED_RED_POMS_RED,"ABCD_HORM_METH_SEM_MEASUREMEMT_4.24.26.dat", inpfile =T)

## DESCRIPTIVES 
write.csv(DESCRIPS, "ABCD_HORMESIS_METHODS_DESCRIPTIVES_FINAL_4.24.26.csv", row.names = F)
write.csv(qc_table_complete, "ABCD_HORMESIS_METHODS_IMAGING_QC_FINAL_4.24.26.csv", row.names = F)

