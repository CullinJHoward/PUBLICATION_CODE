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

df <- read.csv("ABCD_HORM_METH_PREP_12.31.25.csv")

options(max.print = 10000)

################################################################################
###################### SCALE COMPUTATION AND LONGITUDINAL CHANGE VISUALIZATION


## SET UP A FUNCTION TO HELP COMPUTE

make_scale_means <- function(data, item_sets, new_names = names(item_sets),
                             max_missing_prop = 0.80, digits = 3,
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
    
    # Missing counts for row scoring
    n_missing_items <- rowSums(is.na(scale_df))
    prop_missing <- n_missing_items / n_items
    
    all_missing <- n_missing_items == n_items
    partial_rule <- prop_missing > max_missing_prop & !all_missing
    
    # Compute scale mean
    scale_mean <- rowMeans(scale_df, na.rm = TRUE)
    scale_mean[all_missing] <- NA
    scale_mean[partial_rule] <- NA
    scale_mean[is.nan(scale_mean)] <- NA
    
    df[[new_var]] <- scale_mean
    
    # Descriptive Ns
    N_valid <- sum(!is.na(scale_mean))
    N_missing_total <- sum(is.na(scale_mean))
    N_missing_all_items <- sum(all_missing)
    N_missing_partial_rule <- sum(partial_rule)
    
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
    
    # Reliability flag: X if all available reliability estimates are < .70
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
      N_missing_partial_rule = N_missing_partial_rule,
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

names(df)

item_sets <- list(
  mHOMEsf_11 = c("ErlHs1P_11","ErlHs2P_11","ErlHs3P_11","ErlHs4P_11","ErlHs5P_11",
               "ErlHs6P_11","ErlHs7P_11","ErlHs8P_11","ErlHs9P_11","ErlHs10P_11"),
  mFconP_1 = c("Fcon1P_1","Fcon2rP_1","Fcon3P_1","Fcon4rP_1","Fcon5P_1","Fcon6P_1",
               "Fcon7rP_1" ,"Fcon8P_1" ,"Fcon8P.1_1"),
  mFconP_3 = c("Fcon1P_3","Fcon2rP_3","Fcon3P_3","Fcon4rP_3","Fcon5P_3","Fcon6P_3",
               "Fcon7rP_3" ,"Fcon8P_3" ,"Fcon8P.1_3"),
  mFconP_5 = c("Fcon1P_5","Fcon2rP_5","Fcon3P_5","Fcon4rP_5","Fcon5P_5","Fcon6P_5",
               "Fcon7rP_5" ,"Fcon8P_5" ,"Fcon8P.1_5"),
  mFconP_7 = c("Fcon1P_7","Fcon2rP_7","Fcon3P_7","Fcon4rP_7","Fcon5P_7","Fcon6P_7",
               "Fcon7rP_7" ,"Fcon8P_7" ,"Fcon8P.1_7"),
  mFconP_9 = c("Fcon1P_9","Fcon2rP_9","Fcon3P_9","Fcon4rP_9","Fcon5P_9","Fcon6P_9",
               "Fcon7rP_9" ,"Fcon8P_9" ,"Fcon8P.1_9"),
  mFconP_11 = c("Fcon1P_11","Fcon2rP_11","Fcon3P_11","Fcon4rP_11","Fcon5P_11","Fcon6P_11",
               "Fcon7rP_11" ,"Fcon8P_11" ,"Fcon8P.1_11"),
  mFconY_1 = c("Fcon1Y_1","Fcon2rY_1","Fcon3Y_1","Fcon4rY_1","Fcon5Y_1","Fcon6Y_1",
               "Fcon7rY_1" ,"Fcon8Y_1" ,"Fcon9Y_1"),
  mFconY_3 = c("Fcon1Y_3","Fcon2rY_3","Fcon3Y_3","Fcon4rY_3","Fcon5Y_3","Fcon6Y_3",
               "Fcon7rY_3" ,"Fcon8Y_3" ,"Fcon9Y_3"),
  mFconY_5 = c("Fcon1Y_5","Fcon2rY_5","Fcon3Y_5","Fcon4rY_5","Fcon5Y_5","Fcon6Y_5",
               "Fcon7rY_5" ,"Fcon8Y_5" ,"Fcon9Y_5"),
  mFconY_7 = c("Fcon1Y_7","Fcon2rY_7","Fcon3Y_7","Fcon4rY_7","Fcon5Y_7","Fcon6Y_7",
               "Fcon7rY_7" ,"Fcon8Y_7" ,"Fcon9Y_7"),
  mFconY_9 = c("Fcon1Y_9","Fcon2rY_9","Fcon3Y_9","Fcon4rY_9","Fcon5Y_9","Fcon6Y_9",
               "Fcon7rY_9" ,"Fcon8Y_9" ,"Fcon9Y_9"),
  mFconY_11 = c("Fcon1Y_11","Fcon2rY_11","Fcon3Y_11","Fcon4rY_11","Fcon5Y_11","Fcon6Y_11",
                "Fcon7rY_11" ,"Fcon8Y_11" ,"Fcon9Y_11"),
  mNBsaf_1 = c("NBHsaf1P_1", "NBHsaf2P_1", "NBHsaf3P_1"),
  mNBsaf_3 = c("NBHsaf1P_3", "NBHsaf2P_3", "NBHsaf3P_3"),
  mNBsaf_5 = c("NBHsaf1P_5", "NBHsaf2P_5", "NBHsaf3P_5"),
  mNBsaf_7 = c("NBHsaf1P_7", "NBHsaf2P_7", "NBHsaf3P_7"),
  mNBsaf_9 = c("NBHsaf1P_9", "NBHsaf2P_9", "NBHsaf3P_9"),
  mNBsaf_11 = c("NBHsaf1P_11", "NBHsaf2P_11", "NBHsaf3P_11"),
  mPedsup_7 = c("PNGedu1Y_7", "PNGedu2Y_7", "PNGedu3Y_7"),
  mPsupr_7 = c("PNGsup1Y_7", "PNGsup2Y_7", "PNGsup3Y_7", "PNGsup4Y_7", "PNGsup5Y_7"),
  mSchEnv_1 = c("SCHenv1Y_1", "SCHenv2Y_1", "SCHenv3Y_1", "SCHenv4Y_1", "SCHenv5Y_1", "SCHenv6Y_1"),
  mSchEnv_3 = c("SCHenv1Y_3", "SCHenv2Y_3", "SCHenv3Y_3", "SCHenv4Y_3", "SCHenv5Y_3", "SCHenv6Y_3"),
  mSchEnv_5 = c("SCHenv1Y_5", "SCHenv2Y_5", "SCHenv3Y_5", "SCHenv4Y_5", "SCHenv5Y_5", "SCHenv6Y_5"),
  mSchEnv_7 = c("SCHenv1Y_7", "SCHenv2Y_7", "SCHenv3Y_7", "SCHenv4Y_7", "SCHenv5Y_7", "SCHenv6Y_7"),
  mSchEnv_9 = c("SCHenv1Y_9", "SCHenv2Y_9", "SCHenv3Y_9", "SCHenv4Y_9", "SCHenv5Y_9", "SCHenv6Y_9"),
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
  mReApp_11 = c("yReAp1Y_11", "yReAp2Y_11", "yReAp3Y_11")
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
"FAMCONP1",      "FAMCONP3",      "FAMCONP5",      "FAMCONP7",      "FAMCONP9",      "FAMCONP11",     "FAMCONY1",      "FAMCONY3", 
"FAMCONY5",      "FAMCONY7",     "FAMCONY9",      "FAMCONY11",      "ErlHs1PR",      "ErlHs2PR",     
"ErlHs3PR",      "ErlHs4PR",      "ErlHs5PR",      "ErlHs6PR",      "NBHsaf1P_1R",   "NBHsaf2P_1R",   "NBHsaf3P_1R",   "NBHsaf1P_3R",  
"NBHsaf2P_3R",   "NBHsaf3P_3R",   "NBHsaf1P_5R",   "NBHsaf2P_5R",   "NBHsaf3P_5R",   "NBHsaf1P_7R",   "NBHsaf2P_7R",   "NBHsaf3P_7R",  
"NBHsaf1P_9R",   "NBHsaf2P_9R",   "NBHsaf3P_9R",   "NBHsaf1P_11R",  "NBHsaf2P_11R",  "NBHsaf3P_11R",  "SCHenv1Y_1R",   "SCHenv2Y_1R",  
"SCHenv3Y_1R",   "SCHenv4Y_1R",   "SCHenv5Y_1R",   "SCHenv6Y_1R",   "SCHenv1Y_3R",   "SCHenv2Y_3R",   "SCHenv3Y_3R",   "SCHenv4Y_3R",  
"SCHenv5Y_3R",   "SCHenv6Y_3R",   "SCHenv1Y_5R",   "SCHenv2Y_5R",   "SCHenv3Y_5R",   "SCHenv4Y_5R",   "SCHenv5Y_5R",   "SCHenv6Y_5R",  
"SCHenv1Y_7R",   "SCHenv2Y_7R",   "SCHenv3Y_7R",   "SCHenv4Y_7R",   "SCHenv5Y_7R",   "SCHenv6Y_7R",   "SCHenv1Y_9R",   "SCHenv2Y_9R",  
"SCHenv3Y_9R",   "SCHenv4Y_9R",   "SCHenv5Y_9R",   "SCHenv6Y_9R",   "SocMob_1R",     "HiParEdu_1R",   "FAMCONP1",      "FAMCONP3",     
"FAMCONP5",      "FAMCONP7",      "FAMCONP9",      "FAMCONP11",     "FAMCONY1",      "FAMCONY3",      "FAMCONY5",      "FAMCONY7",     
"FAMCONY9",      "FAMCONY11"))


## REVERSE SCORE SOME SCALES 

reverse_items <- function(data, items, prefix = "R_") {
  
  df <- data
  
  for (v in items) {
    
    if (!v %in% names(df)) {
      stop(paste("Variable not found:", v))
    }
    
    min_val <- min(df[[v]], na.rm = TRUE)
    max_val <- max(df[[v]], na.rm = TRUE)
    
    new_name <- paste0(prefix, v)
    
    df[[new_name]] <- (max_val + min_val) - df[[v]]
    
  }
  
  return(df)
}

## APPLIED 
names(df_NEWVARS_RED)

reverse_vars <- c(
  ## HOME SCALE 
  "mHOMEsf_11",
  ## NEIGHBORHOOD SAFETY 
  "mNBsaf_1",
  "mNBsaf_3",
  "mNBsaf_5",
  "mNBsaf_7",
  "mNBsaf_9",
  "mNBsaf_11",
  ## SCHOOL ENVIRONMENT 
  "mSchEnv_1",
  "mSchEnv_3",
  "mSchEnv_5",
  "mSchEnv_7",    
  "mSchEnv_9",
  ## SOCIAL MOBILITY
  "SocMob_1",
  ## EDUCATION 
  "HiParEdu_1",
  "HiParEdu_3",
  "HiParEdu_5",
  "HiParEdu_7",
  "HiParEdu_9",
  "HiParEdu_11"
)

df_NEWVARS_RED <- reverse_items(df_NEWVARS_RED, reverse_vars)

hist(df_NEWVARS_RED$R_mSchEnv_1)
hist(df_NEWVARS_RED$mSchEnv_1)


## CENTER AND PERSON-MEAN CENTER STUFF

names(df_NEWVARS_RED)

## WITHIN PERSON CENTER AGE 
age_vars <- grep("^Y_AGE_\\d+$", names(df_NEWVARS_RED), value = TRUE)

df_NEWVARS_RED <- df_NEWVARS_RED %>%
  rowwise() %>%
  mutate(
    AGE_PM = mean(c_across(all_of(age_vars)), na.rm = TRUE)
  ) %>%
  ungroup() %>%
  mutate(
    across(
      all_of(age_vars),
      ~ .x - AGE_PM,
      .names = "PMC_{.col}"
    )
  )

# DEMOGRAPHICS 

df_NEWVARS_RED$C_INCOME6L_1 <- (df_NEWVARS_RED$INCOME6L_1 - mean(df_NEWVARS_RED$INCOME6L_1, na.rm = TRUE))
df_NEWVARS_RED$C_HiParEdu_1<- (df_NEWVARS_RED$HiParEdu_1 - mean(df_NEWVARS_RED$HiParEdu_1, na.rm = TRUE))

## SPLIT UP THE DATA 

# ADD A NUMERIC ID 

df_NEWVARS_RED$NUMID <- row.names(df_NEWVARS_RED)

SEM_DF <- df_NEWVARS_RED %>%
  select(c("subID", "SiteID",  "FamilyID", "Y_HISP", "Y_SEX", "ppensity", "NUMID",
           "HiParEdu_1", "C_HiParEdu_1",  "INCOME6L_1", "C_INCOME6L_1", "MarWrkSt_1",
           "PRBlPov_1", "PRCrowd_1", "PRLowEd_1",  "PRUnEmp_1", "R_SocMob_1",
           "Y_AGE_1", "Y_AGE_3", "Y_AGE_5", "Y_AGE_7", "Y_AGE_9", "Y_AGE_11",
           "PMC_Y_AGE_1", "PMC_Y_AGE_3", "PMC_Y_AGE_5", "PMC_Y_AGE_7", "PMC_Y_AGE_9",
           "PMC_Y_AGE_11", "mHOMEsf_11",    "mFconP_1",      "mFconP_3",     
           "mFconP_5",      "mFconP_7",      "mFconP_9",      "mFconP_11",     
           "mFconY_1",      "mFconY_3",      "mFconY_5",      "mFconY_7",      
           "mFconY_9",      "mFconY_11",     "mNBsaf_1",     
           "mNBsaf_3",      "mNBsaf_5",      "mNBsaf_7",      "mNBsaf_9",      
           "mNBsaf_11",     "mPedsup_7",     "mPsupr_7",      "mSchEnv_1",     
           "mSchEnv_3",     "mSchEnv_5",     "mSchEnv_7",    
           "mSchEnv_9",     "mNgUrg_1",      "mNgUrg_5",      "mNgUrg_9",      
           "mPers_1",       "mPers_5",       "mPers_9",       "mPlan_1",       
           "mPlan_5",       "mPlan_9",       "mPsUrg_1",     
           "mPsUrg_5",      "mPsUrg_9",      "mSenSe_1",      "mSenSe_5",      
           "mSenSe_9",      "mBISy_1",       "mBISy_5",       "mBISy_9",       
           "mDRVy_1",       "mDRVy_5",       "mDRVy_9",      
           "mFnSky_1",      "mFnSky_5",      "mFnSky_9",      "mRwRspy_1",     
           "mRwRspy_5",     "mRwRspy_9",     "mAttune_7",     "mAttune_9",     
           "mAttune_11",    "mDist_7",       "mDist_9",      
           "mDist_11",      "mNgSec_7",      "mNgSec_9",      "mNgSec_11",     
           "mCata_7",       "mCata_9",       "mCata_11",      "mEmoSup_7",     
           "mEmoSup_9",     "mEmoSup_11",    "mReApp_7",     
           "mReApp_9",      "mReApp_11",     "R_mHOMEsf_11",  "R_mNBsaf_1",    
           "R_mNBsaf_3",    "R_mNBsaf_5",    "R_mNBsaf_7",    "R_mNBsaf_9",    
           "R_mNBsaf_11",   "R_mSchEnv_1",   "R_mSchEnv_3",  
           "R_mSchEnv_5",   "R_mSchEnv_7",   "R_mSchEnv_9",      
           "R_HiParEdu_1",  "R_HiParEdu_3",  "R_HiParEdu_5",  "R_HiParEdu_7",  
           "R_HiParEdu_9",  "R_HiParEdu_11"))

######################## SAVE 

#FULL
write.csv(df_NEWVARS_RED,"ABCD_HORMESIS_METHODS_FULL_3.9.26.csv", row.names = F)
prepareMplusData(df_NEWVARS_RED,"ABCD_HORMESIS_METHODS_FULL_3.9.26.dat", inpfile =T)
## REDUCED SEM 

write.csv(SEM_DF,"ABCD_HORMESIS_METHODS_REDUCED_3.9.26.csv", row.names = F)
prepareMplusData(SEM_DF,"ABCD_HORMESIS_METHODS_REDUCED_3.9.26.dat", inpfile =T)

## DESCRIPTIVES 
write.csv(DESCRIPS, "ABCD_HORMESIS_METHODS_SCALE_DESCRIPTIVES_3.9.26.csv", row.names = F)
