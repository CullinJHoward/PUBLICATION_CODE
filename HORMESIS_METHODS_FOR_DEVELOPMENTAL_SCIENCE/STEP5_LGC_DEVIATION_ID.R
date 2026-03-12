#Library
library(dplyr)
library(lavaan)
library(survey)
library(lavaan.survey)
library(MplusAutomation)

# Set working directory 

setwd("C:\\Users\\cjh37695\\Dropbox\\HORMESIS_METHODS\\ANALYSIS\\")

# LOAD DF 

df <- read.csv("ABCD_HORMESIS_METHODS_REDUCED_STRUC_STEP2_V2.csv")

##### FIX FAMILY ID

# get missing 
n_missing <- sum(is.na(df$FamilyID))
#fill it 
df$FamilyID[is.na(df$FamilyID)] <-
  seq(from = max(df$FamilyID, na.rm = TRUE) + 1,
      length.out = n_missing)

###############################################################################
########################## SEM MODELS 

## REWARD RESPOSIIVTY 

RWDRSP_model <- '

  # Latent growth factors

  I =~ 1*mRwRspy_1 + 1*mRwRspy_5 + 1*mRwRspy_9
  S =~ 0*mRwRspy_1 + 3*mRwRspy_5 + 5*mRwRspy_9


  # Latent means

  I ~ 1
  S ~ 1


  # Latent variances and covariance

  I ~~ I
  S ~~ S
  I ~~ S

  # Residual variances of observed outcomes

  mRwRspy_1 ~~ mRwRspy_1
  mRwRspy_5 ~~ mRwRspy_5
  mRwRspy_9 ~~ mRwRspy_9

  # Time-varying predictors

  mRwRspy_1 ~ PMC_Y_AGE_1
  mRwRspy_5 ~ PMC_Y_AGE_5
  mRwRspy_9 ~ PMC_Y_AGE_9

   
  # Time-invariant covariates predicting growth factors

  I ~ Y_SEX_EFF 
  S ~ Y_SEX_EFF 
'


# FIT MODEL

RWDRSP_fit <- growth(
  RWDRSP_model,
  data = df,
  sampling.weights = "ppensity",
  cluster = "FamilyID",
  estimator = "MLR",
  missing = "fiml"
)

## SENSATION SEEKING

SENS_model <- '

  # Latent growth factors

  I =~ 1*mSenSe_1 + 1*mSenSe_5 + 1*mSenSe_9
  S =~ 0*mSenSe_1 + 3*mSenSe_5 + 5*mSenSe_9


  # Latent means

  I ~ 1
  S ~ 1


  # Latent variances and covariance

  I ~~ I
  S ~~ S
  I ~~ S

  # Residual variances of observed outcomes

  mSenSe_1 ~~ mSenSe_1
  mSenSe_5 ~~ mSenSe_5
  mSenSe_9 ~~ mSenSe_9

  # Time-varying predictors

  mSenSe_1 ~ PMC_Y_AGE_1
  mSenSe_5 ~ PMC_Y_AGE_5
  mSenSe_9 ~ PMC_Y_AGE_9

   
  # Time-invariant covariates predicting growth factors

  I ~ Y_SEX_EFF 
  S ~ Y_SEX_EFF 
'


# FIT MODEL 
SESNS_fit <- growth(
  SENS_model,
  data = df,
  sampling.weights = "ppensity",
  cluster = "FamilyID",
  estimator = "MLR",
  missing = "fiml"
)


## PLANNING

PLAN_model <- '

  # Latent growth factors

  I =~ 1*mPlan_1 + 1*mPlan_5 + 1*mPlan_9
  S =~ 0*mPlan_1 + 3*mPlan_5 + 5*mPlan_9


  # Latent means

  I ~ 1
  S ~ 1


  # Latent variances and covariance

  I ~~ I
  S ~~ S
  I ~~ S

  # Residual variances of observed outcomes

  mPlan_1 ~~ mPlan_1
  mPlan_5 ~~ mPlan_5
  mPlan_9 ~~ mPlan_9

  # Time-varying predictors

  mPlan_1 ~ PMC_Y_AGE_1
  mPlan_5 ~ PMC_Y_AGE_5
  mPlan_9 ~ PMC_Y_AGE_9

   
  # Time-invariant covariates predicting growth factors

  I ~ Y_SEX_EFF 
  S ~ Y_SEX_EFF 
'


# FIT MODEL 
PLAN_fit <- growth(
  PLAN_model,
  data = df,
  sampling.weights = "ppensity",
  cluster = "FamilyID",
  estimator = "MLR",
  missing = "fiml"
)

## PERS

PERS_model <- '

  # Latent growth factors

  I =~ 1*mPers_1 + 1*mPers_5 + 1*mPers_9
  S =~ 0*mPers_1 + 3*mPers_5 + 5*mPers_9


  # Latent means

  I ~ 1
  S ~ 1


  # Latent variances and covariance

  I ~~ I
  S ~~ S
  I ~~ S

  # Residual variances of observed outcomes

  mPers_1 ~~ mPers_1
  mPers_5 ~~ mPers_5
  mPers_9 ~~ mPers_9

  # Time-varying predictors

  mPers_1 ~ PMC_Y_AGE_1
  mPers_5 ~ PMC_Y_AGE_5
  mPers_9 ~ PMC_Y_AGE_9

   
  # Time-invariant covariates predicting growth factors

  I ~ Y_SEX_EFF 
  S ~ Y_SEX_EFF 
'


# FIT MODEL 
PERS_fit <- growth(
  PERS_model,
  data = df,
  sampling.weights = "ppensity",
  cluster = "FamilyID",
  estimator = "MLR",
  missing = "fiml"
)

## ReAp

ReAp_model <- '

  # Latent growth factors

  I =~ 1*mReApp_7 + 1*mReApp_9 + 1*mReApp_11
  S =~ 0*mReApp_7 + 1*mReApp_9 + 2*mReApp_11


  # Latent means

  I ~ 1
  S ~ 1


  # Latent variances and covariance

  I ~~ I
  S ~~ S
  I ~~ S

  # Residual variances of observed outcomes

  mReApp_7 ~~ mReApp_7
  mReApp_9 ~~ mReApp_9
  mReApp_11 ~~ mReApp_11

  # Time-varying predictors

  mReApp_7 ~ PMC_Y_AGE_7
  mReApp_9 ~ PMC_Y_AGE_9
  mReApp_11 ~ PMC_Y_AGE_11

   
  # Time-invariant covariates predicting growth factors

  I ~ Y_SEX_EFF 
  S ~ Y_SEX_EFF 
'


# FIT MODEL 
ReAp_fit <- growth(
  ReAp_model,
  data = df,
  sampling.weights = "ppensity",
  cluster = "FamilyID",
  estimator = "MLR",
  missing = "fiml"
)

## EmoSup

EmoSup_model <- '

  # Latent growth factors

  I =~ 1*mEmoSup_7 + 1*mEmoSup_9 + 1*mEmoSup_11
  S =~ 0*mEmoSup_7 + 1*mEmoSup_9 + 2*mEmoSup_11


  # Latent means

  I ~ 1
  S ~ 1


  # Latent variances and covariance

  I ~~ I
  S ~~ S
  I ~~ S

  # Residual variances of observed outcomes

  mEmoSup_7 ~~ mEmoSup_7
  mEmoSup_9 ~~ mEmoSup_9
  mEmoSup_11 ~~ mEmoSup_11

  # Time-varying predictors

  mEmoSup_7 ~ PMC_Y_AGE_7
  mEmoSup_9 ~ PMC_Y_AGE_9
  mEmoSup_11 ~ PMC_Y_AGE_11

   
  # Time-invariant covariates predicting growth factors

  I ~ Y_SEX_EFF 
  S ~ Y_SEX_EFF 
'


# FIT MODEL 
EmoSup_fit <- growth(
  EmoSup_model,
  data = df,
  sampling.weights = "ppensity",
  cluster = "FamilyID",
  estimator = "MLR",
  missing = "fiml"
)

## Cata

Cata_model <- '

  # Latent growth factors

  I =~ 1*mCata_7 + 1*mCata_9 + 1*mCata_11
  S =~ 0*mCata_7 + 1*mCata_9 + 2*mCata_11


  # Latent means

  I ~ 1
  S ~ 1


  # Latent variances and covariance

  I ~~ I
  S ~~ S
  I ~~ S

  # Residual variances of observed outcomes

  mCata_7 ~~ mCata_7
  mCata_9 ~~ mCata_9
  mCata_11 ~~ mCata_11

  # Time-varying predictors

  mCata_7 ~ PMC_Y_AGE_7
  mCata_9 ~ PMC_Y_AGE_9
  mCata_11 ~ PMC_Y_AGE_11

   
  # Time-invariant covariates predicting growth factors

  I ~ Y_SEX_EFF 
  S ~ Y_SEX_EFF 
'


# FIT MODEL 
Cata_fit <- growth(
  Cata_model,
  data = df,
  sampling.weights = "ppensity",
  cluster = "FamilyID",
  estimator = "MLR",
  missing = "fiml"
)

## NegSec

NegSec_model <- '

  # Latent growth factors

  I =~ 1*mNgSec_7 + 1*mNgSec_9 + 1*mNgSec_11
  S =~ 0*mNgSec_7 + 1*mNgSec_9 + 2*mNgSec_11


  # Latent means

  I ~ 1
  S ~ 1


  # Latent variances and covariance

  I ~~ I
  S ~~ S
  I ~~ S

  # Residual variances of observed outcomes

  mNgSec_7 ~~ mNgSec_7
  mNgSec_9 ~~ mNgSec_9
  mNgSec_11 ~~ mNgSec_11

  # Time-varying predictors

  mNgSec_7 ~ PMC_Y_AGE_7
  mNgSec_9 ~ PMC_Y_AGE_9
  mNgSec_11 ~ PMC_Y_AGE_11

   
  # Time-invariant covariates predicting growth factors

  I ~ Y_SEX_EFF 
  S ~ Y_SEX_EFF 
'


# FIT MODEL 
NegSec_fit <- growth(
  NegSec_model,
  data = df,
  sampling.weights = "ppensity",
  cluster = "FamilyID",
  estimator = "MLR",
  missing = "fiml"
)

## Dist

Dist_model <- '

  # Latent growth factors

  I =~ 1*mDist_7 + 1*mDist_9 + 1*mDist_11
  S =~ 0*mDist_7 + 1*mDist_9 + 2*mDist_11


  # Latent means

  I ~ 1
  S ~ 1


  # Latent variances and covariance

  I ~~ I
  S ~~ S
  I ~~ S

  # Residual variances of observed outcomes

  mDist_7 ~~ mDist_7
  mDist_9 ~~ mDist_9
  mDist_11 ~~ mDist_11

  # Time-varying predictors

  mDist_7 ~ PMC_Y_AGE_7
  mDist_9 ~ PMC_Y_AGE_9
  mDist_11 ~ PMC_Y_AGE_11

   
  # Time-invariant covariates predicting growth factors

  I ~ Y_SEX_EFF 
  S ~ Y_SEX_EFF 
'


# FIT MODEL 
Dist_fit <- growth(
  Dist_model,
  data = df,
  sampling.weights = "ppensity",
  cluster = "FamilyID",
  estimator = "MLR",
  missing = "fiml"
)

## Attune

Attune_model <- '

  # Latent growth factors

  I =~ 1*mAttune_7 + 1*mAttune_9 + 1*mAttune_11
  S =~ 0*mAttune_7 + 1*mAttune_9 + 2*mAttune_11


  # Latent means

  I ~ 1
  S ~ 1


  # Latent variances and covariance

  I ~~ I
  S ~~ S
  I ~~ S

  # Residual variances of observed outcomes

  mAttune_7 ~~ mAttune_7
  mAttune_9 ~~ mAttune_9
  mAttune_11 ~~ mAttune_11

  # Time-varying predictors

  mAttune_7 ~ PMC_Y_AGE_7
  mAttune_9 ~ PMC_Y_AGE_9
  mAttune_11 ~ PMC_Y_AGE_11

   
  # Time-invariant covariates predicting growth factors

  I ~ Y_SEX_EFF 
  S ~ Y_SEX_EFF 
'


# FIT MODEL 
Attune_fit <- growth(
  Attune_model,
  data = df,
  sampling.weights = "ppensity",
  cluster = "FamilyID",
  estimator = "MLR",
  missing = "fiml"
)

summary(EmoSup_fit, standardized = TRUE, fit.measures = TRUE, ci = TRUE)
modindices(fit, sort. = TRUE)
parameterEstimates(fit, standardized = TRUE, ci = TRUE)


# ----------------------------
# function to summarize one LGC
# ----------------------------

summarize_lgc <- function(fit,
                          construct_name = "UNKNOWN",
                          growth_terms,
                          intercept_term = growth_terms[1],
                          slope_term = if (length(growth_terms) >= 2) growth_terms[2],
                          summary_df_name = "summary_df") {
  
  # -----------------------------
  # checks
  # -----------------------------
  if (missing(growth_terms) || length(growth_terms) == 0) {
    stop("You must provide growth_terms, e.g., c('I', 'S').")
  }
  
  if (!inherits(fit, "lavaan")) {
    stop("fit must be a lavaan object.")
  }
  
  pe <- lavaan::parameterEstimates(fit, standardized = FALSE, ci = TRUE)
  
  # -----------------------------
  # fit indices
  # -----------------------------
  fm <- lavaan::fitMeasures(
    fit,
    c("chisq.scaled", "df.scaled", "pvalue.scaled",
      "rmsea.scaled", "cfi.scaled", "tli.scaled", "srmr")
  )
  
  # helper to pull one parameter row
  get_param <- function(lhs, op, rhs = "") {
    row <- pe[pe$lhs == lhs & pe$op == op & pe$rhs == rhs, , drop = FALSE]
    if (nrow(row) == 0) {
      return(list(est = NA_real_, p = NA_real_))
    }
    list(est = row$est[1], p = row$pvalue[1])
  }
  
  # -----------------------------
  # extract means and variances
  # -----------------------------
  out_list <- list(
    Construct = construct_name,
    ChiSq = unname(fm["chisq.scaled"]),
    DF = unname(fm["df.scaled"]),
    ChiSq_p = unname(fm["pvalue.scaled"]),
    RMSEA = unname(fm["rmsea.scaled"]),
    CFI = unname(fm["cfi.scaled"]),
    TLI = unname(fm["tli.scaled"]),
    SRMR = unname(fm["srmr"])
  )
  
  for (gt in growth_terms) {
    mean_par <- get_param(gt, "~1")
    var_par  <- get_param(gt, "~~", gt)
    
    out_list[[paste0(gt, "_Est")]]   <- mean_par$est
    out_list[[paste0(gt, "_p")]]     <- mean_par$p
    out_list[[paste0(gt, "_Var")]]   <- var_par$est
    out_list[[paste0(gt, "_Var_p")]] <- var_par$p
  }
  
  # -----------------------------
  # dynamic wave-specific predictions
  # -----------------------------
  int_rows <- pe[pe$lhs == intercept_term & pe$op == "=~", , drop = FALSE]
  
  if (nrow(int_rows) == 0) {
    stop("Could not find measurement rows for intercept_term = ", intercept_term)
  }
  
  observed_vars <- int_rows$rhs
  intercept_est <- out_list[[paste0(intercept_term, "_Est")]]
  
  # extract numeric wave suffix from variable name
  get_wave_suffix <- function(x) {
    out <- sub(".*_([0-9]+)$", "\\1", x)
    ifelse(grepl("^[0-9]+$", out), out, x)
  }
  
  if (!is.null(slope_term) && !is.na(slope_term)) {
    slope_rows <- pe[pe$lhs == slope_term & pe$op == "=~", , drop = FALSE]
    
    if (nrow(slope_rows) == 0) {
      stop("Could not find measurement rows for slope_term = ", slope_term)
    }
    
    slope_loadings <- slope_rows$est[match(observed_vars, slope_rows$rhs)]
    slope_est <- out_list[[paste0(slope_term, "_Est")]]
    pred_vals <- intercept_est + slope_est * slope_loadings
  } else {
    slope_loadings <- rep(0, length(observed_vars))
    pred_vals <- rep(intercept_est, length(observed_vars))
  }
  
  for (j in seq_along(observed_vars)) {
    wave_suffix <- get_wave_suffix(observed_vars[j])
    out_list[[paste0("Wave_", wave_suffix, "_Mean")]] <- pred_vals[j]
    out_list[[paste0("Wave_", wave_suffix, "_Loading")]] <- slope_loadings[j]
  }
  
  # one-row data frame
  out <- as.data.frame(out_list, check.names = FALSE)
  rownames(out) <- NULL
  
  # -----------------------------
  # append to existing summary_df or create new one
  # -----------------------------
  if (exists(summary_df_name, envir = .GlobalEnv)) {
    old_df <- get(summary_df_name, envir = .GlobalEnv)
    
    all_cols <- union(names(old_df), names(out))
    
    for (nm in setdiff(all_cols, names(old_df))) {
      old_df[[nm]] <- NA
    }
    for (nm in setdiff(all_cols, names(out))) {
      out[[nm]] <- NA
    }
    
    old_df <- old_df[, all_cols, drop = FALSE]
    out    <- out[, all_cols, drop = FALSE]
    
    summary_df <- rbind(old_df, out)
  } else {
    summary_df <- out
  }
  
  assign(summary_df_name, summary_df, envir = .GlobalEnv)
  
  return(out)
}

#
summarize_lgc(
  fit = PERS_fit,
  construct_name = "Pers",
  growth_terms = c("I", "S")
)
summarize_lgc(
  fit = PLAN_fit,
  construct_name = "Plan",
  growth_terms = c("I", "S")
)
summarize_lgc(
  fit = RWDRSP_fit,
  construct_name = "RwRsp",
  growth_terms = c("I", "S")
)
summarize_lgc(
  fit = SESNS_fit,
  construct_name = "SenSee",
  growth_terms = c("I", "S")
)
summarize_lgc(
  fit = Attune_fit,
  construct_name = "Attun",
  growth_terms = c("I", "S")
)
summarize_lgc(
  fit = Cata_fit,
  construct_name = "Cata",
  growth_terms = c("I", "S")
)
summarize_lgc(
  fit = Dist_fit,
  construct_name = "Dist",
  growth_terms = c("I", "S")
)
summarize_lgc(
  fit = EmoSup_fit,
  construct_name = "EmoSup",
  growth_terms = c("I", "S")
)
summarize_lgc(
  fit = NegSec_fit,
  construct_name = "NegSec",
  growth_terms = c("I", "S")
)
summarize_lgc(
  fit = ReAp_fit,
  construct_name = "ReApp",
  growth_terms = c("I", "S")
)


## Reorganize output for later differencing
## Add expected trajectory means directly into df

df$SP_Pers_1 <- summary_df$Wave_1_Mean[summary_df$Construct == "Pers"]
df$SP_Pers_5 <- summary_df$Wave_5_Mean[summary_df$Construct == "Pers"]
df$SP_Pers_9 <- summary_df$Wave_9_Mean[summary_df$Construct == "Pers"]

df$SP_Plan_1 <- summary_df$Wave_1_Mean[summary_df$Construct == "Plan"]
df$SP_Plan_5 <- summary_df$Wave_5_Mean[summary_df$Construct == "Plan"]
df$SP_Plan_9 <- summary_df$Wave_9_Mean[summary_df$Construct == "Plan"]

df$SP_RwRsp_1 <- summary_df$Wave_1_Mean[summary_df$Construct == "RwRsp"]
df$SP_RwRsp_5 <- summary_df$Wave_5_Mean[summary_df$Construct == "RwRsp"]
df$SP_RwRsp_9 <- summary_df$Wave_9_Mean[summary_df$Construct == "RwRsp"]

df$SP_SenSee_1 <- summary_df$Wave_1_Mean[summary_df$Construct == "SenSee"]
df$SP_SenSee_5 <- summary_df$Wave_5_Mean[summary_df$Construct == "SenSee"]
df$SP_SenSee_9 <- summary_df$Wave_9_Mean[summary_df$Construct == "SenSee"]

df$SP_Attun_7  <- summary_df$Wave_7_Mean[summary_df$Construct == "Attun"]
df$SP_Attun_9  <- summary_df$Wave_9_Mean[summary_df$Construct == "Attun"]
df$SP_Attun_11 <- summary_df$Wave_11_Mean[summary_df$Construct == "Attun"]

df$SP_Cata_7  <- summary_df$Wave_7_Mean[summary_df$Construct == "Cata"]
df$SP_Cata_9  <- summary_df$Wave_9_Mean[summary_df$Construct == "Cata"]
df$SP_Cata_11 <- summary_df$Wave_11_Mean[summary_df$Construct == "Cata"]

df$SP_Dist_7  <- summary_df$Wave_7_Mean[summary_df$Construct == "Dist"]
df$SP_Dist_9  <- summary_df$Wave_9_Mean[summary_df$Construct == "Dist"]
df$SP_Dist_11 <- summary_df$Wave_11_Mean[summary_df$Construct == "Dist"]

df$SP_EmoSup_7  <- summary_df$Wave_7_Mean[summary_df$Construct == "EmoSup"]
df$SP_EmoSup_9  <- summary_df$Wave_9_Mean[summary_df$Construct == "EmoSup"]
df$SP_EmoSup_11 <- summary_df$Wave_11_Mean[summary_df$Construct == "EmoSup"]

df$SP_NegSec_7  <- summary_df$Wave_7_Mean[summary_df$Construct == "NegSec"]
df$SP_NegSec_9  <- summary_df$Wave_9_Mean[summary_df$Construct == "NegSec"]
df$SP_NegSec_11 <- summary_df$Wave_11_Mean[summary_df$Construct == "NegSec"]

df$SP_ReApp_7  <- summary_df$Wave_7_Mean[summary_df$Construct == "ReApp"]
df$SP_ReApp_9  <- summary_df$Wave_9_Mean[summary_df$Construct == "ReApp"]
df$SP_ReApp_11 <- summary_df$Wave_11_Mean[summary_df$Construct == "ReApp"]

## CREATING RESIDUALS 

df$RSpers1 <- (df$mPers_1 - df$SP_Pers_1)
df$RSpers5 <- (df$mPers_5 - df$SP_Pers_5)
df$RSpers9 <- (df$mPers_9 - df$SP_Pers_9)          
df$RSplan1 <- (df$mPlan_1 - df$SP_Plan_1)   
df$RSplan5 <- (df$mPlan_5 - df$SP_Plan_5)   
df$RSplan9 <- (df$mPlan_9 - df$SP_Plan_9)                
df$RSsense1 <- (df$mSenSe_1 - df$SP_SenSee_1)  
df$RSsense5 <- (df$mSenSe_5 - df$SP_SenSee_5)  
df$RSsense9 <- (df$mSenSe_9 - df$SP_SenSee_9)  
df$RSrwrspy1 <- (df$mRwRspy_1 - df$SP_RwRsp_1)  
df$RSrwrspy5 <- (df$mRwRspy_5 - df$SP_RwRsp_5)  
df$RSrwrspy9 <- (df$mRwRspy_9 - df$SP_RwRsp_9)  
df$RSattune7 <- (df$mAttune_7 - df$SP_Attun_7)
df$RSattune9 <- (df$mAttune_9 - df$SP_Attun_9)
df$RSattune11 <- (df$mAttune_11 - df$SP_Attun_11)
df$RSdist7 <- (df$mDist_7 - df$SP_Dist_7)
df$RSdist9 <- (df$mDist_9 - df$SP_Dist_9)
df$RSdist11 <- (df$mDist_11 - df$SP_Dist_11)
df$RSngsec7 <- (df$mNgSec_7 - df$SP_NegSec_7)
df$RSngsec9 <- (df$mNgSec_9 - df$SP_NegSec_9)
df$RSngsec11  <- (df$mNgSec_11 - df$SP_NegSec_11)        
df$RScata7 <- (df$mCata_7 - df$SP_Cata_7)
df$RScata9 <- (df$mCata_9 - df$SP_Cata_9)
df$RScata11 <- (df$mCata_11 - df$SP_Cata_11)
df$RSemosup7 <- (df$mEmoSup_7 - df$SP_EmoSup_7)
df$RSemosup9  <- (df$mEmoSup_9 - df$SP_EmoSup_9)         
df$RSemosup11 <- (df$mEmoSup_11 - df$SP_EmoSup_11)
df$RSreapp7 <- (df$mReApp_7 - df$SP_ReApp_7)
df$RSreapp9 <- (df$mReApp_9 - df$SP_ReApp_9)
df$RSreapp11 <- (df$mReApp_11 - df$SP_ReApp_11)


###############################################################################
################### MERGE AND SAVE 

#ONLY COMPLETE 
write.csv(df, "ABCD_HORMESIS_METHODS_REDUCED_STRUC_STEP3.csv", row.names = FALSE)

prepareMplusData(df, "ABCD_HORMESIS_METHODS_REDUCED_STRUC_STEP3.dat")

