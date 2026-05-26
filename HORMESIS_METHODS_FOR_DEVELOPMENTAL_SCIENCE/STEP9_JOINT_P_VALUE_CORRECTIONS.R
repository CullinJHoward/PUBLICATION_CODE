# library
library(dplyr)


#set up environment 

setwd("C:\\Users\\cjh37695\\Dropbox\\HORMESIS_METHODS\\ANALYSIS\\")

df <- read.csv("lgc_p_value_correction_table.csv")

## split up contrast families 

CR <- subset(df, Famly == "Reappraisal ")

ES <- subset(df, Famly == "Suppression")

head(CR)


## use Cauchy approach for combinig p-value probabilties

# -------------------------------
# Cauchy combination test function
# -------------------------------

cauchy_combine <- function(p) {
  
  p <- p[!is.na(p)]
  
  if (length(p) == 0) {
    return(NA_real_)
  }
  
  if (any(p <= 0 | p >= 1)) {
    stop("All p-values must be strictly between 0 and 1 for the Cauchy combination test.")
  }
  
  T_stat <- mean(tan((0.5 - p) * pi))
  p_combined <- 0.5 - atan(T_stat) / pi
  
  return(p_combined)
}


# -----------------------------------------
# Create joint S1/S2 p-values + BH correction
# -----------------------------------------

make_joint_cauchy_df <- function(df, family_var = "Famly") {
  
  # Clean and prepare data
  df_clean <- df %>%
    mutate(
      Observed.p = as.numeric(Observed.p),
      
      # Mplus often prints p < .001 as 0.000.
      # This approximates printed zeroes as .0005.
      Observed.p = ifelse(Observed.p == 0, 0.0005, Observed.p),
      
      Effect = as.character(Effect),
      Predictor = as.character(Predictor),
      Outcome = as.character(Outcome),
      Famly = as.character(Famly)
    ) %>%
    filter(Effect %in% c("S1", "S2"))
  
  # Pair S1 and S2 within predictor-outcome paths
  joint_df <- df_clean %>%
    group_by(Famly, Predictor, Outcome) %>%
    summarise(
      S1_p = Observed.p[Effect == "S1"][1],
      S2_p = Observed.p[Effect == "S2"][1],
      n_S1 = sum(Effect == "S1"),
      n_S2 = sum(Effect == "S2"),
      .groups = "drop"
    ) %>%
    mutate(
      complete_pair = n_S1 == 1 & n_S2 == 1 & !is.na(S1_p) & !is.na(S2_p),
      nonpair_detected = !complete_pair,
      
      S1S2_p = ifelse(
        complete_pair,
        mapply(
          function(p1, p2) cauchy_combine(c(p1, p2)),
          S1_p,
          S2_p
        ),
        NA_real_
      ),
      
      path = paste0(Outcome, " ON ", Predictor)
    ) %>%
    
    # Order p-values BEFORE BH correction for transparency
    group_by(Famly) %>%
    arrange(S1S2_p, .by_group = TRUE) %>%
    mutate(
      BH_rank = ifelse(!is.na(S1S2_p), row_number(), NA_integer_),
      BH_family_n = sum(!is.na(S1S2_p)),
      S1S2_p_BH = p.adjust(S1S2_p, method = "BH"),
      S1S2_sig_BH = S1S2_p_BH < .05
    ) %>%
    ungroup() %>%
    select(
      Famly,
      BH_rank,
      BH_family_n,
      path,
      Predictor,
      Outcome,
      S1_p,
      S2_p,
      S1S2_p,
      S1S2_p_BH,
      S1S2_sig_BH,
      n_S1,
      n_S2,
      complete_pair,
      nonpair_detected
    )
  
  # -------------------------------
  # Printed audit report
  # -------------------------------
  
  report_df <- joint_df %>%
    group_by(Famly) %>%
    summarise(
      total_predictor_outcome_sets = n(),
      complete_pairs_found = sum(complete_pair, na.rm = TRUE),
      nonpairs_detected = sum(nonpair_detected, na.rm = TRUE),
      duplicated_or_missing_S1 = sum(n_S1 != 1, na.rm = TRUE),
      duplicated_or_missing_S2 = sum(n_S2 != 1, na.rm = TRUE),
      contrasts_in_BH_family = first(BH_family_n),
      significant_after_BH = sum(S1S2_sig_BH, na.rm = TRUE),
      .groups = "drop"
    )
  
  cat("\n--------------------------------------------\n")
  cat("Joint S1/S2 p-value report\n")
  cat("--------------------------------------------\n")
  cat("Joint p-value method: Cauchy combination test\n")
  cat("Multiple-testing correction: Benjamini-Hochberg FDR\n")
  cat("BH correction applied within: Famly\n")
  cat("P-values printed as 0.000 were recoded to: 0.0005\n")
  cat("--------------------------------------------\n\n")
  
  print(report_df)
  
  cat("\nNote: 'nonpairs_detected' means the path did not contain exactly one S1 and exactly one S2 p-value.\n")
  cat("Rows with nonpairs are retained but receive NA for S1S2_p and BH-adjusted p-values.\n")
  cat("--------------------------------------------\n\n")
  
  return(joint_df)
}


# -----------------------------------------
# applied 
# -----------------------------------------

# cognitive reappraisal 
CR_joint <- make_joint_cauchy_df(CR)

CR_joint

#emotional suppression
ES_joint <- make_joint_cauchy_df(ES)

ES_joint


# -----------------------------------------
# save jointly p-value adjusted tables 
# -----------------------------------------

write.csv(CR_joint, "CR_JOINT_P_VAL_BH_ADJ.csv", row.names = F)
write.csv(ES_joint, "ES_JOINT_P_VAL_BH_ADJ.csv", row.names = F)
