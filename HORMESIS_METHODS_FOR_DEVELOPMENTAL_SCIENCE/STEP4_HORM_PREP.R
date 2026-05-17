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

df <- read.csv("ABCD_HORM_METH_SEM_MEASUREMEMT_4.24.26.csv")


################################################################################
##################### MERGE IN HSES VARIABLE


## LOAD IN SAVED FACTORS
HSES <- read.table("HARSH_SES_FACTOR.txt", header = FALSE, na.strings = "*")

##Subset the data frame to keep only NUMID, and the latent scores you want

HSES_FACTORS <- HSES[, c("V7", "V10")]

# Give names to the variables
colnames(HSES_FACTORS) <- c("fHSES", "NUMID")

## WINSORIZE THE OUTLIERS 
HSES_FACTORS$fHSES <- winsor(HSES_FACTORS$fHSES, trim = 0.01)


## MERGE IN 
FULL_df <- full_join(df, HSES_FACTORS, by = "NUMID") 

## POMS SCORE HSES factor 

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


## applied 

FULL_df_POMS <- FULL_df %>%
  POMS(c(fHSES), wins = TRUE, trim = .01,  scale_to = 10, plot = TRUE)
  

################################################################################
#####################  CENTER AND POLYNOMIAL VARIABLS 

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


## APPLIED 
names(FULL_df_POMS)

age = c("Y_AGE_1", "Y_AGE_3", "Y_AGE_5", "Y_AGE_7", "Y_AGE_9",  "Y_AGE_11" ) # pooled longitudinal grand mean centered
#inv_adv = c("WP_fHSES", "WP_R_mHOMEsf_11" ) # grand mean centered 
#NBH_SAF = c("WP_R_mNBsaf_1", "WP_R_mNBsaf_3", "WP_R_mNBsaf_5", "WP_R_mNBsaf_7", "WP_R_mNBsaf_9", "WP_R_mNBsaf_11") # pooled longitudinal grand mean centered
#SCH_ENV = c("WP_R_mSchEnv_1", "WP_R_mSchEnv_3", "WP_R_mSchEnv_5", "WP_R_mSchEnv_7", "WP_R_mSchEnv_9", "WP_mFconP_1") # pooled longitudinal grand mean centered
#FAM_CONp = c("WP_mFconP_1", "WP_mFconP_3", "WP_mFconP_5", "WP_mFconP_7", "WP_mFconP_9", "WP_mFconP_11") # pooled longitudinal grand mean centered
#FAM_CONy = c("WP_mFconY_1", "WP_mFconY_3", "WP_mFconY_5", "WP_mFconY_7", "WP_mFconY_9", "WP_mFconY_11") # pooled longitudinal grand mean centered
mot = c("WP_dtiMnMOT_1", "WP_dtiMnMOT_5", "WP_dtiMnMOT_9") # grand mean centered 
GLB_WHT = c("WP_dts_FA_ALL_1", "WP_dts_FA_ALL_5", "WP_dts_FA_ALL_9") # pooled longitudinal grand mean centered
demo = c("INCOME6L", "HiParEdu_1") # grand mean centered 
pub = c("PUBlev_1", "PUBlev_3", "PUBlev_5", "PUBlev_7", "PUBlev_9", "PUBlev_11") # pooled longitudinal grand mean centered
#slp_exsmn = c("ExSmnS_1", "ExSmnS_3", "ExSmnS_5", "ExSmnS_7", "ExSmnS_9", "ExSmnS_11") # pooled longitudinal grand mean centered
#slp_Inmn = c("InMnSlpS_1", "InMnSlpS_3", "InMnSlpS_5", "InMnSlpS_7", "InMnSlpS_9", "InMnSlpS_11") # pooled longitudinal grand mean centered
#slp_Dis = c("SlpDisS_1", "SlpDisS_3", "SlpDisS_5", "SlpDisS_7", "SlpDisS_9", "SlpDisS_11") # pooled longitudinal grand mean centered
BP_Dys = c("DysBpM_5", "DysBpM_7", "DysBpM_9", "DysBpM_11") # pooled longitudinal grand mean centered
BP_Sys = c("SysBpM_5", "SysBpM_7", "SysBpM_9", "SysBpM_11") # pooled longitudinal grand mean centered
Pmon = c("WP_mPmon_1", "WP_mPmon_3", "WP_mPmon_5", "WP_mPmon_7", "WP_mPmon_9", "WP_mPmon_11") # pooled longitudinal grand mean centered
SupFrd = c("WP_mPNH_5", "WP_mPNH_7", "WP_mPNH_9") # pooled longitudinal grand mean centered
#PhsAct = c("WP_PHYACTy_1", "WP_PHYACTy_5", "WP_PHYACTy_7", "WP_PHYACTy_9") # pooled longitudinal grand mean centered

FULL_df_POMS_CENT <- CenterForge(
  data = FULL_df_POMS,
  vars = age,
  types = c("LMC"),
  save_means = FALSE,
  set_name = "age",
  poly = "none"
)

# FULL_df_POMS_CENT <- CenterForge(
#   data = FULL_df_POMS_CENT,
#   vars = inv_adv,
#   types = c("GMC"),
#   save_means = FALSE,
#   set_name = "inv_adv",
#   poly = "quadratic"
# )

# FULL_df_POMS_CENT <- CenterForge(
#   data = FULL_df_POMS_CENT,
#   vars = NBH_SAF,
#   types = c("LMC"),
#   save_means = FALSE,
#   set_name = "NBH_SAF",
#   poly = "quadratic"
# )

# FULL_df_POMS_CENT <- CenterForge(
#   data = FULL_df_POMS_CENT,
#   vars = SCH_ENV,
#   types = c("LMC"),
#   save_means = FALSE,
#   set_name = "SCH_ENV",
#   poly = "quadratic"
# )

# FULL_df_POMS_CENT <- CenterForge(
#   data = FULL_df_POMS_CENT,
#   vars = SCH_ENV,
#   types = c("LMC"),
#   save_means = FALSE,
#   set_name = "SCH_ENV",
#   poly = "quadratic"
# )

# FULL_df_POMS_CENT <- CenterForge(
#   data = FULL_df_POMS_CENT,
#   vars = FAM_CONy,
#   types = c("LMC"),
#   save_means = FALSE,
#   set_name = "FAM_CONy",
#   poly = "quadratic"
# )

FULL_df_POMS_CENT <- CenterForge(
  data = FULL_df_POMS_CENT,
  vars = mot,
  types = c("GMC"),
  save_means = FALSE,
  set_name = "mot",
  poly = "none"
)

FULL_df_POMS_CENT <- CenterForge(
  data = FULL_df_POMS_CENT,
  vars = GLB_WHT,
  types = c("LMC"),
  save_means = FALSE,
  set_name = "GLB_WHT",
  poly = "none"
)

FULL_df_POMS_CENT <- CenterForge(
  data = FULL_df_POMS_CENT,
  vars = demo,
  types = c("GMC"),
  save_means = FALSE,
  set_name = "demo",
  poly = "none"
)

FULL_df_POMS_CENT <- CenterForge(
  data = FULL_df_POMS_CENT,
  vars = pub,
  types = c("LMC"),
  save_means = FALSE,
  set_name = "pub",
  poly = "none"
)

# FULL_df_POMS_CENT <- CenterForge(
#   data = FULL_df_POMS_CENT,
#   vars = slp_exsmn,
#   types = c("LMC"),
#   save_means = FALSE,
#   set_name = "slp_exsmn",
#   poly = "none"
# )

# FULL_df_POMS_CENT <- CenterForge(
#   data = FULL_df_POMS_CENT,
#   vars = slp_Inmn,
#   types = c("LMC"),
#   save_means = FALSE,
#   set_name = "slp_Inmn",
#   poly = "none"
# )

# FULL_df_POMS_CENT <- CenterForge(
#   data = FULL_df_POMS_CENT,
#   vars = slp_Dis,
#   types = c("LMC"),
#   save_means = FALSE,
#   set_name = "slp_Dis",
#   poly = "none"
# )

FULL_df_POMS_CENT <- CenterForge(
  data = FULL_df_POMS_CENT,
  vars = BP_Dys,
  types = c("LMC"),
  save_means = FALSE,
  set_name = "BP_Dys",
  poly = "none"
)

FULL_df_POMS_CENT <- CenterForge(
  data = FULL_df_POMS_CENT,
  vars = BP_Sys,
  types = c("LMC"),
  save_means = FALSE,
  set_name = "BP_Sys",
  poly = "none"
)

FULL_df_POMS_CENT <- CenterForge(
  data = FULL_df_POMS_CENT,
  vars = Pmon,
  types = c("LMC"),
  save_means = FALSE,
  set_name = "Pmon",
  poly = "none"
)

FULL_df_POMS_CENT <- CenterForge(
  data = FULL_df_POMS_CENT,
  vars = SupFrd,
  types = c("LMC"),
  save_means = FALSE,
  set_name = "SupFrd",
  poly = "none"
)

# FULL_df_POMS_CENT <- CenterForge(
#   data = FULL_df_POMS_CENT,
#   vars = PhsAct,
#   types = c("LMC"),
#   save_means = FALSE,
#   set_name = "PhsAct",
#   poly = "none"
# )

###############################################################################
#################### REDUCE VARIABLE SET
names(FULL_df_POMS_CENT)

FULL_df_POMS_CENT_RED  <- FULL_df_POMS_CENT %>%
  select(-c("DysBpM_5",             "DysBpM_7",             "DysBpM_9",             "DysBpM_11",      
            "SysBpM_5",             "SysBpM_7",             "SysBpM_9",             "SysBpM_11",            
            "PUBlev_1",             "PUBlev_3",             "PUBlev_5",             "PUBlev_7",             "PUBlev_9",            
            "PUBlev_11",            
            "WP_dtiMnMOT_1",        "WP_dtiMnMOT_5",       
            "WP_dtiMnMOT_9",        
            "WP_mPmon_1",           "WP_mPmon_3",           "WP_mPmon_5",           "WP_mPmon_7",           "WP_mPmon_9",          
            "WP_mPmon_11",          "WP_mPNH_5",            "WP_mPNH_7",            "WP_mPNH_9"
           ))

names(FULL_df_POMS_CENT_RED)
###############################################################################
#################### SAVE DF

write.csv(FULL_df_POMS_CENT_RED,"ABCD_HORM_METH_SEM_STRUCTRUAL_4.24.26.csv", row.names = F)
#prepareMplusData(FULL_df_POMS_CENT_RED,"ABCD_HORM_METH_SEM_STRUCTRUAL_4.23.26.dat", inpfile =T)


## getting the updated 4-item parental monitoring 

PMON <- FULL_df_POMS_CENT_RED %>%
  select("subID", "LMC_WP_mPmon_1", "LMC_WP_mPmon_3" , "LMC_WP_mPmon_5", "LMC_WP_mPmon_7",     
         "LMC_WP_mPmon_9", "LMC_WP_mPmon_11")

write.csv(PMON, "PARENTAL_MONITORING_4_ITEM_SUB_4.27.26.csv", row.names = F)
