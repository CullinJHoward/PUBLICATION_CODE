#Library

library(dplyr)
library(psych)
library(stringr)
library(tidyr)
library(purrr)
library(MplusAutomation)
library(tidyverse)

# Set working directory 

setwd("/home/cjh37695/ABCD_PROJECTS/HORMESIS_METHODS/")

# LOAD DF 

df <- read.csv("ABCD_HORMESIS_METHODS_REDUCED_PREP_STEP1.csv")


################################################################################
################ LOAD IN HARSH SES FACTO SCORES 

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

################################################################################
############### POMS, CENTER, AND WINSORIZE STUFF

## POMS FUNCTION

POMS <- function(data, cols,
                 min_val = NULL, max_val = NULL,
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
      
      # tolerance-based comparison
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
    x_poms <- 100 * (x_work - this_min) / (this_max - this_min)
    message("[", var, "] POMS=Y; formula=100*(x-min)/(max-min)")
    
    # Step 4: centering
    if (center) {
      if (is.null(center_at)) {
        center_value <- mean(x_poms, na.rm = TRUE)
        message("[", var, "] CENTER=Y; at mean of transformed values=", round(center_value, 4))
      } else {
        center_value <- center_at
        message("[", var, "] CENTER=Y; at user-specified POMS value=", center_value)
      }
      
      x_final <- x_poms - center_value
    } else {
      x_final <- x_poms
      message("[", var, "] CENTER=N")
    }
    
    new_name <- paste0(prefix, var)
    out_data[[new_name]] <- x_final
    message("[", var, "] OUTPUT=", new_name)
    
    # Step 5: print diagnostic plot directly
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
      
      # explicit label placement to avoid dotted-box issue
      x_rng <- range(x_final, na.rm = TRUE)
      y_max <- max(hist(x_final, plot = FALSE, breaks = bins)$counts, na.rm = TRUE)
      
      label_df <- data.frame(
        panel = factor("Transformed", levels = c("Original", "Transformed")),
        x = x_rng[2] - 0.02 * diff(x_rng),
        y = y_max * 0.95,
        lab = paste0("# winsorized = ", n_wins)
      )
      
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
          title = paste("POMS Transformation:", var),
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


FULL_df <- FULL_df %>%
  POMS(c(R_mNBsaf_1, R_mNBsaf_3, R_mNBsaf_5, R_mNBsaf_7, R_mNBsaf_9, R_mNBsaf_11), min_val = 1, max_val = 5, wins = TRUE, trim = .01, plot = TRUE) %>%
  POMS(c(R_mSchEnv_1, R_mSchEnv_3, R_mSchEnv_5, R_mSchEnv_7, R_mSchEnv_9), min_val = 1, max_val = 4, wins = TRUE, trim = .01, plot = TRUE)%>%
  POMS(c(mFconP_1, mFconP_3, mFconP_5, mFconP_7, mFconP_9, mFconP_11, mFconY_1,
         mFconY_3, mFconY_5, mFconY_7, mFconY_9, mFconY_11), min_val = 0, max_val = 1, wins = TRUE, trim = .01, plot = TRUE)%>%
  POMS(c(R_mHOMEsf_11), min_val = 0, max_val = 4.8, wins = TRUE, trim = .01, plot = TRUE)%>%
  POMS(c(fHSES), wins = TRUE, trim = .01, plot = TRUE)


################################################################################
################# PLOT BIVARIATE ASSOCIATIONS 


### PLOT FUNCTION 
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggExtra)

plot_bivariate_array <- function(data,
                                 id_var,
                                 data_format = c("wide", "long"),
                                 time_var = NULL,
                                 x_static = NULL,
                                 x_long = NULL,
                                 y_static = NULL,
                                 y_long = NULL,
                                 x_vars = NULL,
                                 y_vars = NULL,
                                 names_pattern = "^(.*)_(\\d+)$",
                                 pair_mode = c("all", "same_wave"),
                                 sample_n = 800,
                                 point_alpha = 0.35,
                                 point_size = 1.5,
                                 print_plots = FALSE,
                                 save_file = NULL,
                                 width = 11,
                                 height = 8.5,
                                 plots_per_page = 4) {
  
  data_format <- match.arg(data_format)
  pair_mode   <- match.arg(pair_mode)
  
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("dplyr required")
  if (!requireNamespace("tidyr", quietly = TRUE)) stop("tidyr required")
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 required")
  if (!requireNamespace("ggExtra", quietly = TRUE)) stop("ggExtra required")
  if (!requireNamespace("gridExtra", quietly = TRUE)) stop("gridExtra required")
  
  #-----------------------------------
  # helper functions
  #-----------------------------------
  
  check_vars_exist <- function(vars, data_names, arg_name) {
    if (is.null(vars)) return(invisible(NULL))
    missing_vars <- setdiff(vars, data_names)
    if (length(missing_vars) > 0) {
      stop(arg_name, " not found: ", paste(missing_vars, collapse = ", "))
    }
  }
  
  check_stems_exist <- function(stems, data_names, arg_name) {
    if (is.null(stems)) return(invisible(NULL))
    for (s in stems) {
      patt <- paste0("^", s, "_\\d+$")
      hits <- grep(patt, data_names, value = TRUE)
      if (length(hits) == 0) {
        stop("No columns found for ", arg_name, " stem: ", s)
      }
    }
  }
  
  get_repeated_cols <- function(stems, data_names) {
    if (is.null(stems)) return(character(0))
    patt <- paste0("^(", paste(stems, collapse = "|"), ")_\\d+$")
    grep(patt, data_names, value = TRUE)
  }
  
  #-----------------------------------
  # basic checks
  #-----------------------------------
  
  if (!id_var %in% names(data)) {
    stop("id_var not found in data")
  }
  
  if (data_format == "wide") {
    check_vars_exist(x_static, names(data), "x_static")
    check_vars_exist(y_static, names(data), "y_static")
    check_stems_exist(x_long, names(data), "x_long")
    check_stems_exist(y_long, names(data), "y_long")
  }
  
  #-----------------------------------
  # reshape data
  #-----------------------------------
  
  if (data_format == "wide") {
    
    x_rep_cols <- get_repeated_cols(x_long, names(data))
    y_rep_cols <- get_repeated_cols(y_long, names(data))
    
    static_cols   <- unique(c(x_static, y_static))
    repeated_cols <- unique(c(x_rep_cols, y_rep_cols))
    
    keep_cols <- unique(c(id_var, static_cols, repeated_cols))
    keep_cols <- keep_cols[keep_cols %in% names(data)]
    
    df_sub <- data[, keep_cols, drop = FALSE]
    
    if (length(repeated_cols) > 0) {
      df_rep <- df_sub %>%
        dplyr::select(dplyr::all_of(c(id_var, repeated_cols))) %>%
        tidyr::pivot_longer(
          cols = -dplyr::all_of(id_var),
          names_to = c(".value", "wave"),
          names_pattern = names_pattern
        ) %>%
        dplyr::mutate(wave = suppressWarnings(as.numeric(wave)))
    } else {
      df_rep <- df_sub %>%
        dplyr::distinct(dplyr::across(dplyr::all_of(id_var))) %>%
        dplyr::mutate(wave = NA)
    }
    
    if (length(static_cols) > 0) {
      df_static <- df_sub %>%
        dplyr::select(dplyr::all_of(c(id_var, static_cols))) %>%
        dplyr::distinct()
      
      df_long <- df_rep %>%
        dplyr::left_join(df_static, by = id_var)
    } else {
      df_long <- df_rep
    }
    
    x_names <- unique(c(x_static, x_long))
    y_names <- unique(c(y_static, y_long))
    
    x_df <- df_long %>%
      dplyr::select(dplyr::all_of(c(id_var, "wave", x_names))) %>%
      tidyr::pivot_longer(
        cols = dplyr::all_of(x_names),
        names_to = "x_var",
        values_to = "x_value"
      )
    
    y_df <- df_long %>%
      dplyr::select(dplyr::all_of(c(id_var, "wave", y_names))) %>%
      tidyr::pivot_longer(
        cols = dplyr::all_of(y_names),
        names_to = "y_var",
        values_to = "y_value"
      )
    
    if (pair_mode == "same_wave") {
      plot_data <- x_df %>%
        dplyr::inner_join(y_df, by = c(id_var, "wave"))
    } else {
      plot_data <- x_df %>%
        dplyr::inner_join(
          y_df,
          by = id_var,
          relationship = "many-to-many"
        ) %>%
        dplyr::mutate(wave = wave.y)
    }
  }
  
  #-----------------------------------
  # plot generator
  #-----------------------------------
  
  make_joint_plot <- function(df, x_name, y_name) {
    
    df_sub <- df %>%
      dplyr::filter(
        x_var == x_name,
        y_var == y_name,
        !is.na(x_value),
        !is.na(y_value)
      )
    
    if (nrow(df_sub) == 0) return(NULL)
    
    if (nrow(df_sub) > sample_n) {
      set.seed(123)
      df_sub <- df_sub %>% dplyr::sample_n(sample_n)
    }
    
    df_sub$wave <- as.factor(df_sub$wave)
    
    p <- ggplot2::ggplot(
      df_sub,
      ggplot2::aes(x = x_value, y = y_value, color = wave)
    ) +
      ggplot2::geom_point(alpha = point_alpha, size = point_size) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        legend.position = "bottom",
        plot.title = ggplot2::element_text(hjust = 0.5)
      ) +
      ggplot2::labs(
        title = paste(y_name, "vs", x_name),
        x = x_name,
        y = y_name,
        color = "Wave"
      )
    
    ggExtra::ggMarginal(
      p,
      type = "histogram",
      margins = "both",
      groupColour = FALSE,
      groupFill = FALSE
    )
  }
  
  combos <- expand.grid(
    x_var = x_names,
    y_var = y_names,
    stringsAsFactors = FALSE
  )
  
  plot_list <- vector("list", nrow(combos))
  
  for (i in seq_len(nrow(combos))) {
    plot_list[[i]] <- make_joint_plot(
      plot_data,
      combos$x_var[i],
      combos$y_var[i]
    )
  }
  
  plot_list <- plot_list[!sapply(plot_list, is.null)]
  
  #-----------------------------------
  # printing logic
  #-----------------------------------
  
  print_pages <- function(plots) {
    
    n <- length(plots)
    idx <- seq(1, n, by = plots_per_page)
    
    for (i in idx) {
      
      page_plots <- plots[i:min(i + plots_per_page - 1, n)]
      
      gridExtra::grid.arrange(
        grobs = page_plots,
        ncol = 2,
        nrow = 2
      )
    }
  }
  
  if (!is.null(save_file)) {
    
    grDevices::pdf(save_file, width = width, height = height)
    print_pages(plot_list)
    grDevices::dev.off()
    
  } else if (print_plots) {
    
    print_pages(plot_list)
    
  }
  
  invisible(list(
    data = plot_data,
    plots = plot_list
  ))
}

## APPLIED 

res <- plot_bivariate_array(
  data = FULL_df,
  id_var = "subID",
  data_format = "wide",
  x_static = c("WP_R_mHOMEsf_11", "WP_fHSES"),
  x_long = c("WP_R_mSchEnv", "WP_mFconP", "WP_mFconY"),
  y_long = c("mNgUrg", "mPers", "mPlan", "mPsUrg", "mSenSe",
             "mBISy", "mDRVy", "mFnSky", "mRwRspy", "mAttune",
             "mDist", "mNgSec", "mCata", "mEmoSup", "mReApp"),
  pair_mode = "all",
  sample_n = 450,
  save_file = "BIVARIATE_STRESS_DISTRIBUTIONS.pdf"
)

################################################################################
###################### CENTERING PREDICTORS AND MAKING QUADRATIC TERMS 

names(FULL_df)

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
mFconY_vars = c("WP_mFconY_1", "WP_mFconY_3", "WP_mFconY_5",
           "WP_mFconY_7", "WP_mFconY_9", "WP_mFconY_11")


FULL_df_CENT <- CenterForge(
  data = FULL_df,
  vars = mFconY_vars,
  types = c("GMC", "LMC"),
  save_means = FALSE,
  set_name = "mFconY",
  poly = "quadratic"
)


mFconP_vars = c( "WP_mFconP_1", "WP_mFconP_3", "WP_mFconP_5", 
                 "WP_mFconP_7","WP_mFconP_9", "WP_mFconP_11")


FULL_df_CENT <- CenterForge(
  data = FULL_df_CENT,
  vars = mFconP_vars,
  types = c("GMC", "LMC"),
  save_means = FALSE,
  set_name = "mFconP",
  poly = "quadratic"
)

HSES = c("WP_fHSES")

FULL_df_CENT <- CenterForge(
  data = FULL_df_CENT,
  vars = HSES,
  types = c("GMC"),
  save_means = FALSE,
  set_name = "fHSES",
  poly = "quadratic"
)


NBH_SAF = c("WP_R_mNBsaf_1", "WP_R_mNBsaf_3", "WP_R_mNBsaf_5", 
                "WP_R_mNBsaf_7", "WP_R_mNBsaf_9", "WP_R_mNBsaf_11")


FULL_df_CENT <- CenterForge(
  data = FULL_df_CENT,
  vars = NBH_SAF,
  types = c("GMC", "LMC"),
  save_means = FALSE,
  set_name = "NBH_SAF",
  poly = "quadratic"
)

SCHENV = c("WP_R_mSchEnv_1", "WP_R_mSchEnv_3", "WP_R_mSchEnv_5", 
           "WP_R_mSchEnv_7", "WP_R_mSchEnv_9")


FULL_df_CENT <- CenterForge(
  data = FULL_df_CENT,
  vars = SCHENV,
  types = c("GMC", "LMC"),
  save_means = FALSE,
  set_name = "SCHENV",
  poly = "quadratic"
)


################################################################################
################# SAVE DATA 

## FIX SEX TO EFFECTS CODING RATHER THAT CATEGOCIAL 
FULL_df_CENT$Y_SEX_EFF <- ifelse(FULL_df_CENT$Y_SEX == 2,1,
                                 ifelse(FULL_df_CENT$Y_SEX == 1,-1,NA))

table(FULL_df_CENT$Y_SEX)
  
  
## REDUCE FOR SEM 
names(FULL_df_CENT)

FULL_df_CENT_RED <- FULL_df_CENT %>% 
  select(-c("mHOMEsf_11",         
            "mFconP_1", "mFconP_3", "mFconP_5", "mFconP_7", "mFconP_9",  "mFconP_11",          
            "mFconY_1", "mFconY_3", "mFconY_5", "mFconY_7",  "mFconY_9", "mFconY_11",          
            "mNBsaf_1", "mNBsaf_3", "mNBsaf_5", "mNBsaf_7", "mNBsaf_9", "mNBsaf_11",          
            "mPedsup_7", "mPsupr_7", "mSchEnv_1", "mSchEnv_3",  "mSchEnv_5",  "mSchEnv_7",          
            "mSchEnv_9", "R_mNBsaf_1",         
            "R_mNBsaf_3", "R_mNBsaf_5", "R_mNBsaf_7", "R_mNBsaf_9", "R_mNBsaf_11", "R_mSchEnv_1",        
            "R_mSchEnv_3", "R_mSchEnv_5", "R_mSchEnv_7", "R_mSchEnv_9",
            "WP_R_mNBsaf_1", "WP_R_mNBsaf_3", "WP_R_mNBsaf_5", "WP_R_mNBsaf_7", "WP_R_mNBsaf_9", "WP_R_mNBsaf_11",     
            "WP_R_mSchEnv_1", "WP_R_mSchEnv_3", "WP_R_mSchEnv_5", "WP_R_mSchEnv_7",  "WP_R_mSchEnv_9", "WP_mFconP_1",        
            "WP_mFconP_3", "WP_mFconP_5", "WP_mFconP_7", "WP_mFconP_9",  "WP_mFconP_11",  "WP_mFconY_1",        
            "WP_mFconY_3", "WP_mFconY_5",  "WP_mFconY_7", "WP_mFconY_9",  "WP_mFconY_11",  "WP_R_mHOMEsf_11",    
            "WP_fHSES"))
#FULL
write.csv(FULL_df_CENT,"ABCD_HORMESIS_METHODS_FULL_STRUC_STEP2_V2.csv", row.names = F)
prepareMplusData(FULL_df_CENT,"ABCD_HORMESIS_METHODS_FULL_STRUC_STEP2_V2.dat", inpfile =T)

## REDUCED SEM 

write.csv(FULL_df_CENT_RED,"ABCD_HORMESIS_METHODS_REDUCED_STRUC_STEP2_V2.csv", row.names = F)
prepareMplusData(FULL_df_CENT_RED,"ABCD_HORMESIS_METHODS_REDUCED_STRUC_STEP2_V2.dat", inpfile =T)


################################################################################
################# ARCHIVE 


###### PIVOT TO LONG FORMAT 
names(FULL_df_CENT)

# ID INVARIATE VARIABLES 
id_vars <- c(
  "subID",
  "SiteID",
  "FamilyID",
  "Y_HISP",
  "Y_SEX",
  "ppensity",
  "NUMID",          
  "HiParEdu_1",
  "R_HiParEdu_1",
  "C_HiParEdu_1",
  "INCOME6L_1",
  "C_INCOME6L_1",
  "MarWrkSt_1",
  "PRBlPov_1",
  "PRCrowd_1",      
  "PRLowEd_1",
  "PRUnEmp_1",
  "R_SocMob_1",
  "mPedsup_7",      
  "mPsupr_7",
  "mHOMEsf_11",
  "R_mHOMEsf_11",
  "WP_R_mHOMEsf_11",
  "fHSES",
  "WP_fHSES"
)

# ID VARYING VARIABLES (SUFFIX, NOT FULL NAME)

repeated_stems <- c(
  "Y_AGE",
  "PMC_Y_AGE",
  "mFconP",
  "WP_mFconP",
  "mFconY",
  "WP_mFconY",
  "mNBsaf",
  "R_mNBsaf",
  "mSchEnv",
  "R_mSchEnv",
  "WP_R_mSchEnv",
  "mNgUrg",
  "mPers",
  "mPlan",
  "mPsUrg",
  "mSenSe",
  "mBISy",
  "mDRVy",
  "mFnSky",
  "mRwRspy",
  "mAttune",
  "mDist",
  "mNgSec",
  "mCata",
  "mEmoSup",
  "mReApp"
)

# KEEP ONLY WHAT IS WANTED 
pattern <- paste0("^(", paste(repeated_stems, collapse = "|"), ")_\\d+$")

df_sub <- FULL_df %>%
  select(all_of(id_vars), matches(pattern))

# PIVOT TO LONG 
df_long <- df_sub %>%
  pivot_longer(
    cols = matches(pattern),
    names_to = c(".value", "wave"),
    names_pattern = "^(.*)_(\\d+)$"
  ) %>%
  mutate(
    wave = as.numeric(wave)
  ) %>%
  arrange(subID, wave)



