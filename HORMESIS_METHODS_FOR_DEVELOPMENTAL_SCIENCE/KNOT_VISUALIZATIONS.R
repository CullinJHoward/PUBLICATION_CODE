## library

library(dplyr)


## set up environment

setwd("C:\\Users\\cjh37695\\Dropbox\\HORMESIS_METHODS\\ANALYSIS\\")

inf_rec <- read.csv("HORMETIC_INFLECTIONS_BS100_WAVE_GMC_MOD.csv")

df <- read.csv("ABCD_HORM_METH_SEM_STRUCTRUAL_4.27.26.csv")

##-----------------------------------------------------------------------------
## helper functions  
##-----------------------------------------------------------------------------

############ transform data from long to wide for each pairing 

make_xy_long <- function(data,
                         id_var,
                         x_stems,
                         y_stems,
                         pair_mode = c("same_wave", "time_invariant_x"),
                         names_pattern = "^(.*)_(\\d+)$",
                         wave_as_numeric = TRUE,
                         drop_missing = TRUE,
                         print_diagnostics = TRUE) {
  
  require(dplyr)
  require(tidyr)
  require(stringr)
  
  pair_mode <- match.arg(pair_mode)
  
  if (!id_var %in% names(data)) {
    stop("id_var not found in data: ", id_var)
  }
  
  # Escape regex-sensitive characters in stems
  escape_regex <- function(x) {
    stringr::str_replace_all(x, "([\\^\\$\\.\\|\\?\\*\\+\\(\\)\\[\\]\\{\\}])", "\\\\\\1")
  }
  
  # Find time-varying columns from stems, e.g., WP_mNBdang -> WP_mNBdang_1, WP_mNBdang_3
  get_repeated_cols <- function(stems, data_names) {
    stem_pattern <- paste(escape_regex(stems), collapse = "|")
    pattern <- paste0("^(", stem_pattern, ")_\\d+$")
    grep(pattern, data_names, value = TRUE)
  }
  
  # Find single columns from exact names
  get_exact_cols <- function(stems, data_names) {
    stems[stems %in% data_names]
  }
  
  # Summarize repeated stem matches
  summarize_stem_matches <- function(stems, cols, label) {
    
    out <- lapply(stems, function(stem) {
      
      escaped_stem <- escape_regex(stem)
      
      matched_cols <- grep(
        paste0("^", escaped_stem, "_\\d+$"),
        cols,
        value = TRUE
      )
      
      waves <- stringr::str_match(
        matched_cols,
        paste0("^", escaped_stem, "_(\\d+)$")
      )[, 2]
      
      if (wave_as_numeric) {
        waves <- suppressWarnings(as.numeric(waves))
      }
      
      tibble(
        set = label,
        stem = stem,
        n_columns = length(matched_cols),
        waves = ifelse(
          length(waves) == 0,
          NA_character_,
          paste(sort(waves), collapse = ", ")
        ),
        columns = ifelse(
          length(matched_cols) == 0,
          NA_character_,
          paste(matched_cols, collapse = ", ")
        )
      )
    })
    
    bind_rows(out)
  }
  
  # -----------------------------
  # Y is always time-varying here
  # -----------------------------
  
  y_cols <- get_repeated_cols(y_stems, names(data))
  y_summary <- summarize_stem_matches(y_stems, y_cols, "Y")
  
  if (length(y_cols) == 0) {
    stop("No Y columns found for requested y_stems.")
  }
  
  y_long <- data %>%
    select(all_of(c(id_var, y_cols))) %>%
    pivot_longer(
      cols = -all_of(id_var),
      names_to = c("y_var", "wave"),
      names_pattern = names_pattern,
      values_to = "y_value"
    )
  
  if (isTRUE(wave_as_numeric)) {
    y_long <- y_long %>%
      mutate(wave = suppressWarnings(as.numeric(wave)))
  }
  
  if (isTRUE(drop_missing)) {
    y_long <- y_long %>%
      filter(!is.na(y_value))
  }
  
  # -----------------------------
  # Case 1: time-varying X
  # -----------------------------
  
  if (pair_mode == "same_wave") {
    
    x_cols <- get_repeated_cols(x_stems, names(data))
    x_summary <- summarize_stem_matches(x_stems, x_cols, "X")
    
    if (length(x_cols) == 0) {
      stop("No time-varying X columns found for requested x_stems.")
    }
    
    x_long <- data %>%
      select(all_of(c(id_var, x_cols))) %>%
      pivot_longer(
        cols = -all_of(id_var),
        names_to = c("x_var", "wave"),
        names_pattern = names_pattern,
        values_to = "x_value"
      )
    
    if (isTRUE(wave_as_numeric)) {
      x_long <- x_long %>%
        mutate(wave = suppressWarnings(as.numeric(wave)))
    }
    
    if (isTRUE(drop_missing)) {
      x_long <- x_long %>%
        filter(!is.na(x_value))
    }
    
    wave_pairing_summary <- x_long %>%
      distinct(x_var, wave) %>%
      inner_join(
        y_long %>% distinct(y_var, wave),
        by = "wave"
      ) %>%
      group_by(x_var, y_var) %>%
      summarise(
        paired_waves = paste(sort(unique(wave)), collapse = ", "),
        n_paired_waves = n_distinct(wave),
        .groups = "drop"
      )
    
    paired_long <- x_long %>%
      inner_join(y_long, by = c(id_var, "wave")) %>%
      arrange(.data[[id_var]], wave, x_var, y_var)
    
    if (isTRUE(print_diagnostics)) {
      cat("\n-----------------------------\n")
      cat("PAIR MODE: same_wave\n")
      cat("Joining X and Y by:", id_var, "+ wave\n")
      
      cat("\nX variable matches:\n")
      print(x_summary, n = Inf)
      
      cat("\nY variable matches:\n")
      print(y_summary, n = Inf)
      
      cat("\nSame-wave pairing summary:\n")
      print(wave_pairing_summary, n = Inf)
      
      cat("\nFinal paired long dataframe:\n")
      cat("Rows:", nrow(paired_long), "\n")
      cat("Unique IDs:", dplyr::n_distinct(paired_long[[id_var]]), "\n")
      cat("Observed paired waves:", paste(sort(unique(paired_long$wave)), collapse = ", "), "\n")
      cat("-----------------------------\n\n")
    }
  }
  
  # -----------------------------
  # Case 2: time-invariant X
  # -----------------------------
  
  if (pair_mode == "time_invariant_x") {
    
    x_cols <- get_exact_cols(x_stems, names(data))
    
    if (length(x_cols) == 0) {
      stop(
        "No time-invariant X columns found. For pair_mode = 'time_invariant_x', ",
        "x_stems should be exact column names, e.g., x_stems = c('fHSES', 'WP_fHSES')."
      )
    }
    
    missing_x <- setdiff(x_stems, x_cols)
    
    if (length(missing_x) > 0) {
      warning(
        "These requested time-invariant X columns were not found: ",
        paste(missing_x, collapse = ", ")
      )
    }
    
    x_long <- data %>%
      select(all_of(c(id_var, x_cols))) %>%
      pivot_longer(
        cols = -all_of(id_var),
        names_to = "x_var",
        values_to = "x_value"
      )
    
    if (isTRUE(drop_missing)) {
      x_long <- x_long %>%
        filter(!is.na(x_value))
    }
    
    paired_long <- x_long %>%
      inner_join(y_long, by = id_var) %>%
      arrange(.data[[id_var]], x_var, y_var, wave)
    
    time_invariant_pairing_summary <- paired_long %>%
      group_by(x_var, y_var) %>%
      summarise(
        y_waves_stacked = paste(sort(unique(wave)), collapse = ", "),
        n_y_waves = n_distinct(wave),
        n_rows = n(),
        .groups = "drop"
      )
    
    if (isTRUE(print_diagnostics)) {
      cat("\n-----------------------------\n")
      cat("PAIR MODE: time_invariant_x\n")
      cat("Joining X and Y by:", id_var, "only\n")
      cat("Interpretation: each time-invariant X value is repeated across all available Y waves within the same person.\n")
      
      cat("\nX variables found:\n")
      print(x_cols)
      
      cat("\nY variable matches:\n")
      print(y_summary, n = Inf)
      
      cat("\nTime-invariant X by time-varying Y pairing summary:\n")
      print(time_invariant_pairing_summary, n = Inf)
      
      cat("\nFinal paired long dataframe:\n")
      cat("Rows:", nrow(paired_long), "\n")
      cat("Unique IDs:", dplyr::n_distinct(paired_long[[id_var]]), "\n")
      cat("Observed Y waves:", paste(sort(unique(paired_long$wave)), collapse = ", "), "\n")
      cat("-----------------------------\n\n")
    }
  }
  
  paired_long
}

############ center y-values within waves
center_y_by_wave <- function(dat) {
  dat %>%
    dplyr::group_by(wave) %>%
    dplyr::mutate(
      y_value = y_value - mean(y_value, na.rm = TRUE)
    ) %>%
    dplyr::ungroup()
}

############ extract plotting information 

extract_spline_info <- function(inf_df,
                                x_target,
                                y_target,
                                print_info = TRUE) {
  
  row <- inf_df %>%
    dplyr::filter(
      x_variable == x_target,
      y_variable == y_target
    )
  
  if (nrow(row) == 0) {
    stop("No matching X-Y association found.")
  }
  
  if (nrow(row) > 1) {
    stop("More than one matching X-Y association found. Check x_variable/y_variable names.")
  }
  
  needed_cols <- c(
    "x_variable", "y_variable",
    "horm_vrt",
    "horm_ref", "horm_ref_x", "horm_inf_x",
    "slope1_b", "slope2_b"
  )
  
  missing_cols <- setdiff(needed_cols, names(row))
  
  if (length(missing_cols) > 0) {
    stop(
      "The following needed columns are missing from inf_df: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  get_col <- function(col, default = NA_real_) {
    if (!col %in% names(row)) return(default)
    val <- row[[col]][1]
    if (length(val) == 0) return(default)
    val
  }
  
  knot_x <- row$horm_vrt
  b1     <- row$slope1_b
  b2     <- row$slope2_b
  
  ref_y  <- row$horm_ref
  ref_x  <- row$horm_ref_x
  inf_x  <- row$horm_inf_x
  
  # Backward-compatible fallback for older ScatterSplatter CSVs.
  # New CSVs should contain knot_y_left and knot_y_right directly.
  fallback_knot_y <- ref_y + b1 * (knot_x - ref_x)
  
  knot_y_left  <- get_col("knot_y_left", fallback_knot_y)
  knot_y_right <- get_col("knot_y_right", knot_y_left)
  knot_y       <- get_col("knot_y", knot_y_left)
  knot_y_gap   <- get_col("knot_y_gap", knot_y_right - knot_y_left)
  
  out <- list(
    x_variable = row$x_variable,
    y_variable = row$y_variable,
    
    # Core plotting parameters
    knot_x = knot_x,
    knot_y = knot_y,
    knot_y_left = knot_y_left,
    knot_y_right = knot_y_right,
    knot_y_gap = knot_y_gap,
    b1 = b1,
    b2 = b2,
    
    # Landmark values
    ref_x = ref_x,
    ref_y = ref_y,
    inf_x = inf_x,
    
    # Optional extra values if present
    slope1_p = if ("slope1_p" %in% names(row)) row$slope1_p else NA_real_,
    slope2_p = if ("slope2_p" %in% names(row)) row$slope2_p else NA_real_,
    slope1_se = if ("slope1_se" %in% names(row)) row$slope1_se else NA_real_,
    slope2_se = if ("slope2_se" %in% names(row)) row$slope2_se else NA_real_,
    glm1_intercept = get_col("glm1_intercept", NA_real_),
    glm2_intercept = get_col("glm2_intercept", NA_real_),
    glm1_high_jump = get_col("glm1_high_jump", NA_real_),
    glm2_high_jump = get_col("glm2_high_jump", NA_real_)
  )
  
  if (isTRUE(print_info)) {
    cat("\n-----------------------------\n")
    cat("Spline information extracted\n")
    cat("-----------------------------\n")
    cat("X variable:", out$x_variable, "\n")
    cat("Y variable:", out$y_variable, "\n\n")
    
    cat("Reference x:", round(out$ref_x, 4), "\n")
    cat("Reference y:", round(out$ref_y, 4), "\n")
    cat("Knot x:", round(out$knot_x, 4), "\n")
    cat("Knot y, left:", round(out$knot_y_left, 4), "\n")
    cat("Knot y, right:", round(out$knot_y_right, 4), "\n")
    cat("Knot y gap, right - left:", round(out$knot_y_gap, 6), "\n")
    cat("Inflection/toxic threshold x:", round(out$inf_x, 4), "\n\n")
    
    cat("Slope 1 b:", round(out$b1, 4), "\n")
    cat("Slope 2 b:", round(out$b2, 4), "\n")
    
    if (!is.na(out$slope1_p)) cat("Slope 1 p:", round(out$slope1_p, 4), "\n")
    if (!is.na(out$slope2_p)) cat("Slope 2 p:", round(out$slope2_p, 4), "\n")
    
    cat("-----------------------------\n\n")
  }
  
  out
}

##-----------------------------------------------------------------------------
## subset the wanted data & within-wave center the outcome 
##-----------------------------------------------------------------------------

names(df)

HSES_SUPP <- make_xy_long(
  data = df,
  id_var = "subID",
  x_stems = c("WP_fHSES"),
  y_stems = c("mEmoSup"),
  pair_mode = "time_invariant_x"
) %>%
  center_y_by_wave()

NBDANG_SUPP <- make_xy_long(
  data = df,
  id_var = "subID",
  x_stems = c("WP_mNBdang"),
  y_stems = c("mEmoSup"),
  pair_mode = "same_wave"
) %>%
  center_y_by_wave()

SCHADV_SUPP <- make_xy_long(
  data = df,
  id_var = "subID",
  x_stems = c("WP_mSchAdv"),
  y_stems = c("mEmoSup"),
  pair_mode = "same_wave"
) %>%
  center_y_by_wave()

FCON_SUPP <- make_xy_long(
  data = df,
  id_var = "subID",
  x_stems = c("WP_FAMCONY"),
  y_stems = c("mEmoSup"),
  pair_mode = "same_wave"
) %>%
  center_y_by_wave()

HSES_REAP <- make_xy_long(
  data = df,
  id_var = "subID",
  x_stems = c("WP_fHSES"),
  y_stems = c("mReApp"),
  pair_mode = "time_invariant_x"
) %>%
  center_y_by_wave()

NBDANG_REAP <- make_xy_long(
  data = df,
  id_var = "subID",
  x_stems = c("WP_mNBdang"),
  y_stems = c("mReApp"),
  pair_mode = "same_wave"
) %>%
  center_y_by_wave()

SCHADV_REAP <- make_xy_long(
  data = df,
  id_var = "subID",
  x_stems = c("WP_mSchAdv"),
  y_stems = c("mReApp"),
  pair_mode = "same_wave"
) %>%
  center_y_by_wave()

FCON_REAP <- make_xy_long(
  data = df,
  id_var = "subID",
  x_stems = c("WP_FAMCONY"),
  y_stems = c("mReApp"),
  pair_mode = "same_wave"
) %>%
  center_y_by_wave()

##-----------------------------------------------------------------------------
## extract two-lines parameter estimates for plotting  
##----------------------------------------------------------------------------- 

names(inf_rec)
table(inf_rec$x_variable)
table(inf_rec$y_variable)

HSES_SUPP_par <- extract_spline_info(
  inf_df = inf_rec,
  x_target = "WP_fHSES",
  y_target = "mEmoSup"
)

FCON_SUPP_par <- extract_spline_info(
  inf_df = inf_rec,
  x_target = "WP_FAMCONY",
  y_target = "mEmoSup"
)

NBDANG_SUPP_par <- extract_spline_info(
  inf_df = inf_rec,
  x_target = "WP_mNBdang",
  y_target = "mEmoSup"
)
SCHADV_SUPP_par <- extract_spline_info(
  inf_df = inf_rec,
  x_target = "WP_mSchAdv",
  y_target = "mEmoSup"
)

HSES_REAP_par <- extract_spline_info(
  inf_df = inf_rec,
  x_target = "WP_fHSES",
  y_target = "mReApp"
)

FCON_REAP_par <- extract_spline_info(
  inf_df = inf_rec,
  x_target = "WP_FAMCONY",
  y_target = "mReApp"
)

NBDANG_REAP_par <- extract_spline_info(
  inf_df = inf_rec,
  x_target = "WP_mNBdang",
  y_target = "mReApp"
)
SCHADV_REAP_par <- extract_spline_info(
  inf_df = inf_rec,
  x_target = "WP_mSchAdv",
  y_target = "mReApp"
)

##-----------------------------------------------------------------------------
## plotting 
##-----------------------------------------------------------------------------

make_spline_panel_grid <- function(panel_specs,
                                   domain_order = c("HSES", "Neighborhood Danger", "School Adversity", "Family Conflict"),
                                   outcome_order = c("Emotional Suppression", "Reappraisal"),
                                   
                                   # x scaling
                                   x_multiplier = 10,
                                   x_range = NULL,
                                   x_ranges = NULL,
                                   
                                   # y scaling
                                   y_ranges = NULL,
                                   y_pad_lower = 0.08,
                                   y_pad_upper = 0.18,
                                   
                                   # optional polynomial overlays
                                   polynomial = c("none", "quadratic", "cubic", "both"),
                                   quadratic_color = "#A78BFA",
                                   cubic_color = "#E69F00",
                                   polynomial_linewidth = 1.00,
                                   polynomial_linetype = "solid",
                                   
                                   # observed points
                                   show_points = FALSE,
                                   point_mode = c("gray", "colored"),
                                   n_points = NULL,
                                   point_alpha = 0.25,
                                   point_size = 0.75,
                                   seed = 1234,
                                   gray_point_color = "gray55",
                                   left_point_color = "#0072B2",
                                   right_point_color = "#C65A5A",
                                   
                                   # spline
                                   slope1_color = "#0072B2",
                                   slope2_color = "#C65A5A",
                                   slope_linewidth = 1.2,
                                   knot_line_color = "black",
                                   knot_linewidth = 0.8,
                                   
                                   # knot label
                                   show_knot_label = TRUE,
                                   knot_label_digits = 1,
                                   knot_label_size = 3.4,
                                   knot_label_top_nudge = 0.04,
                                   knot_label_fill = "white",
                                   
                                   # bottom histogram
                                   show_x_hist = TRUE,
                                   hist_bins = 35,
                                   hist_fill = "gray70",
                                   hist_color = "white",
                                   hist_alpha = 0.90,
                                   hist_height = 0.16,
                                   
                                   # legend
                                   show_legend = TRUE,
                                   legend_position = "top",
                                   legend_text_size = 10.5,
                                   legend_linewidth = 1.2,
                                   legend_height = 0.10,
                                   
                                   # appearance
                                   base_size = 12,
                                   base_family = "Arial",
                                   axis_title_size = 12,
                                   axis_text_size = 10,
                                   strip_title_size = 13,
                                   show_y_axis_text_all = FALSE,
                                   plot_margin = ggplot2::margin(5, 6, 4, 6),
                                   
                                   print_diagnostics = TRUE) {
  
  require(dplyr)
  require(ggplot2)
  require(cowplot)
  require(tibble)
  
  point_mode <- match.arg(point_mode)
  polynomial <- match.arg(polynomial)
  
  plot_quadratic <- polynomial %in% c("quadratic", "both")
  plot_cubic     <- polynomial %in% c("cubic", "both")
  
  needed_spec_cols <- c(
    "domain", "outcome",
    "domain_label", "outcome_label",
    "plot_df", "spline_par"
  )
  
  missing_spec_cols <- setdiff(needed_spec_cols, names(panel_specs))
  
  if (length(missing_spec_cols) > 0) {
    stop("panel_specs is missing: ", paste(missing_spec_cols, collapse = ", "))
  }
  
  if (!is.null(x_range) && (!is.numeric(x_range) || length(x_range) != 2)) {
    stop("x_range must be NULL or a numeric vector of length 2.")
  }
  
  if (!is.null(x_ranges) && !is.list(x_ranges)) {
    stop("x_ranges must be NULL or a named list, e.g., list('HSES' = c(0, 80)).")
  }
  
  get_named_range <- function(range_list, key1, key2 = NULL) {
    
    if (is.null(range_list)) {
      return(NULL)
    }
    
    if (!is.null(key1) && key1 %in% names(range_list)) {
      return(range_list[[key1]])
    }
    
    if (!is.null(key2) && key2 %in% names(range_list)) {
      return(range_list[[key2]])
    }
    
    NULL
  }
  
  # ------------------------------------------------------------
  # Domain-specific observed X ranges
  # Default behavior: each adversity domain gets its own observed range.
  # Top and bottom rows share the same X range within domain.
  # ------------------------------------------------------------
  
  domain_x_ranges <- list()
  
  for (dom_lab in domain_order) {
    
    dom_specs <- panel_specs %>%
      dplyr::filter(.data$domain_label == dom_lab)
    
    if (nrow(dom_specs) == 0) {
      stop("No panels found for domain_label: ", dom_lab)
    }
    
    manual_dom_range <- get_named_range(
      range_list = x_ranges,
      key1 = dom_lab,
      key2 = dom_specs$domain[1]
    )
    
    if (!is.null(manual_dom_range)) {
      
      if (!is.numeric(manual_dom_range) || length(manual_dom_range) != 2) {
        stop("Each manual x-range in x_ranges must be numeric length 2.")
      }
      
      domain_x_ranges[[dom_lab]] <- manual_dom_range
      
    } else if (!is.null(x_range)) {
      
      domain_x_ranges[[dom_lab]] <- x_range
      
    } else {
      
      observed_x <- unlist(
        lapply(dom_specs$plot_df, function(df) {
          if (!"x_value" %in% names(df)) {
            stop("Each plot_df must contain x_value.")
          }
          df$x_value
        })
      )
      
      observed_x <- observed_x[is.finite(observed_x)]
      
      if (length(observed_x) == 0) {
        stop("No observed x_value data found for domain_label: ", dom_lab)
      }
      
      domain_x_ranges[[dom_lab]] <- range(observed_x * x_multiplier, na.rm = TRUE)
    }
  }
  
  # ------------------------------------------------------------
  # Polynomial helper
  # Models are estimated using mean-centered raw X.
  # Predictions are plotted back on the original 0-100-ish metric.
  # ------------------------------------------------------------
  
  build_polynomial_line <- function(dat, degree, x_min_raw, x_max_raw, x_multiplier) {
    
    poly_dat <- dat %>%
      dplyr::filter(
        !is.na(x_raw),
        !is.na(y_plot),
        x_raw >= x_min_raw,
        x_raw <= x_max_raw
      )
    
    if (nrow(poly_dat) < 5) {
      return(tibble::tibble())
    }
    
    x_mean <- mean(poly_dat$x_raw, na.rm = TRUE)
    
    poly_dat <- poly_dat %>%
      dplyr::mutate(
        x_c = x_raw - x_mean,
        x_c2 = x_c^2,
        x_c3 = x_c^3
      )
    
    if (degree == 2) {
      fit <- stats::lm(y_plot ~ x_c + x_c2, data = poly_dat)
    }
    
    if (degree == 3) {
      fit <- stats::lm(y_plot ~ x_c + x_c2 + x_c3, data = poly_dat)
    }
    
    pred_x_raw <- seq(
      from = x_min_raw,
      to = x_max_raw,
      length.out = 300
    )
    
    pred_dat <- tibble::tibble(
      x_raw = pred_x_raw,
      x_plot = pred_x_raw * x_multiplier,
      x_c = pred_x_raw - x_mean
    )
    
    pred_dat <- pred_dat %>%
      dplyr::mutate(
        x_c2 = x_c^2,
        x_c3 = x_c^3
      )
    
    pred_dat$y_plot <- as.numeric(
      stats::predict(fit, newdata = pred_dat)
    )
    
    pred_dat <- pred_dat %>%
      dplyr::mutate(
        model = ifelse(degree == 2, "Quadratic", "Cubic"),
        x_mean_raw = x_mean,
        x_mean_0_100 = x_mean * x_multiplier
      )
    
    pred_dat
  }
  
  # ------------------------------------------------------------
  # Build data for each panel
  # ------------------------------------------------------------
  
  build_panel_data <- function(plot_df, spline_par, domain_label) {
    
    needed_par <- c("knot_x", "knot_y", "b1", "b2", "ref_x", "ref_y")
    missing_par <- setdiff(needed_par, names(spline_par))
    
    if (length(missing_par) > 0) {
      stop("spline_par is missing: ", paste(missing_par, collapse = ", "))
    }
    
    if (!all(c("x_value", "y_value") %in% names(plot_df))) {
      stop("Each plot_df must contain x_value and y_value.")
    }
    
    x_range_this <- domain_x_ranges[[domain_label]]
    x_min_raw <- x_range_this[1] / x_multiplier
    x_max_raw <- x_range_this[2] / x_multiplier
    
    dat <- plot_df %>%
      dplyr::mutate(
        x_raw = x_value,
        x_plot = x_value * x_multiplier,
        y_plot = y_value,
        knot_side = dplyr::case_when(
          x_raw <= spline_par$knot_x ~ "Left of knot",
          x_raw >  spline_par$knot_x ~ "Right of knot",
          TRUE ~ NA_character_
        )
      ) %>%
      dplyr::filter(
        !is.na(x_raw),
        !is.na(x_plot),
        !is.na(y_plot)
      )
    
    if (nrow(dat) == 0) {
      stop("One plot_df has no complete x_value/y_value rows.")
    }
    
    knot_x_raw   <- spline_par$knot_x
    knot_y       <- spline_par$knot_y
    knot_y_left  <- if ("knot_y_left" %in% names(spline_par)) spline_par$knot_y_left else knot_y
    knot_y_right <- if ("knot_y_right" %in% names(spline_par)) spline_par$knot_y_right else knot_y
    knot_x_plot  <- knot_x_raw * x_multiplier
    
    left_line <- tibble::tibble()
    right_line <- tibble::tibble()
    
    if (knot_x_raw > x_min_raw) {
      
      left_x_raw <- seq(
        from = x_min_raw,
        to = min(knot_x_raw, x_max_raw),
        length.out = 150
      )
      
      left_line <- tibble::tibble(
        x_raw = left_x_raw,
        x_plot = left_x_raw * x_multiplier,
        y_plot = knot_y_left + spline_par$b1 * (left_x_raw - knot_x_raw),
        segment = "Low-range two-lines"
      )
    }
    
    if (knot_x_raw < x_max_raw) {
      
      right_x_raw <- seq(
        from = max(knot_x_raw, x_min_raw),
        to = x_max_raw,
        length.out = 150
      )
      
      right_line <- tibble::tibble(
        x_raw = right_x_raw,
        x_plot = right_x_raw * x_multiplier,
        y_plot = knot_y_right + spline_par$b2 * (right_x_raw - knot_x_raw),
        segment = "High-range two-lines"
      )
    }
    
    line_df <- dplyr::bind_rows(left_line, right_line)
    
    if (nrow(line_df) == 0) {
      stop("Could not reconstruct one spline line. Check knot_x and observed x-range.")
    }
    
    quadratic_df <- tibble::tibble()
    cubic_df <- tibble::tibble()
    
    if (isTRUE(plot_quadratic)) {
      quadratic_df <- build_polynomial_line(
        dat = dat,
        degree = 2,
        x_min_raw = x_min_raw,
        x_max_raw = x_max_raw,
        x_multiplier = x_multiplier
      )
    }
    
    if (isTRUE(plot_cubic)) {
      cubic_df <- build_polynomial_line(
        dat = dat,
        degree = 3,
        x_min_raw = x_min_raw,
        x_max_raw = x_max_raw,
        x_multiplier = x_multiplier
      )
    }
    
    list(
      dat = dat,
      line_df = line_df,
      quadratic_df = quadratic_df,
      cubic_df = cubic_df,
      knot_x_raw = knot_x_raw,
      knot_x_plot = knot_x_plot,
      knot_y = knot_y,
      knot_y_left = knot_y_left,
      knot_y_right = knot_y_right,
      knot_y_gap = knot_y_right - knot_y_left,
      x_range = x_range_this,
      x_min_raw = x_min_raw,
      x_max_raw = x_max_raw
    )
  }
  
  panel_objects <- vector("list", nrow(panel_specs))
  
  for (i in seq_len(nrow(panel_specs))) {
    
    panel_objects[[i]] <- build_panel_data(
      plot_df = panel_specs$plot_df[[i]],
      spline_par = panel_specs$spline_par[[i]],
      domain_label = panel_specs$domain_label[i]
    )
  }
  
  panel_specs$.panel_id <- seq_len(nrow(panel_specs))
  
  # ------------------------------------------------------------
  # Shared row-specific y-ranges
  # Uses model lines only, not observed points.
  # Manual y_ranges overrides this.
  # ------------------------------------------------------------
  
  row_y_ranges <- list()
  
  for (out_lab in outcome_order) {
    
    row_specs <- panel_specs %>%
      dplyr::filter(.data$outcome_label == out_lab)
    
    if (nrow(row_specs) == 0) {
      stop("No panels found for outcome_label: ", out_lab)
    }
    
    manual_range <- get_named_range(
      range_list = y_ranges,
      key1 = out_lab,
      key2 = row_specs$outcome[1]
    )
    
    if (!is.null(manual_range)) {
      
      if (!is.numeric(manual_range) || length(manual_range) != 2) {
        stop("Each manual y-range must be numeric length 2.")
      }
      
      row_y_ranges[[out_lab]] <- manual_range
      
    } else {
      
      row_model_y <- unlist(
        lapply(row_specs$.panel_id, function(id) {
          c(
            panel_objects[[id]]$line_df$y_plot,
            panel_objects[[id]]$quadratic_df$y_plot,
            panel_objects[[id]]$cubic_df$y_plot
          )
        })
      )
      
      y_rng <- range(row_model_y, na.rm = TRUE)
      y_span <- diff(y_rng)
      
      if (!is.finite(y_span) || y_span == 0) {
        y_span <- 0.10
      }
      
      row_y_ranges[[out_lab]] <- c(
        y_rng[1] - y_span * y_pad_lower,
        y_rng[2] + y_span * y_pad_upper
      )
    }
  }
  
  # ------------------------------------------------------------
  # Top legend
  # ------------------------------------------------------------
  
  build_top_legend <- function() {
    
    legend_levels <- c("Low-range two-lines", "High-range two-lines")
    legend_colors <- c(
      "Low-range two-lines" = slope1_color,
      "High-range two-lines" = slope2_color
    )
    
    if (isTRUE(plot_quadratic)) {
      legend_levels <- c(legend_levels, "Quadratic function")
      legend_colors <- c(
        legend_colors,
        "Quadratic function" = quadratic_color
      )
    }
    
    if (isTRUE(plot_cubic)) {
      legend_levels <- c(legend_levels, "Cubic function")
      legend_colors <- c(
        legend_colors,
        "Cubic function" = cubic_color
      )
    }
    
    legend_df <- tibble::tibble(
      x = rep(c(0, 1), length(legend_levels)),
      y = rep(seq_along(legend_levels), each = 2),
      model = factor(rep(legend_levels, each = 2), levels = legend_levels)
    )
    
    p_leg <- ggplot2::ggplot(
      legend_df,
      ggplot2::aes(x = x, y = y, color = model)
    ) +
      ggplot2::geom_line(linewidth = legend_linewidth) +
      ggplot2::scale_color_manual(values = legend_colors, drop = FALSE) +
      ggplot2::guides(
        color = ggplot2::guide_legend(
          title = NULL,
          nrow = 1,
          byrow = TRUE,
          override.aes = list(linewidth = legend_linewidth)
        )
      ) +
      ggplot2::theme_void(base_family = base_family) +
      ggplot2::theme(
        legend.position = legend_position,
        legend.text = ggplot2::element_text(size = legend_text_size),
        legend.key.width = grid::unit(1.35, "cm"),
        legend.margin = ggplot2::margin(0, 0, 0, 0)
      )
    
    cowplot::get_legend(p_leg)
  }
  
  # ------------------------------------------------------------
  # Main panel plot
  # ------------------------------------------------------------
  
  build_panel_plot <- function(spec_row, panel_obj, row_index, col_index) {
    
    y_range <- row_y_ranges[[spec_row$outcome_label]]
    x_range_this <- panel_obj$x_range
    
    dat <- panel_obj$dat
    line_df <- panel_obj$line_df
    
    dat_points <- dat %>%
      dplyr::filter(
        x_plot >= x_range_this[1],
        x_plot <= x_range_this[2],
        y_plot >= y_range[1],
        y_plot <= y_range[2]
      )
    
    if (!is.null(n_points) && isTRUE(show_points)) {
      set.seed(seed)
      dat_points <- dat_points %>%
        dplyr::slice_sample(n = min(n_points, nrow(dat_points)))
    }
    
    show_left_axis <- col_index == 1
    
    y_lab <- if (show_left_axis) spec_row$outcome_label else ""
    panel_title <- if (row_index == 1) spec_row$domain_label else ""
    
    knot_label_y <- y_range[2] - diff(y_range) * knot_label_top_nudge
    knot_label <- round(panel_obj$knot_x_plot, knot_label_digits)
    
    p_main <- ggplot2::ggplot() +
      ggplot2::geom_vline(
        xintercept = panel_obj$knot_x_plot,
        linewidth = knot_linewidth,
        color = knot_line_color
      )
    
    if (isTRUE(show_points) && nrow(dat_points) > 0 && point_mode == "gray") {
      p_main <- p_main +
        ggplot2::geom_point(
          data = dat_points,
          ggplot2::aes(x = x_plot, y = y_plot),
          color = gray_point_color,
          alpha = point_alpha,
          size = point_size
        )
    }
    
    if (isTRUE(show_points) && nrow(dat_points) > 0 && point_mode == "colored") {
      p_main <- p_main +
        ggplot2::geom_point(
          data = dat_points,
          ggplot2::aes(x = x_plot, y = y_plot, color = knot_side),
          alpha = point_alpha,
          size = point_size
        ) +
        ggplot2::scale_color_manual(
          values = c(
            "Left of knot" = left_point_color,
            "Right of knot" = right_point_color
          )
        )
    }
    
    p_main <- p_main +
      ggplot2::geom_line(
        data = dplyr::filter(line_df, segment == "Low-range two-lines"),
        ggplot2::aes(x = x_plot, y = y_plot),
        color = slope1_color,
        linewidth = slope_linewidth
      ) +
      ggplot2::geom_line(
        data = dplyr::filter(line_df, segment == "High-range two-lines"),
        ggplot2::aes(x = x_plot, y = y_plot),
        color = slope2_color,
        linewidth = slope_linewidth
      )
    
    if (isTRUE(plot_quadratic) && nrow(panel_obj$quadratic_df) > 0) {
      p_main <- p_main +
        ggplot2::geom_line(
          data = panel_obj$quadratic_df,
          ggplot2::aes(x = x_plot, y = y_plot),
          color = quadratic_color,
          linewidth = polynomial_linewidth,
          linetype = polynomial_linetype
        )
    }
    
    if (isTRUE(plot_cubic) && nrow(panel_obj$cubic_df) > 0) {
      p_main <- p_main +
        ggplot2::geom_line(
          data = panel_obj$cubic_df,
          ggplot2::aes(x = x_plot, y = y_plot),
          color = cubic_color,
          linewidth = polynomial_linewidth,
          linetype = polynomial_linetype
        )
    }
    
    p_main <- p_main +
      ggplot2::coord_cartesian(
        xlim = x_range_this,
        ylim = y_range,
        clip = "off"
      ) +
      ggplot2::labs(
        x = "",
        y = y_lab,
        title = panel_title
      ) +
      ggplot2::theme_minimal(
        base_size = base_size,
        base_family = base_family
      ) +
      ggplot2::theme(
        panel.grid.minor = ggplot2::element_blank(),
        panel.grid.major = ggplot2::element_line(
          linewidth = 0.22,
          color = "gray90"
        ),
        axis.title.x = ggplot2::element_blank(),
        axis.title.y = ggplot2::element_text(
          size = axis_title_size,
          face = "bold"
        ),
        axis.text = ggplot2::element_text(
          size = axis_text_size,
          color = "black"
        ),
        axis.ticks = ggplot2::element_line(
          color = "black",
          linewidth = 0.30
        ),
        plot.title = ggplot2::element_text(
          size = strip_title_size,
          hjust = 0.5,
          face = "bold"
        ),
        legend.position = "none",
        plot.margin = plot_margin
      )
    
    if (!show_left_axis && !isTRUE(show_y_axis_text_all)) {
      p_main <- p_main +
        ggplot2::theme(
          axis.text.y = ggplot2::element_blank(),
          axis.ticks.y = ggplot2::element_blank()
        )
    }
    
    if (isTRUE(show_knot_label)) {
      p_main <- p_main +
        ggplot2::annotate(
          "label",
          x = panel_obj$knot_x_plot,
          y = knot_label_y,
          label = knot_label,
          size = knot_label_size,
          family = base_family,
          fontface = "bold",
          fill = knot_label_fill,
          color = "black",
          label.size = 0.20
        )
    }
    
    p_main
  }
  
  # ------------------------------------------------------------
  # Single bottom-row histogram per adversity domain
  # ------------------------------------------------------------
  
  build_hist_plot <- function(domain_label_value) {
    
    spec_match <- panel_specs %>%
      dplyr::filter(.data$domain_label == domain_label_value) %>%
      dplyr::slice(1)
    
    panel_id <- spec_match$.panel_id[1]
    x_range_this <- panel_objects[[panel_id]]$x_range
    
    dat_hist <- panel_objects[[panel_id]]$dat %>%
      dplyr::filter(
        x_plot >= x_range_this[1],
        x_plot <= x_range_this[2]
      )
    
    p_hist <- ggplot2::ggplot(dat_hist, ggplot2::aes(x = x_plot)) +
      ggplot2::geom_histogram(
        bins = hist_bins,
        fill = hist_fill,
        color = hist_color,
        alpha = hist_alpha
      ) +
      ggplot2::coord_cartesian(
        xlim = x_range_this,
        clip = "off"
      ) +
      ggplot2::labs(x = "", y = "") +
      ggplot2::theme_minimal(
        base_size = base_size,
        base_family = base_family
      ) +
      ggplot2::theme(
        panel.grid = ggplot2::element_blank(),
        axis.title = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_text(
          size = axis_text_size,
          color = "black"
        ),
        axis.text.y = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_line(
          color = "black",
          linewidth = 0.30
        ),
        axis.ticks.y = ggplot2::element_blank(),
        plot.margin = ggplot2::margin(0, 6, 2, 6)
      )
    
    p_hist
  }
  
  # ------------------------------------------------------------
  # Assemble 4 columns x 2 outcome rows + one bottom histogram row
  # ------------------------------------------------------------
  
  row_plots <- list()
  panel_plot_list <- list()
  
  for (r in seq_along(outcome_order)) {
    
    out_lab <- outcome_order[r]
    this_row <- list()
    
    for (c in seq_along(domain_order)) {
      
      dom_lab <- domain_order[c]
      
      spec_match <- panel_specs %>%
        dplyr::filter(
          .data$outcome_label == out_lab,
          .data$domain_label == dom_lab
        )
      
      if (nrow(spec_match) != 1) {
        stop(
          "Expected exactly one panel for domain_label = '",
          dom_lab,
          "' and outcome_label = '",
          out_lab,
          "'. Found: ",
          nrow(spec_match)
        )
      }
      
      panel_id <- spec_match$.panel_id[1]
      
      p <- build_panel_plot(
        spec_row = spec_match[1, ],
        panel_obj = panel_objects[[panel_id]],
        row_index = r,
        col_index = c
      )
      
      this_row[[c]] <- p
      panel_plot_list[[paste(out_lab, dom_lab, sep = " | ")]] <- p
    }
    
    row_plots[[r]] <- cowplot::plot_grid(
      plotlist = this_row,
      nrow = 1,
      align = "hv",
      axis = "tblr"
    )
  }
  
  if (isTRUE(show_x_hist)) {
    
    hist_row <- cowplot::plot_grid(
      plotlist = lapply(domain_order, build_hist_plot),
      nrow = 1,
      align = "h",
      axis = "tb"
    )
    
    body_plot <- cowplot::plot_grid(
      plotlist = c(row_plots, list(hist_row)),
      ncol = 1,
      align = "v",
      axis = "lr",
      rel_heights = c(rep(1, length(row_plots)), hist_height)
    )
    
  } else {
    
    body_plot <- cowplot::plot_grid(
      plotlist = row_plots,
      ncol = 1,
      align = "v",
      axis = "lr",
      rel_heights = rep(1, length(row_plots))
    )
  }
  
  if (isTRUE(show_legend)) {
    
    top_legend <- build_top_legend()
    
    final_plot <- cowplot::plot_grid(
      top_legend,
      body_plot,
      ncol = 1,
      rel_heights = c(legend_height, 1)
    )
    
  } else {
    
    final_plot <- body_plot
  }
  
  # ------------------------------------------------------------
  # Diagnostics
  # ------------------------------------------------------------
  
  if (isTRUE(print_diagnostics)) {
    
    cat("\n-----------------------------\n")
    cat("Spline panel-grid diagnostics\n")
    cat("-----------------------------\n")
    cat("X multiplier:", x_multiplier, "\n")
    cat("Polynomial overlay:", polynomial, "\n\n")
    
    cat("Domain-specific X ranges, plotted metric:\n")
    print(domain_x_ranges)
    
    cat("\nRow-specific Y ranges:\n")
    print(row_y_ranges)
    
    cat("\nKnot values plotted on observed 0-100-style metric:\n")
    
    knot_diag <- dplyr::bind_rows(
      lapply(seq_len(nrow(panel_specs)), function(i) {
        tibble::tibble(
          domain = panel_specs$domain_label[i],
          outcome = panel_specs$outcome_label[i],
          knot_raw = panel_objects[[i]]$knot_x_raw,
          knot_plot_metric = panel_objects[[i]]$knot_x_plot
        )
      })
    )
    
    print(knot_diag)
    
    if (polynomial != "none") {
      
      cat("\nPolynomial centering values, shown on plotted metric:\n")
      
      poly_diag <- dplyr::bind_rows(
        lapply(seq_len(nrow(panel_specs)), function(i) {
          
          q_mean <- NA_real_
          c_mean <- NA_real_
          
          if (nrow(panel_objects[[i]]$quadratic_df) > 0) {
            q_mean <- unique(panel_objects[[i]]$quadratic_df$x_mean_0_100)[1]
          }
          
          if (nrow(panel_objects[[i]]$cubic_df) > 0) {
            c_mean <- unique(panel_objects[[i]]$cubic_df$x_mean_0_100)[1]
          }
          
          tibble::tibble(
            domain = panel_specs$domain_label[i],
            outcome = panel_specs$outcome_label[i],
            quadratic_center_plot_metric = q_mean,
            cubic_center_plot_metric = c_mean
          )
        })
      )
      
      print(poly_diag)
    }
    
    cat("-----------------------------\n\n")
  }
  
  list(
    plot = final_plot,
    panels = panel_plot_list,
    panel_specs = panel_specs,
    panel_objects = panel_objects,
    x_ranges = domain_x_ranges,
    y_ranges = row_y_ranges
  )
}

##-----------------------------------------------------------------------------
## make the plot 
##-----------------------------------------------------------------------------

panel_specs <- tibble::tibble(
  domain = c(
    "HSES", "Neighborhood Danger", "School Adversity", "Family Conflict",
    "HSES", "Neighborhood Danger", "School Adversity", "Family Conflict"
  ),
  
  outcome = c(
    "mEmoSup", "mEmoSup", "mEmoSup", "mEmoSup",
    "mReApp",  "mReApp",  "mReApp",  "mReApp"
  ),
  
  domain_label = c(
    "HSES", "Neighborhood Danger", "School Adversity", "Family Conflict",
    "HSES", "Neighborhood Danger", "School Adversity", "Family Conflict"
  ),
  
  outcome_label = c(
    "Emotional Suppression", "Emotional Suppression",
    "Emotional Suppression", "Emotional Suppression",
    "Reappraisal", "Reappraisal",
    "Reappraisal", "Reappraisal"
  ),
  
  plot_df = list(
    HSES_SUPP,
    NBDANG_SUPP,
    SCHADV_SUPP,
    FCON_SUPP,
    HSES_REAP,
    NBDANG_REAP,
    SCHADV_REAP,
    FCON_REAP
  ),
  
  spline_par = list(
    HSES_SUPP_par,
    NBDANG_SUPP_par,
    SCHADV_SUPP_par,
    FCON_SUPP_par,
    HSES_REAP_par,
    NBDANG_REAP_par,
    SCHADV_REAP_par,
    FCON_REAP_par
  )
)


SplineGrid <- make_spline_panel_grid(
  panel_specs = panel_specs,
  show_points = FALSE,
  polynomial = "none"
)

SplineGrid$plot

## fit quadratic effect 
SplineGrid_quad <- make_spline_panel_grid(
  panel_specs = panel_specs,
  show_points = FALSE,
  polynomial = "quadratic"
)

SplineGrid_quad$plot

## fit cubic effect 
SplineGrid_both <- make_spline_panel_grid(
  panel_specs = panel_specs,
  show_points = FALSE,
  polynomial = "both"
)

SplineGrid_both$plot

## save plot 

ggsave(
  filename = "two_lines_panel_plot.png",
  plot = SplineGrid$plot,
  width = 13,
  height = 7,
  dpi = 600,
  bg = "white"
)
