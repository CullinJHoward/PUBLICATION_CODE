## liBRARY 
library(ggplot2)

## set environment 
setwd("C:\\Users\\cjh37695\\Dropbox\\HORMESIS_METHODS\\ANALYSIS\\") #work

# pick-a-point data 
df <- read.csv("PICK_A_POINT_DATA.csv")
# participant data 

df_vars <- read.csv("ABCD_HORM_METH_SEM_STRUCTRUAL_MOD_5.25.26.csv")


##------------------------------------------------------------------------------
                          # plotting interactions # 
## -----------------------------------------------------------------------------

## trajectories 

PieceOfTraj <- function(
    data,
    traj_group,
    traj_int,
    traj_slope,
    uncon_traj,
    comp_traj = NULL,
    comp_colors = NULL,
    sig_name = NULL,
    mark_sig = FALSE,
    add_sig_leg = FALSE,
    sig_text_color = "#B22222",
    sig_text_size = 4,
    sig_nudge_x = NULL,
    sig_nudge_y = NULL,
    mod_name = NULL,
    mod_traj = NULL,
    reference_mod_level = NULL,
    conditional_low_label = "Conditional Low",
    time_values,
    x_label = "Time",
    y_label = "Predicted Outcome",
    color_label = "Trajectory",
    title = "Predicted Trajectories",
    show_points = TRUE,
    point_size = 4,
    uncon_linewidth = 2.5,
    comp_linewidth = 2.5,
    show_uncon_legend = TRUE,
    base_size = 12,
    title_size = 13,
    axis_title_size = 12,
    axis_text_size = 15,
    legend_title_size = 11,
    legend_text_size = 10,
    strip_text_size = 11,
    make_individual_plots = TRUE,
    split_moderated_individual_plots = TRUE,
    facet_integrated_plot = TRUE,
    print_integrated_plot = TRUE,
    print_individual_plots = TRUE,
    return_plot_data = TRUE,
    transparent_background = TRUE,
    show_gridlines = TRUE,
    gridline_color = "grey78",
    gridline_linewidth = 0.35,
    gridline_linetype = "solid",
    show_background_bands = TRUE,
    background_band_n = 8,
    background_band_low = "grey50",
    background_band_high = "grey98",
    background_band_alpha = 0.75,
    y_axis_digits = 2
) {
  
  if (!is.data.frame(data)) stop("`data` must be a data frame.")
  
  required_name_args <- c("traj_group", "traj_int", "traj_slope")
  arg_vals <- list(traj_group = traj_group, traj_int = traj_int, traj_slope = traj_slope)
  bad_args <- required_name_args[
    !vapply(arg_vals, function(x) is.character(x) && length(x) == 1 && !is.na(x), logical(1))
  ]
  if (length(bad_args) > 0) {
    stop("These arguments must be single character strings naming columns in `data`: ",
         paste(bad_args, collapse = ", "))
  }
  
  if (!is.null(mod_name) &&
      (!is.character(mod_name) || length(mod_name) != 1 || is.na(mod_name))) {
    stop("`mod_name` must be NULL or a single character string naming a column in `data`.")
  }
  
  if (!is.null(reference_mod_level) &&
      (!is.character(reference_mod_level) || length(reference_mod_level) != 1 || is.na(reference_mod_level))) {
    stop("`reference_mod_level` must be NULL or a single character string.")
  }
  
  if (!is.character(conditional_low_label) || length(conditional_low_label) != 1 || is.na(conditional_low_label)) {
    stop("`conditional_low_label` must be a single character string.")
  }
  
  if (!is.null(comp_colors)) {
    if (!is.character(comp_colors) || is.null(names(comp_colors)) || any(names(comp_colors) == "")) {
      stop("`comp_colors` must be a named character vector of colors.")
    }
  }
  
  if (!is.null(sig_name) &&
      (!is.character(sig_name) || length(sig_name) != 1 || is.na(sig_name))) {
    stop("`sig_name` must be NULL or a single character string naming a numeric column in `data`.")
  }
  
  logical_args <- list(
    mark_sig = mark_sig,
    add_sig_leg = add_sig_leg,
    show_uncon_legend = show_uncon_legend,
    show_points = show_points,
    make_individual_plots = make_individual_plots,
    split_moderated_individual_plots = split_moderated_individual_plots,
    facet_integrated_plot = facet_integrated_plot,
    print_integrated_plot = print_integrated_plot,
    print_individual_plots = print_individual_plots,
    return_plot_data = return_plot_data,
    transparent_background = transparent_background,
    show_gridlines = show_gridlines,
    show_background_bands = show_background_bands
  )
  bad_logical <- names(logical_args)[
    !vapply(logical_args, function(x) is.logical(x) && length(x) == 1 && !is.na(x), logical(1))
  ]
  if (length(bad_logical) > 0) {
    stop("These arguments must be TRUE or FALSE: ", paste(bad_logical, collapse = ", "))
  }
  
  needed_cols <- c(traj_group, traj_int, traj_slope, mod_name)
  if (!is.null(sig_name)) needed_cols <- c(needed_cols, sig_name)
  needed_cols <- needed_cols[!is.na(needed_cols) & !sapply(needed_cols, is.null)]
  missing_cols <- setdiff(needed_cols, names(data))
  if (length(missing_cols) > 0) {
    stop("Missing required column(s) in `data`: ", paste(missing_cols, collapse = ", "))
  }
  
  if (!is.null(sig_name) && !is.numeric(data[[sig_name]])) stop("`sig_name` must name a numeric column.")
  if (mark_sig && is.null(sig_name)) stop("When `mark_sig = TRUE`, `sig_name` must be provided.")
  
  if (!is.character(uncon_traj) || length(uncon_traj) != 1 || is.na(uncon_traj)) {
    stop("`uncon_traj` must be a single character string matching one value in `traj_group`.")
  }
  
  if (!is.null(comp_traj)) {
    if (!is.character(comp_traj) || length(comp_traj) < 1 || any(is.na(comp_traj))) {
      stop("`comp_traj` must be NULL or a character vector of trajectory labels.")
    }
    if (anyDuplicated(comp_traj)) stop("`comp_traj` contains duplicate labels.")
    if (uncon_traj %in% comp_traj) stop("`uncon_traj` should not also appear in `comp_traj`.")
  }
  
  if (is.null(mod_name)) {
    if (!is.null(mod_traj)) stop("`mod_traj` should only be supplied when `mod_name` is used.")
    if (!is.null(reference_mod_level)) stop("`reference_mod_level` should only be supplied when `mod_name` is used.")
  } else {
    if (is.null(mod_traj)) stop("`mod_traj` is required when `mod_name` is supplied.")
    if (!is.character(mod_traj) || length(mod_traj) < 1 || any(is.na(mod_traj))) {
      stop("`mod_traj` must be a character vector of moderator-group labels.")
    }
    if (anyDuplicated(mod_traj)) stop("`mod_traj` contains duplicate labels.")
  }
  
  if (missing(time_values)) stop("`time_values` is required.")
  if (!is.numeric(time_values) || length(time_values) < 2 || any(is.na(time_values))) {
    stop("`time_values` must be a numeric vector with at least two non-missing values.")
  }
  
  numeric_single_args <- list(
    point_size = point_size,
    uncon_linewidth = uncon_linewidth,
    comp_linewidth = comp_linewidth,
    base_size = base_size,
    title_size = title_size,
    axis_title_size = axis_title_size,
    axis_text_size = axis_text_size,
    legend_title_size = legend_title_size,
    legend_text_size = legend_text_size,
    strip_text_size = strip_text_size,
    sig_text_size = sig_text_size,
    gridline_linewidth = gridline_linewidth,
    background_band_alpha = background_band_alpha,
    y_axis_digits = y_axis_digits
  )
  bad_numeric <- names(numeric_single_args)[
    !vapply(numeric_single_args, function(x) is.numeric(x) && length(x) == 1 && !is.na(x), logical(1))
  ]
  if (length(bad_numeric) > 0) {
    stop("These arguments must be single non-missing numeric values: ", paste(bad_numeric, collapse = ", "))
  }
  
  if (background_band_alpha < 0 || background_band_alpha > 1) {
    stop("`background_band_alpha` must be between 0 and 1.")
  }
  if (is.null(background_band_n)) background_band_n <- 8
  if (!is.numeric(background_band_n) || length(background_band_n) != 1 ||
      is.na(background_band_n) || background_band_n < 2) {
    stop("`background_band_n` must be NULL or a single numeric value >= 2.")
  }
  background_band_n <- as.integer(round(background_band_n))
  y_axis_digits <- as.integer(round(y_axis_digits))
  if (y_axis_digits < 0) stop("`y_axis_digits` must be >= 0.")
  
  if (!is.character(background_band_low) || length(background_band_low) != 1 || is.na(background_band_low)) {
    stop("`background_band_low` must be a single color string.")
  }
  if (!is.character(background_band_high) || length(background_band_high) != 1 || is.na(background_band_high)) {
    stop("`background_band_high` must be a single color string.")
  }
  if (!is.null(sig_nudge_x) && (!is.numeric(sig_nudge_x) || length(sig_nudge_x) != 1 || is.na(sig_nudge_x))) {
    stop("`sig_nudge_x` must be NULL or a single numeric value.")
  }
  if (!is.null(sig_nudge_y) && (!is.numeric(sig_nudge_y) || length(sig_nudge_y) != 1 || is.na(sig_nudge_y))) {
    stop("`sig_nudge_y` must be NULL or a single numeric value.")
  }
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required. Please install it.", call. = FALSE)
  }
  
  .sig_to_symbol <- function(p) {
    ifelse(is.na(p), "",
           ifelse(p <= 0.001, "***",
                  ifelse(p <= 0.01, "**",
                         ifelse(p <= 0.05, "*",
                                ifelse(p <= 0.10, "^", "")))))
  }
  
  .format_axis_labels <- function(x, digits = 2) {
    x <- round(x, digits)
    out <- formatC(x, format = "f", digits = digits)
    out <- sub("0+$", "", out)
    out <- sub("\\.$", "", out)
    out
  }
  
  .safe_plot_name <- function(x) {
    x <- gsub("\\+", "plus", x)
    x <- gsub("-", "minus", x)
    x <- gsub("[^A-Za-z0-9]+", "_", x)
    x <- gsub("_+", "_", x)
    x <- gsub("^_|_$", "", x)
    x
  }
  
  .resolve_named_colors <- function(labels, color_vec, color_arg_name = "colors") {
    if (is.null(labels) || length(labels) == 0) return(NULL)
    if (is.null(color_vec)) return(NULL)
    exact_idx <- match(labels, names(color_vec))
    if (all(!is.na(exact_idx))) return(stats::setNames(unname(color_vec[exact_idx]), labels))
    lower_idx <- match(tolower(labels), tolower(names(color_vec)))
    if (all(!is.na(lower_idx))) return(stats::setNames(unname(color_vec[lower_idx]), labels))
    missing_labels <- labels[is.na(lower_idx)]
    stop("These trajectory label(s) are missing from `", color_arg_name, "`: ",
         paste(missing_labels, collapse = ", "), call. = FALSE)
  }
  
  .make_y_grid_and_bands <- function(y_limits,
                                     background_band_n = 8,
                                     background_band_low = "grey50",
                                     background_band_high = "grey98",
                                     background_band_alpha = 0.75) {
    if (!is.numeric(y_limits) || length(y_limits) != 2 || any(is.na(y_limits))) {
      stop("Internal plotting error: invalid `y_limits`.")
    }
    if (y_limits[1] == y_limits[2]) y_limits <- y_limits + c(-0.5, 0.5)
    background_band_n <- as.integer(round(background_band_n))
    y_grid_breaks <- seq(from = y_limits[1], to = y_limits[2], length.out = background_band_n + 1)
    n_bands <- length(y_grid_breaks) - 1
    band_cols <- grDevices::colorRampPalette(c(background_band_low, background_band_high))(n_bands)
    band_cols <- grDevices::adjustcolor(band_cols, alpha.f = background_band_alpha)
    band_df <- data.frame(
      ymin = y_grid_breaks[-length(y_grid_breaks)],
      ymax = y_grid_breaks[-1],
      .band_fill = band_cols,
      stringsAsFactors = FALSE
    )
    list(y_grid_breaks = y_grid_breaks, band_df = band_df)
  }
  
  .sig_key <- "Significance key: ^ p <= .10, * p <= .05, ** p <= .01, *** p <= .001"
  bg_fill <- if (transparent_background) "transparent" else "white"
  
  df <- data
  df$.traj_group <- as.character(df[[traj_group]])
  df$.traj_int <- df[[traj_int]]
  df$.traj_slope <- df[[traj_slope]]
  df$.moderator <- if (is.null(mod_name)) "Overall" else as.character(df[[mod_name]])
  df$.sig <- if (is.null(sig_name)) NA_real_ else df[[sig_name]]
  
  if (!is.numeric(df$.traj_int)) stop("`traj_int` must name a numeric column.")
  if (!is.numeric(df$.traj_slope)) stop("`traj_slope` must name a numeric column.")
  if (any(is.na(df$.traj_group))) stop("`traj_group` contains missing values.")
  if (any(is.na(df$.traj_int))) stop("`traj_int` contains missing values.")
  if (any(is.na(df$.traj_slope))) stop("`traj_slope` contains missing values.")
  
  if (is.null(mod_name)) {
    df$.moderator[is.na(df$.moderator)] <- "Overall"
  } else if (any(is.na(df$.moderator))) {
    stop("Moderator labels contain missing values. For moderated plots, every plotted row needs a moderator label.")
  }
  
  raw_uncon_rows <- df[df$.traj_group == uncon_traj, , drop = FALSE]
  if (nrow(raw_uncon_rows) == 0) {
    stop("`uncon_traj` was not found in `traj_group`. Available values are: ",
         paste(unique(df$.traj_group), collapse = ", "))
  }
  
  if (nrow(raw_uncon_rows) == 1) {
    uncon_row <- raw_uncon_rows[1, , drop = FALSE]
    conditional_low_rows <- raw_uncon_rows[FALSE, , drop = FALSE]
  } else {
    if (is.null(mod_name)) {
      stop("More than one row matched `uncon_traj`. Exactly one unconditional trajectory is allowed when `mod_name` is NULL.")
    }
    if (is.null(reference_mod_level)) {
      if ("Mean" %in% raw_uncon_rows$.moderator) {
        reference_mod_level <- "Mean"
      } else {
        stop(
          "More than one row matched `uncon_traj`. Supply `reference_mod_level` to identify the single true reference row."
        )
      }
    }
    ref_idx <- raw_uncon_rows$.moderator == reference_mod_level
    if (sum(ref_idx) != 1) {
      stop(
        "`reference_mod_level` must identify exactly one `uncon_traj` row. Found ",
        sum(ref_idx), " row(s) for moderator level: ", reference_mod_level
      )
    }
    uncon_row <- raw_uncon_rows[ref_idx, , drop = FALSE]
    conditional_low_rows <- raw_uncon_rows[!ref_idx, , drop = FALSE]
    conditional_low_rows$.traj_group <- conditional_low_label
    message(
      "Using `", uncon_traj, "` at moderator level `", reference_mod_level,
      "` as the single true reference. Other `", uncon_traj,
      "` rows were relabeled as `", conditional_low_label, "` for plotting."
    )
  }
  
  df_non_uncon <- df[df$.traj_group != uncon_traj, , drop = FALSE]
  df_comp <- rbind(df_non_uncon, conditional_low_rows)
  
  if (!is.null(mod_name)) {
    available_mod <- unique(df_comp$.moderator)
    missing_mod <- setdiff(mod_traj, available_mod)
    mod_traj_use <- intersect(mod_traj, available_mod)
    if (length(missing_mod) > 0) {
      message("Dropping missing `mod_traj` label(s): ", paste(missing_mod, collapse = ", "))
    }
    if (length(mod_traj_use) == 0) stop("None of the requested `mod_traj` labels were found in the comparison rows.")
    df_comp <- df_comp[df_comp$.moderator %in% mod_traj_use, , drop = FALSE]
    df_comp$.moderator <- factor(df_comp$.moderator, levels = mod_traj_use)
    moderator_levels <- mod_traj_use
  } else {
    df_comp$.moderator <- factor("Overall", levels = "Overall")
    moderator_levels <- "Overall"
  }
  
  if (!is.null(comp_traj)) {
    available_comp <- unique(as.character(df_comp$.traj_group))
    if (conditional_low_label %in% available_comp && !(conditional_low_label %in% comp_traj)) {
      comp_traj <- c(conditional_low_label, comp_traj)
      message("Added `", conditional_low_label, "` to `comp_traj` because moderated low-adversity rows were present.")
    }
    missing_comp <- setdiff(comp_traj, available_comp)
    comp_traj_use <- intersect(comp_traj, available_comp)
    if (length(missing_comp) > 0) {
      message("Dropping missing `comp_traj` label(s): ", paste(missing_comp, collapse = ", "))
    }
    if (length(comp_traj_use) == 0) {
      message("None of the requested `comp_traj` labels were found. Plotting only `uncon_traj`.")
      comp_traj <- NULL
    } else {
      comp_traj <- comp_traj_use
    }
  } else {
    available_comp <- unique(as.character(df_comp$.traj_group))
    if (length(available_comp) > 0) comp_traj <- available_comp
  }
  
  desired_traj_order <- if (is.null(comp_traj)) uncon_traj else c(uncon_traj, comp_traj)
  desired_traj_order <- unique(desired_traj_order)
  
  uncon_row <- uncon_row[, c(".traj_group", ".traj_int", ".traj_slope", ".sig"), drop = FALSE]
  uncon_row$.traj_group <- factor(uncon_row$.traj_group, levels = desired_traj_order)
  uncon_row$.is_unconditional <- TRUE
  uncon_row$.linewidth <- uncon_linewidth
  uncon_row$.sig_symbol <- .sig_to_symbol(uncon_row$.sig)
  
  if (!is.null(mod_name)) {
    uncon_plot_df <- do.call(rbind, lapply(moderator_levels, function(mod) {
      data.frame(
        .moderator = factor(mod, levels = moderator_levels),
        .traj_group = uncon_row$.traj_group,
        .traj_int = uncon_row$.traj_int,
        .traj_slope = uncon_row$.traj_slope,
        .sig = uncon_row$.sig,
        .sig_symbol = uncon_row$.sig_symbol,
        .is_unconditional = TRUE,
        .linewidth = uncon_linewidth,
        stringsAsFactors = FALSE
      )
    }))
  } else {
    uncon_plot_df <- data.frame(
      .moderator = factor("Overall", levels = moderator_levels),
      .traj_group = uncon_row$.traj_group,
      .traj_int = uncon_row$.traj_int,
      .traj_slope = uncon_row$.traj_slope,
      .sig = uncon_row$.sig,
      .sig_symbol = uncon_row$.sig_symbol,
      .is_unconditional = TRUE,
      .linewidth = uncon_linewidth,
      stringsAsFactors = FALSE
    )
  }
  
  if (is.null(comp_traj)) {
    comp_plot_df <- df_comp[FALSE, c(".moderator", ".traj_group", ".traj_int", ".traj_slope", ".sig"), drop = FALSE]
    comp_plot_df$.sig_symbol <- character(0)
    comp_plot_df$.is_unconditional <- logical(0)
    comp_plot_df$.linewidth <- numeric(0)
  } else {
    comp_plot_df <- df_comp[
      df_comp$.traj_group %in% comp_traj,
      c(".moderator", ".traj_group", ".traj_int", ".traj_slope", ".sig"),
      drop = FALSE
    ]
    comp_plot_df$.traj_group <- factor(as.character(comp_plot_df$.traj_group), levels = desired_traj_order)
    comp_plot_df$.sig_symbol <- .sig_to_symbol(comp_plot_df$.sig)
    comp_plot_df$.is_unconditional <- FALSE
    comp_plot_df$.linewidth <- comp_linewidth
  }
  
  plot_df <- rbind(
    uncon_plot_df[, c(".moderator", ".traj_group", ".traj_int", ".traj_slope", ".sig", ".sig_symbol", ".is_unconditional", ".linewidth")],
    comp_plot_df[, c(".moderator", ".traj_group", ".traj_int", ".traj_slope", ".sig", ".sig_symbol", ".is_unconditional", ".linewidth")]
  )
  plot_df$.moderator <- factor(as.character(plot_df$.moderator), levels = moderator_levels)
  plot_df$.traj_group <- factor(as.character(plot_df$.traj_group), levels = desired_traj_order)
  
  traj_list <- vector("list", nrow(plot_df))
  for (i in seq_len(nrow(plot_df))) {
    traj_list[[i]] <- data.frame(
      .moderator = factor(as.character(plot_df$.moderator[i]), levels = moderator_levels),
      .traj_group = factor(as.character(plot_df$.traj_group[i]), levels = desired_traj_order),
      .sig = plot_df$.sig[i],
      .sig_symbol = plot_df$.sig_symbol[i],
      .is_unconditional = plot_df$.is_unconditional[i],
      .linewidth = plot_df$.linewidth[i],
      time = time_values,
      y = plot_df$.traj_int[i] + plot_df$.traj_slope[i] * time_values,
      stringsAsFactors = FALSE
    )
  }
  traj_df <- do.call(rbind, traj_list)
  traj_df$.moderator <- factor(as.character(traj_df$.moderator), levels = moderator_levels)
  traj_df$.traj_group <- factor(as.character(traj_df$.traj_group), levels = desired_traj_order)
  
  x_limits_all <- range(traj_df$time, na.rm = TRUE)
  y_limits_all <- range(traj_df$y, na.rm = TRUE)
  y_pad <- 0.06 * diff(y_limits_all)
  if (is.na(y_pad) || y_pad == 0) y_pad <- 0.10
  y_limits_all <- c(y_limits_all[1] - y_pad, y_limits_all[2] + y_pad)
  
  grid_band_info <- .make_y_grid_and_bands(
    y_limits = y_limits_all,
    background_band_n = background_band_n,
    background_band_low = background_band_low,
    background_band_high = background_band_high,
    background_band_alpha = background_band_alpha
  )
  background_band_df <- grid_band_info$band_df
  foreground_grid_y <- grid_band_info$y_grid_breaks
  foreground_grid_x <- time_values
  
  default_palette_named <- c(
    Reference = "#0072B2",
    `Conditional Low` = "#56B4E9",
    Vertex = "#A78BFA",
    Threshold = "#2A9D8F",
    Toxic = "#C65A5A",
    PTG = "forestgreen"
  )
  fallback_palette <- c("#0072B2", "#56B4E9", "#A78BFA", "#2A9D8F", "#C65A5A",
                        "forestgreen", "#8D7B68", "#C2A14A", "#7F8C8D", "#E9C46A")
  
  if (is.null(comp_traj)) {
    resolved_comp_colors <- NULL
  } else if (!is.null(comp_colors)) {
    resolved_comp_colors <- .resolve_named_colors(comp_traj, comp_colors, "comp_colors")
  } else {
    default_matches <- comp_traj[tolower(comp_traj) %in% tolower(names(default_palette_named))]
    non_default_matches <- setdiff(comp_traj, default_matches)
    resolved_comp_colors <- c()
    if (length(default_matches) > 0) {
      resolved_comp_colors <- c(
        resolved_comp_colors,
        .resolve_named_colors(default_matches, default_palette_named, "default_palette_named")
      )
    }
    if (length(non_default_matches) > 0) {
      fallback_available <- fallback_palette[!tolower(fallback_palette) %in% tolower(unname(resolved_comp_colors))]
      if (length(non_default_matches) > length(fallback_available)) fallback_available <- fallback_palette
      resolved_comp_colors <- c(
        resolved_comp_colors,
        stats::setNames(fallback_available[seq_along(non_default_matches)], non_default_matches)
      )
    }
    resolved_comp_colors <- resolved_comp_colors[comp_traj]
  }
  
  uncon_color <- if (tolower(uncon_traj) %in% tolower(names(default_palette_named))) {
    .resolve_named_colors(uncon_traj, default_palette_named, "default_palette_named")[[uncon_traj]]
  } else {
    "black"
  }
  all_traj_colors <- c()
  all_traj_colors[uncon_traj] <- uncon_color
  if (!is.null(resolved_comp_colors)) all_traj_colors <- c(all_traj_colors, resolved_comp_colors)
  
  uncon_plot_df2 <- traj_df[traj_df$.is_unconditional, , drop = FALSE]
  comp_plot_df2 <- traj_df[!traj_df$.is_unconditional, , drop = FALSE]
  
  .base_plot_layers <- function(p_obj) {
    if (show_background_bands) {
      p_obj <- p_obj + ggplot2::geom_rect(
        data = background_band_df,
        ggplot2::aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax, fill = .band_fill),
        inherit.aes = FALSE,
        color = NA,
        show.legend = FALSE
      )
    }
    if (show_gridlines) {
      p_obj <- p_obj +
        ggplot2::geom_hline(
          yintercept = foreground_grid_y,
          color = gridline_color,
          linewidth = gridline_linewidth,
          linetype = gridline_linetype,
          alpha = 0.85
        ) +
        ggplot2::geom_vline(
          xintercept = foreground_grid_x,
          color = gridline_color,
          linewidth = gridline_linewidth,
          linetype = gridline_linetype,
          alpha = 0.85
        )
    }
    p_obj
  }
  
  p <- ggplot2::ggplot()
  p <- .base_plot_layers(p)
  
  if (!is.null(comp_traj) && nrow(comp_plot_df2) > 0) {
    p <- p + ggplot2::geom_line(
      data = comp_plot_df2,
      ggplot2::aes(x = time, y = y, color = .traj_group,
                   group = interaction(.traj_group, .moderator), linewidth = .linewidth),
      lineend = "round"
    )
    if (show_points) {
      p <- p + ggplot2::geom_point(
        data = comp_plot_df2,
        ggplot2::aes(x = time, y = y, color = .traj_group,
                     group = interaction(.traj_group, .moderator)),
        size = point_size
      )
    }
  }
  
  p <- p + ggplot2::geom_line(
    data = uncon_plot_df2,
    ggplot2::aes(x = time, y = y, color = .traj_group, linetype = .traj_group,
                 group = interaction(.traj_group, .moderator), linewidth = .linewidth),
    lineend = "round"
  )
  if (show_points) {
    p <- p + ggplot2::geom_point(
      data = uncon_plot_df2,
      ggplot2::aes(x = time, y = y, color = .traj_group, shape = .traj_group,
                   group = interaction(.traj_group, .moderator)),
      size = point_size
    )
  }
  
  if (mark_sig) {
    message(.sig_key)
    sig_df <- plot_df[plot_df$.sig_symbol != "",
                      c(".moderator", ".traj_group", ".traj_int", ".traj_slope", ".sig_symbol"), drop = FALSE]
    if (nrow(sig_df) > 0) {
      sig_nudge_x_use <- if (is.null(sig_nudge_x)) 0.02 * diff(range(time_values, na.rm = TRUE)) else sig_nudge_x
      sig_nudge_y_use <- if (is.null(sig_nudge_y)) 0.02 * diff(range(traj_df$y, na.rm = TRUE)) else sig_nudge_y
      sig_df$x <- max(time_values, na.rm = TRUE) + sig_nudge_x_use
      sig_df$y <- sig_df$.traj_int + sig_df$.traj_slope * max(time_values, na.rm = TRUE) + sig_nudge_y_use
      sig_df$.moderator <- factor(as.character(sig_df$.moderator), levels = moderator_levels)
      p <- p + ggplot2::geom_text(
        data = sig_df,
        ggplot2::aes(x = x, y = y, label = .sig_symbol),
        inherit.aes = FALSE,
        color = sig_text_color,
        show.legend = FALSE,
        size = sig_text_size,
        fontface = "bold"
      )
    }
  }
  
  p <- p + ggplot2::scale_linewidth_identity()
  if (show_background_bands) p <- p + ggplot2::scale_fill_identity(guide = "none")
  p <- p + ggplot2::scale_color_manual(values = all_traj_colors, breaks = desired_traj_order,
                                       drop = FALSE, name = color_label)
  plot_caption <- if (mark_sig && add_sig_leg) .sig_key else NULL
  p <- p +
    ggplot2::scale_linetype_manual(values = stats::setNames("solid", uncon_traj),
                                   breaks = uncon_traj, labels = uncon_traj, name = NULL)
  if (show_points) {
    p <- p + ggplot2::scale_shape_manual(values = stats::setNames(16, uncon_traj),
                                         breaks = uncon_traj, labels = uncon_traj, name = NULL)
  }
  
  p <- p +
    ggplot2::labs(x = x_label, y = y_label, title = title, caption = plot_caption) +
    ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = title_size),
      plot.caption = ggplot2::element_text(hjust = 0, size = axis_text_size),
      legend.position = "bottom",
      legend.box = "vertical",
      legend.title = ggplot2::element_text(face = "bold", size = legend_title_size),
      legend.text = ggplot2::element_text(size = legend_text_size),
      legend.key = ggplot2::element_rect(fill = bg_fill, color = NA),
      legend.background = ggplot2::element_rect(fill = bg_fill, color = NA),
      legend.box.background = ggplot2::element_rect(fill = bg_fill, color = NA),
      axis.title = ggplot2::element_text(face = "bold", size = axis_title_size),
      axis.text = ggplot2::element_text(size = axis_text_size),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold", size = strip_text_size),
      plot.background = ggplot2::element_rect(fill = bg_fill, color = NA),
      panel.background = ggplot2::element_rect(fill = bg_fill, color = NA),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank()
    )
  
  if (!is.null(mod_name) && facet_integrated_plot) {
    p <- p + ggplot2::facet_wrap(~ .moderator, nrow = 1)
  }
  
  legend_linewidths <- c(uncon_linewidth, rep(comp_linewidth, length(desired_traj_order) - 1))
  if (show_points) {
    p <- p + ggplot2::guides(
      color = ggplot2::guide_legend(order = 1, override.aes = list(linewidth = legend_linewidths, shape = 16)),
      linetype = "none",
      shape = "none"
    )
  } else {
    p <- p + ggplot2::guides(
      color = ggplot2::guide_legend(order = 1, override.aes = list(linewidth = legend_linewidths)),
      linetype = "none"
    )
  }
  
  p <- p +
    ggplot2::scale_x_continuous(breaks = time_values) +
    ggplot2::scale_y_continuous(
      breaks = foreground_grid_y,
      labels = function(x) .format_axis_labels(x, digits = y_axis_digits)
    ) +
    ggplot2::coord_cartesian(xlim = x_limits_all, ylim = y_limits_all)
  
  individual_plots <- NULL
  if (make_individual_plots) {
    if (!is.null(mod_name) && split_moderated_individual_plots) {
      plot_keys <- unique(traj_df[, c(".traj_group", ".moderator"), drop = FALSE])
      plot_keys$.traj_group <- factor(as.character(plot_keys$.traj_group), levels = desired_traj_order)
      plot_keys$.moderator <- factor(as.character(plot_keys$.moderator), levels = moderator_levels)
      plot_keys <- plot_keys[order(plot_keys$.traj_group, plot_keys$.moderator), , drop = FALSE]
    } else {
      plot_keys <- unique(traj_df[, c(".traj_group"), drop = FALSE])
      plot_keys$.moderator <- NA_character_
      plot_keys$.traj_group <- factor(as.character(plot_keys$.traj_group), levels = desired_traj_order)
      plot_keys <- plot_keys[order(plot_keys$.traj_group), , drop = FALSE]
    }
    
    individual_plots <- lapply(seq_len(nrow(plot_keys)), function(i) {
      this_traj <- as.character(plot_keys$.traj_group[i])
      this_mod <- as.character(plot_keys$.moderator[i])
      if (!is.null(mod_name) && split_moderated_individual_plots) {
        this_df <- traj_df[
          as.character(traj_df$.traj_group) == this_traj & as.character(traj_df$.moderator) == this_mod,
          , drop = FALSE
        ]
      } else {
        this_df <- traj_df[as.character(traj_df$.traj_group) == this_traj, , drop = FALSE]
      }
      if (nrow(this_df) == 0) return(NULL)
      
      this_is_uncon <- unique(this_df$.is_unconditional)
      if (length(this_is_uncon) != 1) {
        stop("Unexpected plotting issue: trajectory type is not unique within ", this_traj)
      }
      this_color <- if (this_traj %in% names(all_traj_colors)) all_traj_colors[[this_traj]] else "grey40"
      this_linewidth <- if (isTRUE(this_is_uncon)) uncon_linewidth else comp_linewidth
      
      p_single <- ggplot2::ggplot()
      p_single <- .base_plot_layers(p_single)
      p_single <- p_single +
        ggplot2::geom_line(
          data = this_df,
          ggplot2::aes(x = time, y = y, group = interaction(.traj_group, .moderator)),
          color = this_color,
          linewidth = this_linewidth,
          lineend = "round"
        )
      if (show_points) {
        p_single <- p_single + ggplot2::geom_point(
          data = this_df,
          ggplot2::aes(x = time, y = y, group = interaction(.traj_group, .moderator)),
          color = this_color,
          size = point_size
        )
      }
      
      if (mark_sig) {
        if (!is.null(mod_name) && split_moderated_individual_plots) {
          sig_single_df <- plot_df[
            as.character(plot_df$.traj_group) == this_traj &
              as.character(plot_df$.moderator) == this_mod &
              plot_df$.sig_symbol != "",
            c(".moderator", ".traj_group", ".traj_int", ".traj_slope", ".sig_symbol"),
            drop = FALSE
          ]
        } else {
          sig_single_df <- plot_df[
            as.character(plot_df$.traj_group) == this_traj & plot_df$.sig_symbol != "",
            c(".moderator", ".traj_group", ".traj_int", ".traj_slope", ".sig_symbol"),
            drop = FALSE
          ]
        }
        if (nrow(sig_single_df) > 0) {
          sig_nudge_x_use <- if (is.null(sig_nudge_x)) 0.02 * diff(range(time_values, na.rm = TRUE)) else sig_nudge_x
          sig_nudge_y_use <- if (is.null(sig_nudge_y)) 0.02 * diff(range(traj_df$y, na.rm = TRUE)) else sig_nudge_y
          sig_single_df$x <- max(time_values, na.rm = TRUE) + sig_nudge_x_use
          sig_single_df$y <- sig_single_df$.traj_int +
            sig_single_df$.traj_slope * max(time_values, na.rm = TRUE) + sig_nudge_y_use
          sig_single_df$.moderator <- factor(as.character(sig_single_df$.moderator), levels = moderator_levels)
          p_single <- p_single + ggplot2::geom_text(
            data = sig_single_df,
            ggplot2::aes(x = x, y = y, label = .sig_symbol),
            inherit.aes = FALSE,
            color = sig_text_color,
            show.legend = FALSE,
            size = sig_text_size,
            fontface = "bold"
          )
        }
      }
      
      if (show_background_bands) p_single <- p_single + ggplot2::scale_fill_identity(guide = "none")
      if (!is.null(mod_name) && !split_moderated_individual_plots) {
        p_single <- p_single + ggplot2::facet_wrap(~ .moderator, nrow = 1)
      }
      
      single_title <- if (!is.null(mod_name) && split_moderated_individual_plots) {
        paste0(title, ": ", this_traj, " at ", this_mod)
      } else {
        paste0(title, ": ", this_traj)
      }
      
      p_single <- p_single +
        ggplot2::scale_x_continuous(breaks = time_values) +
        ggplot2::scale_y_continuous(
          breaks = foreground_grid_y,
          labels = function(x) .format_axis_labels(x, digits = y_axis_digits)
        ) +
        ggplot2::coord_cartesian(xlim = x_limits_all, ylim = y_limits_all) +
        ggplot2::labs(x = x_label, y = y_label, title = single_title, caption = NULL) +
        ggplot2::theme_classic(base_size = base_size) +
        ggplot2::theme(
          plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = title_size),
          legend.position = "none",
          axis.title = ggplot2::element_text(face = "bold", size = axis_title_size),
          axis.text = ggplot2::element_text(size = axis_text_size),
          strip.background = ggplot2::element_blank(),
          strip.text = ggplot2::element_text(face = "bold", size = strip_text_size),
          plot.background = ggplot2::element_rect(fill = bg_fill, color = NA),
          panel.background = ggplot2::element_rect(fill = bg_fill, color = NA),
          panel.grid.major = ggplot2::element_blank(),
          panel.grid.minor = ggplot2::element_blank()
        )
      p_single
    })
    
    if (!is.null(mod_name) && split_moderated_individual_plots) {
      plot_names <- paste0(
        .safe_plot_name(as.character(plot_keys$.traj_group)),
        "_",
        .safe_plot_name(as.character(plot_keys$.moderator))
      )
    } else {
      plot_names <- .safe_plot_name(as.character(plot_keys$.traj_group))
    }
    names(individual_plots) <- plot_names
    individual_plots <- individual_plots[!vapply(individual_plots, is.null, logical(1))]
  }
  
  if (print_integrated_plot) print(p)
  if (make_individual_plots && print_individual_plots) {
    for (nm in names(individual_plots)) print(individual_plots[[nm]])
  }
  
  out <- list(integrated_plot = p, individual_plots = individual_plots)
  if (return_plot_data) {
    out$plot_data <- traj_df
    out$plot_table <- plot_df
    out$background_bands <- background_band_df
    out$foreground_grid_y <- foreground_grid_y
    out$foreground_grid_x <- foreground_grid_x
    out$reference_mod_level <- reference_mod_level
  }
  return(out)
}



################################################################################
############ DENSITY HISTOGRAMS

DensityHorm <- function(
    data,
    var_name,
    horm_ver,
    horm_thresh,
    horm_ref = NULL,
    toxic_x = NULL,
    bins = 80,
    x_label = "",
    y_label = "",
    title = NULL,
    x_label_size = 12,
    y_label_size = 12,
    x_text_size = 10,
    y_text_size = 10,
    title_size = 13,
    landmark_value_size = 3.4,
    landmark_digits = 2,
    x_tick_margin_top = -2,
    low_color = "#0072B2",
    ver_color = "#A78BFA",
    thresh_color = "#2A9D8F",
    high_color = "#7F1D1D",
    outline_color = "grey20",
    outline_linewidth = 0.2,
    landmark_width_frac = 0.010,
    landmark_height_frac = 0.018,
    landmark_outline_color = "white",
    landmark_outline_linewidth = 0.25,
    landmark_text_y_nudge = 1.25,
    base_size = 12
) {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required.", call. = FALSE)
  }
  if (!requireNamespace("scales", quietly = TRUE)) {
    stop("Package 'scales' is required.", call. = FALSE)
  }
  
  if (!is.data.frame(data)) {
    stop("`data` must be a data frame.")
  }
  
  if (!is.character(var_name) || length(var_name) != 1 || is.na(var_name)) {
    stop("`var_name` must be a single character string naming a column in `data`.")
  }
  
  if (!var_name %in% names(data)) {
    stop("`var_name` was not found in `data`.")
  }
  
  x <- data[[var_name]]
  if (!is.numeric(x)) {
    stop("`var_name` must name a numeric column.")
  }
  
  x <- x[!is.na(x)]
  if (length(x) == 0) {
    stop("No non-missing values found in `var_name`.")
  }
  
  x_min <- min(x, na.rm = TRUE)
  x_max <- max(x, na.rm = TRUE)
  
  if (is.null(horm_ref)) {
    horm_ref <- x_min
  }
  
  if (is.null(toxic_x)) {
    toxic_x <- x_max
  }
  
  marker_args <- list(
    horm_ref = horm_ref,
    horm_ver = horm_ver,
    horm_thresh = horm_thresh,
    toxic_x = toxic_x
  )
  
  bad_markers <- names(marker_args)[
    !vapply(marker_args, function(z) is.numeric(z) && length(z) == 1 && !is.na(z), logical(1))
  ]
  if (length(bad_markers) > 0) {
    stop(
      "These arguments must be single non-missing numeric values: ",
      paste(bad_markers, collapse = ", ")
    )
  }
  
  numeric_single_args <- list(
    bins = bins,
    x_label_size = x_label_size,
    y_label_size = y_label_size,
    x_text_size = x_text_size,
    y_text_size = y_text_size,
    title_size = title_size,
    landmark_value_size = landmark_value_size,
    landmark_digits = landmark_digits,
    x_tick_margin_top = x_tick_margin_top,
    outline_linewidth = outline_linewidth,
    landmark_width_frac = landmark_width_frac,
    landmark_height_frac = landmark_height_frac,
    landmark_outline_linewidth = landmark_outline_linewidth,
    landmark_text_y_nudge = landmark_text_y_nudge,
    base_size = base_size
  )
  
  bad_numeric <- names(numeric_single_args)[
    !vapply(numeric_single_args, function(z) is.numeric(z) && length(z) == 1 && !is.na(z), logical(1))
  ]
  if (length(bad_numeric) > 0) {
    stop(
      "These arguments must be single non-missing numeric values: ",
      paste(bad_numeric, collapse = ", ")
    )
  }
  
  if (bins < 5) {
    stop("`bins` must be >= 5.")
  }
  
  if (horm_ref > horm_ver) {
    stop("`horm_ref` must be less than or equal to `horm_ver`.")
  }
  if (horm_ver > horm_thresh) {
    stop("`horm_ver` must be less than or equal to `horm_thresh`.")
  }
  if (horm_thresh > toxic_x) {
    stop("`horm_thresh` must be less than or equal to `toxic_x`.")
  }
  
  if (x_max > horm_thresh) {
    mid_post <- horm_thresh + 0.55 * (x_max - horm_thresh)
    anchor_vals <- c(horm_ref, horm_ver, horm_thresh, mid_post, x_max)
    anchor_cols <- c(low_color, ver_color, thresh_color, "#A12F3A", high_color)
  } else {
    anchor_vals <- c(horm_ref, horm_ver, horm_thresh, x_max)
    anchor_cols <- c(low_color, ver_color, thresh_color, high_color)
  }
  
  keep_idx <- !duplicated(anchor_vals)
  anchor_vals <- anchor_vals[keep_idx]
  anchor_cols <- anchor_cols[keep_idx]
  
  ord <- order(anchor_vals)
  anchor_vals <- anchor_vals[ord]
  anchor_cols <- anchor_cols[ord]
  
  anchor_pos <- scales::rescale(anchor_vals, from = c(x_min, x_max))
  
  plot_df <- data.frame(x = x)
  
  hist_obj <- graphics::hist(x, breaks = bins, plot = FALSE)
  y_max <- max(hist_obj$counts, na.rm = TRUE)
  
  x_range <- x_max - x_min
  if (x_range == 0) x_range <- 1
  
  diamond_half_width  <- x_range * landmark_width_frac
  diamond_half_height <- y_max * landmark_height_frac
  
  landmark_df <- data.frame(
    x = c(horm_ref, horm_ver, horm_thresh, toxic_x),
    fill = c(low_color, ver_color, thresh_color, high_color),
    group = c("Reference", "Vertex", "Threshold", "Toxic"),
    stringsAsFactors = FALSE
  )
  
  diamond_df <- do.call(
    rbind,
    lapply(seq_len(nrow(landmark_df)), function(i) {
      xi <- landmark_df$x[i]
      data.frame(
        x = c(xi - diamond_half_width, xi, xi + diamond_half_width, xi),
        y = c(0, diamond_half_height, 0, -diamond_half_height),
        fill = landmark_df$fill[i],
        group = landmark_df$group[i],
        stringsAsFactors = FALSE
      )
    })
  )
  
  label_df <- data.frame(
    x = c(horm_ver, horm_thresh),
    y = c(-diamond_half_height * landmark_text_y_nudge,
          -diamond_half_height * landmark_text_y_nudge),
    label = c(
      format(round(horm_ver, landmark_digits), nsmall = landmark_digits),
      format(round(horm_thresh, landmark_digits), nsmall = landmark_digits)
    ),
    color = c(ver_color, thresh_color),
    stringsAsFactors = FALSE
  )
  
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = x)) +
    ggplot2::geom_histogram(
      ggplot2::aes(fill = after_stat(x)),
      bins = bins,
      color = outline_color,
      linewidth = outline_linewidth
    ) +
    ggplot2::scale_fill_gradientn(
      colours = anchor_cols,
      values = anchor_pos,
      limits = c(x_min, x_max),
      guide = "none",
      oob = scales::squish
    ) +
    ggplot2::geom_hline(
      yintercept = 0,
      linewidth = 0.45,
      color = "black"
    ) +
    ggplot2::geom_polygon(
      data = diamond_df,
      ggplot2::aes(x = x, y = y, group = group),
      inherit.aes = FALSE,
      fill = diamond_df$fill,
      color = landmark_outline_color,
      linewidth = landmark_outline_linewidth
    ) +
    ggplot2::geom_text(
      data = label_df,
      ggplot2::aes(x = x, y = y, label = label),
      inherit.aes = FALSE,
      color = label_df$color,
      size = landmark_value_size,
      fontface = "plain",
      vjust = 1
    ) +
    ggplot2::labs(
      x = x_label,
      y = y_label,
      title = title
    ) +
    ggplot2::coord_cartesian(
      ylim = c(-diamond_half_height * 2.8, y_max * 1.03),
      clip = "off"
    ) +
    ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = title_size),
      axis.title.x = ggplot2::element_text(face = "bold", size = x_label_size),
      axis.title.y = ggplot2::element_text(face = "bold", size = y_label_size),
      axis.text.x = ggplot2::element_text(
        size = x_text_size,
        margin = ggplot2::margin(t = x_tick_margin_top)
      ),
      axis.text.y = ggplot2::element_text(size = y_text_size),
      axis.line.x = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(t = 8, r = 8, b = 20, l = 35)
    )
  
  return(p)
}



## school adversity x pubery 

ReApp_SCHADV_PUB <- subset(df, model == "ReApp_SCHADV_PUB")

ReApp_SCH_PUB_MOD_PLOT <- PieceOfTraj(
  data = ReApp_SCHADV_PUB,
  traj_group = "adv_level",
  traj_int = "Intercept",
  traj_slope = "Slope",
  uncon_traj = "Reference",
  comp_traj = c("Conditional low", "Vertex", "Threshold", "Toxic"),
  mod_name = "mod_level",
  mod_traj = c("Low (-1 SD)", "Mean", "High (+1 SD)"),
  time_values = c(0, 1, 2),
  axis_text_size = 18,
  mark_sig = TRUE,
  sig_name = "p_val",
  sig_nudge_y = .03,
  sig_nudge_x = -.02,
  sig_text_size = 7,
  title = "Predicted Reappraisal Trajectories by Pubertal Development",
  x_label = "Time",
  y_label = "Reappraisal",
  make_individual_plots = TRUE,
  split_moderated_individual_plots = TRUE,
  print_integrated_plot = FALSE,
  print_individual_plots = TRUE
)

# save them


## reference
ggsave(
  filename = "REAP_SCHADV_modTRAJ_Ref_Mod_Low.png",
  plot = ReApp_SCH_PUB_MOD_PLOT$individual_plots$Conditional_Low_Low_minus1_SD,
  width = 5,
  height = 4,
  dpi = 600,
  bg = "transparent"
)

ggsave(
  filename = "REAP_SCHADV_modTRAJ_Ref_Mod_mean.png",
  plot = ReApp_SCH_PUB_MOD_PLOT$individual_plots$Reference_Mean,
  width = 5,
  height = 4,
  dpi = 600,
  bg = "transparent"
)

ggsave(
  filename = "REAP_SCHADV_modTRAJ_Ref_Mod_Hi.png",
  plot = ReApp_SCH_PUB_MOD_PLOT$individual_plots$Conditional_Low_High_plus1_SD,
  width = 5,
  height = 4,
  dpi = 600,
  bg = "transparent"
)

#vertex 
ggsave(
  filename = "REAP_SCHADV_modTRAJ_Vrt_Mod_Low.png",
  plot = ReApp_SCH_PUB_MOD_PLOT$individual_plots$Vertex_Low_minus1_SD,
  width = 5,
  height = 4,
  dpi = 600,
  bg = "transparent"
)

ggsave(
  filename = "REAP_SCHADV_modTRAJ_Vrt_Mod_mean.png",
  plot = ReApp_SCH_PUB_MOD_PLOT$individual_plots$Vertex_Mean,
  width = 5,
  height = 4,
  dpi = 600,
  bg = "transparent"
)

ggsave(
  filename = "REAP_SCHADV_modTRAJ_Vrt_Mod_Hi.png",
  plot = ReApp_SCH_PUB_MOD_PLOT$individual_plots$Vertex_High_plus1_SD,
  width = 5,
  height = 4,
  dpi = 600,
  bg = "transparent"
)

#threshold  
ggsave(
  filename = "REAP_SCHADV_modTRAJ_thr_Mod_Low.png",
  plot = ReApp_SCH_PUB_MOD_PLOT$individual_plots$Threshold_Low_minus1_SD,
  width = 5,
  height = 4,
  dpi = 600,
  bg = "transparent"
)

ggsave(
  filename = "REAP_SCHADV_modTRAJ_thr_Mod_mean.png",
  plot = ReApp_SCH_PUB_MOD_PLOT$individual_plots$Threshold_Mean,
  width = 5,
  height = 4,
  dpi = 600,
  bg = "transparent"
)

ggsave(
  filename = "REAP_SCHADV_modTRAJ_thr_Mod_Hi.png",
  plot = ReApp_SCH_PUB_MOD_PLOT$individual_plots$Threshold_High_plus1_SD,
  width = 5,
  height = 4,
  dpi = 600,
  bg = "transparent"
)

#toxic  
ggsave(
  filename = "REAP_SCHADV_modTRAJ_tox_Mod_Low.png",
  plot = ReApp_SCH_PUB_MOD_PLOT$individual_plots$Toxic_Low_minus1_SD,
  width = 5,
  height = 4,
  dpi = 600,
  bg = "transparent"
)

ggsave(
  filename = "REAP_SCHADV_modTRAJ_tox_Mod_mean.png",
  plot = ReApp_SCH_PUB_MOD_PLOT$individual_plots$Toxic_Mean,
  width = 5,
  height = 4,
  dpi = 600,
  bg = "transparent"
)

ggsave(
  filename = "REAP_SCHADV_modTRAJ_tox_Mod_Hi.png",
  plot = ReApp_SCH_PUB_MOD_PLOT$individual_plots$Toxic_High_plus1_SD,
  width = 5,
  height = 4,
  dpi = 600,
  bg = "transparent"
)


############### density plot and region mapping 

# moderator appears to take effect at about .5 SD above the mean (.255)
names(df_vars)

#rescale to a 0 - 100 scale for plotting
df_vars$M_SchAdvRS <- df_vars$M_SchAdv*10

## locate how many people are: a) below .25 on the moderator (i.e., in ROS)
## and b) at less than 1.58 on school adversity (hormetic inflection)

ROS <- subset(df_vars, GMC_Ptemp113 < .25)
5347/11868
11868-5347
# n = 5347 (45%)
ROS_HORM <- subset(ROS, M_SchAdvRS < 15.8) #i.e., 1.58 x 10 to be on the same scale 
1225/11868
1225/5347
#n = 1,225 (23% ROS, 10% Overall)
ROS_STG <- subset(ROS, M_SchAdvRS < 5.60) #i.e., .556 (from stage 1) x 10  
123/11868
123/5347
# n = 123 (1% overall, 2% of ROS)
ROS_BUFF <- subset(ROS, M_SchAdvRS > 5.60 & M_SchAdvRS < 15.80)  
1102/11868
1102/5347
#n = 741 (9% overall, 21% ROS)


(1102+123)/5347
(5347-(1102+123))/5347
(5347-(1102+123))
###### density histogram


REAPP_SCH_ADV_DENS <- DensityHorm(
  data = df_vars,
  var_name = "M_SchAdvRS",
  # landmarks on the ORIGINAL adversity scale
  horm_ref = 0,
  horm_ver = 5.60,
  horm_thresh = 15.80,
  x_text_size = 15,
  y_text_size = 15,
  landmark_value_size = 4.5,
  x_tick_margin_top = -15,
  landmark_text_y_nudge = 1.5,
  bins = 80,
  #  x_label = "Neighborhood Danger",
  y_label = "Count"#,
  #  title = "Distribution of Neighborhood Danger Across Hormetic Landmarks"
)

REAPP_SCH_ADV_DENS

ggsave(
  filename = "REAPP_SCHADV_MOD_DENSITY.png",
  plot = REAPP_SCH_ADV_DENS,
  width = 12,
  height = 5,
  units = "in",
  dpi = 600,
  bg = "transparent"
)

