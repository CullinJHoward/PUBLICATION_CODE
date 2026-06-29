ScatterSplatter <- function(data,
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
                            sample_n = NULL,
                            point_alpha = 0.35,
                            point_size = 1.5,
                            x_center = FALSE,
                            y_center = FALSE,
                            x_center_method = c("wave_grand_mean", "person_mean", "pooled_grand_mean"),
                            y_center_method = c("wave_grand_mean", "person_mean", "pooled_grand_mean"),
                            horm_thresh = FALSE,
                            horm_trim = 0.10,
                            horm_boot_draws = 100,
                            horm_boot_seed = 123,
                            horm_knot_select = c("median", "mean", "mode"),
                            print_plots = FALSE,
                            save_file = NULL,
                            width = 11,
                            height = 8.5,
                            plots_per_page = 4) {
  
  data_format <- match.arg(data_format)
  pair_mode   <- match.arg(pair_mode)
  x_center_method <- match.arg(x_center_method)
  y_center_method <- match.arg(y_center_method)
  horm_knot_select <- match.arg(horm_knot_select)
  
  req_pkgs <- c("dplyr", "tidyr", "ggplot2", "ggExtra", "gridExtra", "magrittr")
  if (isTRUE(horm_thresh)) {
    req_pkgs <- c(req_pkgs, "mgcv", "sandwich", "lmtest")
  }
  
  for (pkg in req_pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(paste0(pkg, " is required but not installed. Please install it."))
    }
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  }
  
  if (isTRUE(horm_thresh)) {
    if (!requireNamespace("mgcv", quietly = TRUE)) stop("mgcv required when horm_thresh = TRUE")
    if (!requireNamespace("sandwich", quietly = TRUE)) stop("sandwich required when horm_thresh = TRUE")
    if (!requireNamespace("lmtest", quietly = TRUE)) stop("lmtest required when horm_thresh = TRUE")
  }
  
  `%>%` <- magrittr::`%>%`
  
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
  
  .rp <- function(p) {
    if (is.na(p)) return("p = NA")
    if (p < .0001) return("p < .0001")
    paste0("p = ", sub("^(-?)0\\.", "\\1.", sprintf("%.4f", p)))
  }
  
  .derive_csv_name <- function(save_file) {
    if (is.null(save_file)) return("ScatterSplatter_summary.csv")
    tools::file_path_sans_ext(save_file) |> paste0(".csv")
  }
  
  .kernel_mode <- function(x) {
    x <- x[is.finite(x)]
    if (length(x) == 0) return(NA_real_)
    if (length(unique(x)) == 1) return(unique(x)[1])
    d <- stats::density(x, na.rm = TRUE)
    d$x[which.max(d$y)]
  }
  
  .select_final_knot <- function(x, method = "median") {
    x <- x[is.finite(x)]
    if (length(x) == 0) return(NA_real_)
    if (method == "median") return(stats::median(x, na.rm = TRUE))
    if (method == "mean") return(mean(x, na.rm = TRUE))
    if (method == "mode") return(.kernel_mode(x))
    NA_real_
  }
  
  .make_hormesis_summary_row <- function(x_name, y_name, res = NULL, boot = NULL,
                                         horm_extra = NULL, failure_reason = NA_character_) {
    if (is.null(res)) {
      return(data.frame(
        x_variable = x_name,
        y_variable = y_name,
        slope1_b = NA_real_,
        slope1_z = NA_real_,
        slope1_p = NA_real_,
        slope2_b = NA_real_,
        slope2_z = NA_real_,
        slope2_p = NA_real_,
        horm_vrt = NA_real_,
        knot_y = NA_real_,
        knot_y_left = NA_real_,
        knot_y_right = NA_real_,
        knot_y_gap = NA_real_,
        glm1_intercept = NA_real_,
        glm2_intercept = NA_real_,
        glm1_high_jump = NA_real_,
        glm2_high_jump = NA_real_,
        horm_ref = NA_real_,
        horm_ref_x = NA_real_,
        horm_inf_x = NA_real_,
        knot_select_method = if (!is.null(boot)) boot$knot_select_method else NA_character_,
        knot_boot_mean = if (!is.null(boot)) boot$knot_boot_mean else NA_real_,
        knot_boot_median = if (!is.null(boot)) boot$knot_boot_median else NA_real_,
        knot_boot_mode = if (!is.null(boot)) boot$knot_boot_mode else NA_real_,
        knot_boot_sd = if (!is.null(boot)) boot$knot_boot_sd else NA_real_,
        n_boot_success = if (!is.null(boot)) boot$n_boot_success else NA_integer_,
        n_boot_fail = if (!is.null(boot)) boot$n_boot_fail else NA_integer_,
        prop_boot_success = if (!is.null(boot)) boot$prop_boot_success else NA_real_,
        failure_reason = failure_reason,
        stringsAsFactors = FALSE
      ))
    }
    
    get_extra <- function(name, default = NA_real_) {
      if (is.null(horm_extra)) return(default)
      if (!name %in% names(horm_extra)) return(default)
      val <- horm_extra[[name]]
      if (length(val) == 0) return(default)
      val[1]
    }
    
    data.frame(
      x_variable = x_name,
      y_variable = y_name,
      slope1_b = unname(res$b1),
      slope1_z = unname(res$z1),
      slope1_p = unname(res$p1),
      slope2_b = unname(res$b2),
      slope2_z = unname(res$z2),
      slope2_p = unname(res$p2),
      horm_vrt = unname(res$horm_vrt),
      knot_y = get_extra("knot_y"),
      knot_y_left = get_extra("knot_y_left"),
      knot_y_right = get_extra("knot_y_right"),
      knot_y_gap = get_extra("knot_y_gap"),
      glm1_intercept = get_extra("glm1_intercept"),
      glm2_intercept = get_extra("glm2_intercept"),
      glm1_high_jump = get_extra("glm1_high_jump"),
      glm2_high_jump = get_extra("glm2_high_jump"),
      horm_ref = get_extra("horm_ref"),
      horm_ref_x = get_extra("horm_ref_x"),
      horm_inf_x = get_extra("horm_inf_x"),
      knot_select_method = boot$knot_select_method,
      knot_boot_mean = boot$knot_boot_mean,
      knot_boot_median = boot$knot_boot_median,
      knot_boot_mode = boot$knot_boot_mode,
      knot_boot_sd = boot$knot_boot_sd,
      n_boot_success = boot$n_boot_success,
      n_boot_fail = boot$n_boot_fail,
      prop_boot_success = boot$prop_boot_success,
      failure_reason = failure_reason,
      stringsAsFactors = FALSE
    )
  }
  
  .center_static <- function(df, vars) {
    if (is.null(vars) || length(vars) == 0) return(df)
    for (v in vars) {
      if (is.numeric(df[[v]])) {
        df[[v]] <- df[[v]] - mean(df[[v]], na.rm = TRUE)
      }
    }
    df
  }
  
  .center_long <- function(df, stems, method, id_var) {
    if (is.null(stems) || length(stems) == 0) return(df)
    for (v in stems) {
      if (!v %in% names(df)) next
      if (!is.numeric(df[[v]])) next
      
      if (method == "wave_grand_mean") {
        df <- df %>%
          dplyr::group_by(.data$wave) %>%
          dplyr::mutate(!!v := .data[[v]] - mean(.data[[v]], na.rm = TRUE)) %>%
          dplyr::ungroup()
      } else if (method == "person_mean") {
        df <- df %>%
          dplyr::group_by(.data[[id_var]]) %>%
          dplyr::mutate(!!v := .data[[v]] - mean(.data[[v]], na.rm = TRUE)) %>%
          dplyr::ungroup()
      } else if (method == "pooled_grand_mean") {
        df[[v]] <- df[[v]] - mean(df[[v]], na.rm = TRUE)
      }
    }
    df
  }
  
  .reg2_internal <- function(data, f, horm_vrt) {
    vars <- all.vars(f)
    data <- data[, vars, drop = FALSE]
    data <- data[stats::complete.cases(data), , drop = FALSE]
    
    y.f <- all.vars(f)[1]
    x.f <- all.vars(f)[2]
    
    xu <- data[[x.f]]
    yu <- data[[y.f]]
    
    if (nrow(data) < 6) stop("Too few complete rows for piecewise regression.")
    if (!is.numeric(xu) || !is.numeric(yu)) stop("x_value and y_value must both be numeric.")
    if (length(unique(xu)) < 3) stop("x_value has fewer than 3 unique values.")
    if (!is.finite(horm_vrt)) stop("horm_vrt is not finite.")
    
    n_left  <- sum(xu <= horm_vrt, na.rm = TRUE)
    n_right <- sum(xu >  horm_vrt, na.rm = TRUE)
    
    if (n_left < 3 || n_right < 3) {
      stop("Too few observations on one side of the horm_vrt.")
    }
    
    data$xlow1  <- ifelse(xu <= horm_vrt, xu - horm_vrt, 0)
    data$xhigh1 <- ifelse(xu >  horm_vrt, xu - horm_vrt, 0)
    data$high1  <- ifelse(xu >  horm_vrt, 1, 0)
    
    data$xlow2  <- ifelse(xu <  horm_vrt, xu - horm_vrt, 0)
    data$xhigh2 <- ifelse(xu >= horm_vrt, xu - horm_vrt, 0)
    data$high2  <- ifelse(xu >= horm_vrt, 1, 0)
    
    mod1.f <- stats::as.formula(paste0(y.f, " ~ xlow1 + xhigh1 + high1"))
    mod2.f <- stats::as.formula(paste0(y.f, " ~ xlow2 + xhigh2 + high2"))
    
    mod1 <- stats::lm(mod1.f, data = data)
    mod2 <- stats::lm(mod2.f, data = data)
    
    rob1 <- sandwich::vcovHC(mod1, type = "HC1")
    rob2 <- sandwich::vcovHC(mod2, type = "HC1")
    
    ct1 <- lmtest::coeftest(mod1, vcov. = rob1)
    ct2 <- lmtest::coeftest(mod2, vcov. = rob2)
    
    if (!("xlow1" %in% rownames(ct1))) stop("Left-side slope term dropped from model.")
    if (!("xhigh2" %in% rownames(ct2))) stop("Right-side slope term dropped from model.")
    
    b1 <- unname(stats::coef(mod1)["xlow1"])
    z1 <- unname(ct1["xlow1", 3])
    p1 <- unname(ct1["xlow1", 4])
    
    b2 <- unname(stats::coef(mod2)["xhigh2"])
    z2 <- unname(ct2["xhigh2", 3])
    p2 <- unname(ct2["xhigh2", 4])
    
    if (any(!is.finite(c(b1, z1, p1, b2, z2, p2)))) {
      stop("One or more slope statistics were not finite.")
    }
    
    beta1 <- stats::coef(mod1)
    beta2 <- stats::coef(mod2)
    
    list(
      b1 = b1, z1 = z1, p1 = p1,
      b2 = b2, z2 = z2, p2 = p2,
      horm_vrt = horm_vrt,
      glm1_intercept = unname(beta1["(Intercept)"]),
      glm2_intercept = unname(beta2["(Intercept)"]),
      glm1_high_jump = if ("high1" %in% names(beta1)) unname(beta1["high1"]) else NA_real_,
      glm2_high_jump = if ("high2" %in% names(beta2)) unname(beta2["high2"]) else NA_real_,
      glm1 = mod1,
      glm2 = mod2
    )
  }
  
  .twolines_internal <- function(data, f, trim = 0.10) {
    vars <- all.vars(f)
    data <- data[, vars, drop = FALSE]
    data <- data[stats::complete.cases(data), , drop = FALSE]
    
    y.f <- all.vars(f)[1]
    x.f <- all.vars(f)[2]
    xu <- data[[x.f]]
    yu <- data[[y.f]]
    
    if (nrow(data) < 6) stop("Too few complete x-y rows.")
    if (!is.numeric(xu) || !is.numeric(yu)) stop("x_value and y_value must both be numeric.")
    if (length(unique(xu)) < 3) stop("x_value has fewer than 3 unique values.")
    if (stats::sd(xu, na.rm = TRUE) == 0) stop("x_value has zero variance.")
    
    try_full <- function(trim_use) {
      unique.x <- length(unique(xu))
      sx.f <- paste0("s(", x.f, ", bs = 'cr', k = min(10, ", unique.x, "))")
      gam.f <- stats::as.formula(paste0(y.f, " ~ ", sx.f))
      
      gams <- mgcv::gam(gam.f, data = data)
      yhat.smooth <- as.numeric(stats::predict(gams, type = "response"))
      if (any(!is.finite(yhat.smooth))) stop("Smooth predictions were not finite.")
      
      lmq <- stats::lm(stats::as.formula(paste0(y.f, " ~ ", x.f, " + I(", x.f, "^2)")), data = data)
      quad.coef <- stats::coef(lmq)[paste0("I(", x.f, "^2)")]
      quad.dir <- ifelse(is.na(quad.coef), 1, sign(quad.coef))
      if (quad.dir == 0) quad.dir <- 1
      
      if (trim_use > 0) {
        xlo <- stats::quantile(xu, trim_use, na.rm = TRUE)
        xhi <- stats::quantile(xu, 1 - trim_use, na.rm = TRUE)
        keep <- xu > xlo & xu < xhi
      } else {
        keep <- rep(TRUE, length(xu))
      }
      
      x.middle <- xu[keep]
      y.middle <- yhat.smooth[keep]
      if (length(x.middle) < 6) stop("Too few observations after trimming.")
      
      idx.extreme <- if (quad.dir > 0) which.min(y.middle) else which.max(y.middle)
      extreme.y <- y.middle[idx.extreme]
      
      spread <- stats::sd(y.middle, na.rm = TRUE)
      if (!is.finite(spread) || spread == 0) {
        spread <- max(diff(range(y.middle, na.rm = TRUE)) * 0.08, 1e-6)
      }
      
      close.idx <- which(abs(y.middle - extreme.y) <= spread * 0.25)
      if (length(close.idx) == 0) close.idx <- idx.extreme
      
      xflat <- sort(unique(x.middle[close.idx]))
      if (length(xflat) == 0) stop("No candidate horm_vrt region found.")
      
      horm_vrt0 <- stats::median(xflat, na.rm = TRUE)
      rmid <- .reg2_internal(data = data, f = f, horm_vrt = horm_vrt0)
      
      z1a <- abs(rmid$z1)
      z2a <- abs(rmid$z2)
      if (!is.finite(z1a)) z1a <- 1
      if (!is.finite(z2a)) z2a <- 1
      
      prob <- z2a / (z1a + z2a)
      if (!is.finite(prob) || is.na(prob)) prob <- 0.5
      prob <- min(max(prob, 0), 1)
      
      horm_vrt <- as.numeric(stats::quantile(xflat, prob, na.rm = TRUE, names = FALSE))
      .reg2_internal(data = data, f = f, horm_vrt = horm_vrt)
    }
    
    res1 <- tryCatch(try_full(trim), error = function(e) e)
    if (!inherits(res1, "error")) return(res1)
    
    res2 <- tryCatch(try_full(0), error = function(e) e)
    if (!inherits(res2, "error")) return(res2)
    
    res3 <- tryCatch(
      .reg2_internal(data = data, f = f, horm_vrt = stats::median(xu, na.rm = TRUE)),
      error = function(e) e
    )
    if (!inherits(res3, "error")) return(res3)
    
    stop(paste("Attempt 1:", res1$message, "| Attempt 2:", res2$message, "| Attempt 3:", res3$message))
  }
  
  .bootstrap_knot_distribution <- function(data, f, trim, n_draws = 100, seed = 123) {
    vars <- all.vars(f)
    data <- data[, vars, drop = FALSE]
    data <- data[stats::complete.cases(data), , drop = FALSE]
    
    if (nrow(data) < 6) {
      return(list(
        knot_values = numeric(0),
        knot_boot_mean = NA_real_,
        knot_boot_median = NA_real_,
        knot_boot_mode = NA_real_,
        knot_boot_sd = NA_real_,
        n_boot_success = 0L,
        n_boot_fail = as.integer(n_draws),
        prop_boot_success = 0
      ))
    }
    
    set.seed(seed)
    knot_boot <- rep(NA_real_, n_draws)
    
    for (b in seq_len(n_draws)) {
      idx <- sample.int(nrow(data), size = nrow(data), replace = TRUE)
      dat_b <- data[idx, , drop = FALSE]
      
      res_b <- tryCatch(
        .twolines_internal(data = dat_b, f = f, trim = trim),
        error = function(e) NULL
      )
      
      if (!is.null(res_b) && is.finite(res_b$horm_vrt)) {
        knot_boot[b] <- res_b$horm_vrt
      }
    }
    
    knot_ok <- knot_boot[is.finite(knot_boot)]
    n_boot_success <- length(knot_ok)
    n_boot_fail <- n_draws - n_boot_success
    prop_boot_success <- n_boot_success / n_draws
    
    if (n_boot_success == 0) {
      return(list(
        knot_values = knot_ok,
        knot_boot_mean = NA_real_,
        knot_boot_median = NA_real_,
        knot_boot_mode = NA_real_,
        knot_boot_sd = NA_real_,
        n_boot_success = n_boot_success,
        n_boot_fail = n_boot_fail,
        prop_boot_success = prop_boot_success
      ))
    }
    
    list(
      knot_values = knot_ok,
      knot_boot_mean = mean(knot_ok, na.rm = TRUE),
      knot_boot_median = stats::median(knot_ok, na.rm = TRUE),
      knot_boot_mode = .kernel_mode(knot_ok),
      knot_boot_sd = stats::sd(knot_ok, na.rm = TRUE),
      n_boot_success = n_boot_success,
      n_boot_fail = n_boot_fail,
      prop_boot_success = prop_boot_success
    )
  }
  
  compute_knot_plot_anchors <- function(tl_res) {
    horm_vrt <- as.numeric(tl_res$horm_vrt)[1]
    
    out_na <- list(
      knot_y = NA_real_,
      knot_y_left = NA_real_,
      knot_y_right = NA_real_,
      knot_y_gap = NA_real_,
      glm1_intercept = NA_real_,
      glm2_intercept = NA_real_,
      glm1_high_jump = NA_real_,
      glm2_high_jump = NA_real_
    )
    
    if (!is.finite(horm_vrt)) return(out_na)
    
    beta1 <- tryCatch(stats::coef(tl_res$glm1), error = function(e) NULL)
    beta2 <- tryCatch(stats::coef(tl_res$glm2), error = function(e) NULL)
    
    nd_left <- data.frame(
      xlow1  = 0,
      xhigh1 = 0,
      high1  = 0
    )
    
    nd_right <- data.frame(
      xlow2  = 0,
      xhigh2 = 0,
      high2  = 1
    )
    
    knot_y_left <- tryCatch(
      as.numeric(stats::predict(tl_res$glm1, newdata = nd_left))[1],
      error = function(e) NA_real_
    )
    
    knot_y_right <- tryCatch(
      as.numeric(stats::predict(tl_res$glm2, newdata = nd_right))[1],
      error = function(e) NA_real_
    )
    
    if (!is.finite(knot_y_left)) knot_y_left <- NA_real_
    if (!is.finite(knot_y_right)) knot_y_right <- NA_real_
    
    list(
      # Default plotting anchor. This preserves the previous continuous-line reconstruction.
      knot_y = knot_y_left,
      
      # Exact anchors used by ScatterSplatter's two separate segment predictions.
      knot_y_left = knot_y_left,
      knot_y_right = knot_y_right,
      knot_y_gap = knot_y_right - knot_y_left,
      
      # Raw model coefficients that define the anchors.
      glm1_intercept = if (!is.null(beta1) && "(Intercept)" %in% names(beta1)) unname(beta1["(Intercept)"]) else NA_real_,
      glm2_intercept = if (!is.null(beta2) && "(Intercept)" %in% names(beta2)) unname(beta2["(Intercept)"]) else NA_real_,
      glm1_high_jump = if (!is.null(beta1) && "high1" %in% names(beta1)) unname(beta1["high1"]) else NA_real_,
      glm2_high_jump = if (!is.null(beta2) && "high2" %in% names(beta2)) unname(beta2["high2"]) else NA_real_
    )
  }
  
  compute_horm_inflection <- function(tl_res, x_vals) {
    x_vals <- x_vals[is.finite(x_vals)]
    
    knot_extra <- compute_knot_plot_anchors(tl_res)
    
    out_na <- c(
      knot_extra,
      list(
        horm_ref = NA_real_,
        horm_ref_x = NA_real_,
        horm_inf_x = NA_real_
      )
    )
    
    if (length(x_vals) < 2) return(out_na)
    
    horm_ref_x <- min(x_vals, na.rm = TRUE)
    max_x <- max(x_vals, na.rm = TRUE)
    horm_vrt <- as.numeric(tl_res$horm_vrt)[1]
    
    if (!is.finite(horm_ref_x) || !is.finite(max_x) || !is.finite(horm_vrt)) {
      return(out_na)
    }
    
    nd_ref <- data.frame(
      xlow1  = ifelse(horm_ref_x <= horm_vrt, horm_ref_x - horm_vrt, 0),
      xhigh1 = ifelse(horm_ref_x >  horm_vrt, horm_ref_x - horm_vrt, 0),
      high1  = ifelse(horm_ref_x >  horm_vrt, 1, 0)
    )
    
    horm_ref <- tryCatch(
      as.numeric(stats::predict(tl_res$glm1, newdata = nd_ref))[1],
      error = function(e) NA_real_
    )
    
    if (!is.finite(horm_ref)) {
      return(c(
        knot_extra,
        list(
          horm_ref = NA_real_,
          horm_ref_x = horm_ref_x,
          horm_inf_x = NA_real_
        )
      ))
    }
    
    if (!is.finite(tl_res$p2) || as.numeric(tl_res$p2)[1] >= 0.10) {
      return(c(
        knot_extra,
        list(
          horm_ref = horm_ref,
          horm_ref_x = horm_ref_x,
          horm_inf_x = NA_real_
        )
      ))
    }
    
    beta2 <- stats::coef(tl_res$glm2)
    
    needed <- c("(Intercept)", "xhigh2", "high2")
    if (!all(needed %in% names(beta2))) {
      return(c(
        knot_extra,
        list(
          horm_ref = horm_ref,
          horm_ref_x = horm_ref_x,
          horm_inf_x = NA_real_
        )
      ))
    }
    
    b0   <- as.numeric(beta2["(Intercept)"])[1]
    b_x  <- as.numeric(beta2["xhigh2"])[1]
    b_hi <- as.numeric(beta2["high2"])[1]
    
    if (!all(is.finite(c(b0, b_x, b_hi, horm_ref)))) {
      return(c(
        knot_extra,
        list(
          horm_ref = horm_ref,
          horm_ref_x = horm_ref_x,
          horm_inf_x = NA_real_
        )
      ))
    }
    
    if (abs(b_x) < sqrt(.Machine$double.eps)) {
      return(c(
        knot_extra,
        list(
          horm_ref = horm_ref,
          horm_ref_x = horm_ref_x,
          horm_inf_x = NA_real_
        )
      ))
    }
    
    horm_inf_x <- horm_vrt + (horm_ref - b0 - b_hi) / b_x
    horm_inf_x <- as.numeric(horm_inf_x)[1]
    
    if (!is.finite(horm_inf_x)) {
      horm_inf_x <- NA_real_
    }
    
    if (isTRUE(horm_inf_x <= horm_vrt)) {
      horm_inf_x <- NA_real_
    }
    
    if (isTRUE(horm_inf_x > max_x)) {
      horm_inf_x <- NA_real_
    }
    
    c(
      knot_extra,
      list(
        horm_ref = horm_ref,
        horm_ref_x = horm_ref_x,
        horm_inf_x = horm_inf_x
      )
    )
  }
  
  .make_segment_predictions <- function(tl_res, x_vals) {
    horm_vrt <- tl_res$horm_vrt
    
    x_left  <- seq(min(x_vals, na.rm = TRUE), horm_vrt, length.out = 100)
    x_right <- seq(horm_vrt, max(x_vals, na.rm = TRUE), length.out = 100)
    
    nd_left <- data.frame(
      xlow1  = ifelse(x_left <= horm_vrt, x_left - horm_vrt, 0),
      xhigh1 = ifelse(x_left >  horm_vrt, x_left - horm_vrt, 0),
      high1  = ifelse(x_left >  horm_vrt, 1, 0)
    )
    
    nd_right <- data.frame(
      xlow2  = ifelse(x_right <  horm_vrt, x_right - horm_vrt, 0),
      xhigh2 = ifelse(x_right >= horm_vrt, x_right - horm_vrt, 0),
      high2  = ifelse(x_right >= horm_vrt, 1, 0)
    )
    
    left_line <- data.frame(
      x = x_left,
      y = as.numeric(stats::predict(tl_res$glm1, newdata = nd_left)),
      segment = "Slope 1"
    )
    
    right_line <- data.frame(
      x = x_right,
      y = as.numeric(stats::predict(tl_res$glm2, newdata = nd_right)),
      segment = "Slope 2"
    )
    
    dplyr::bind_rows(left_line, right_line)
  }
  
  .make_boot_density_strip <- function(knot_vals, x_rng, y_base, strip_height, n_pts = 200) {
    knot_vals <- knot_vals[is.finite(knot_vals)]
    if (length(knot_vals) < 2) return(NULL)
    if (length(unique(knot_vals)) < 2) return(NULL)
    
    d <- stats::density(knot_vals, na.rm = TRUE, n = n_pts)
    
    keep <- d$x >= x_rng[1] & d$x <= x_rng[2]
    d$x <- d$x[keep]
    d$y <- d$y[keep]
    
    if (length(d$x) < 2 || all(!is.finite(d$y)) || max(d$y, na.rm = TRUE) <= 0) {
      return(NULL)
    }
    
    d$y_scaled <- y_base + (d$y / max(d$y, na.rm = TRUE)) * strip_height
    data.frame(x = d$x, y = d$y_scaled)
  }
  
  .collapse_labels <- function(x) {
    if (is.null(x) || length(x) == 0) return("None")
    paste(x, collapse = ", ")
  }
  
  .print_settings_summary <- function() {
    cat("\nScatterSplatter settings summary\n")
    cat("--------------------------------\n")
    cat("ID variable: ", id_var, "\n", sep = "")
    cat("Data format: ", data_format, "\n", sep = "")
    cat("Pair mode: ", pair_mode, "\n", sep = "")
    cat("X static: ", .collapse_labels(x_static), "\n", sep = "")
    cat("X longitudinal stems: ", .collapse_labels(x_long), "\n", sep = "")
    cat("Y static: ", .collapse_labels(y_static), "\n", sep = "")
    cat("Y longitudinal stems: ", .collapse_labels(y_long), "\n", sep = "")
    cat("Subsampling: ", if (is.null(sample_n)) "No; all available observations used per panel" else paste0("Yes; maximum ", sample_n, " observations per panel"), "\n", sep = "")
    cat("X centering: ", if (x_center) paste0("Yes; static = grand mean, longitudinal = ", x_center_method) else "No", "\n", sep = "")
    cat("Y centering: ", if (y_center) paste0("Yes; static = grand mean, longitudinal = ", y_center_method) else "No", "\n", sep = "")
    cat("Hormesis two-lines estimation: ", if (horm_thresh) paste0("Yes; trim = ", horm_trim) else "No", "\n", sep = "")
    if (isTRUE(horm_thresh)) {
      cat("Bootstrap draws per pairing: ", horm_boot_draws, "\n", sep = "")
      cat("Bootstrap knot selector: ", horm_knot_select, "\n", sep = "")
    }
    cat("Plots requested: ", length(plot_list), "\n", sep = "")
    cat("Hormetic summary table created: ", if (horm_thresh) "Yes (assigned as 'hormetic_thresholds')" else "No", "\n", sep = "")
    if (length(x_long) > 0 || length(y_long) > 0) {
      cat("\nCaution: Longitudinal variables were included. These plots are exploratory only. No clustering or nesting by ID is accounted for, and any apparent effects should be interpreted as pooled population-level patterns rather than formal within-person or cluster-adjusted inference.\n", sep = "")
    }
    cat("\n")
  }
  
  if (!id_var %in% names(data)) stop("id_var not found in data")
  
  if (data_format == "wide") {
    check_vars_exist(x_static, names(data), "x_static")
    check_vars_exist(y_static, names(data), "y_static")
    check_stems_exist(x_long, names(data), "x_long")
    check_stems_exist(y_long, names(data), "y_long")
  }
  
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
      df_long <- df_rep %>% dplyr::left_join(df_static, by = id_var)
    } else {
      df_long <- df_rep
    }
    
    if (isTRUE(x_center)) {
      df_long <- .center_static(df_long, x_static)
      df_long <- .center_long(df_long, x_long, x_center_method, id_var)
    }
    if (isTRUE(y_center)) {
      df_long <- .center_static(df_long, y_static)
      df_long <- .center_long(df_long, y_long, y_center_method, id_var)
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
      plot_data <- x_df %>% dplyr::inner_join(y_df, by = c(id_var, "wave"))
    } else {
      plot_data <- x_df %>%
        dplyr::inner_join(y_df, by = id_var, relationship = "many-to-many") %>%
        dplyr::mutate(wave = wave.y)
    }
  } else {
    stop("data_format = 'long' is not yet implemented.")
  }
  
  hormetic_rows <- list()
  
  make_joint_plot <- function(df, x_name, y_name) {
    df_full <- df %>%
      dplyr::filter(
        x_var == x_name,
        y_var == y_name,
        !is.na(x_value),
        !is.na(y_value)
      )
    
    if (nrow(df_full) == 0) {
      if (isTRUE(horm_thresh)) {
        hormetic_rows[[paste(x_name, y_name, sep = "___")]] <<-
          .make_hormesis_summary_row(
            x_name, y_name, NULL,
            boot = list(
              knot_select_method = horm_knot_select,
              knot_boot_mean = NA_real_,
              knot_boot_median = NA_real_,
              knot_boot_mode = NA_real_,
              knot_boot_sd = NA_real_,
              n_boot_success = NA_integer_,
              n_boot_fail = NA_integer_,
              prop_boot_success = NA_real_
            ),
            horm_extra = list(horm_ref = NA_real_, horm_ref_x = NA_real_, horm_inf_x = NA_real_),
            failure_reason = "No complete x-y rows for this pairing."
          )
      }
      return(NULL)
    }
    
    df_est <- df_full
    
    df_plot <- df_full
    if (!is.null(sample_n)) {
      if (!is.numeric(sample_n) || length(sample_n) != 1 || sample_n <= 0) {
        stop("sample_n must be NULL or a single positive number.")
      }
      if (nrow(df_plot) > sample_n) {
        set.seed(123)
        df_plot <- df_plot %>% dplyr::sample_n(sample_n)
      }
    }
    
    df_est$wave  <- as.factor(df_est$wave)
    df_plot$wave <- as.factor(df_plot$wave)
    
    if (!isTRUE(horm_thresh)) {
      p <- ggplot2::ggplot(df_plot, ggplot2::aes(x = x_value, y = y_value, color = wave)) +
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
      
      return(
        ggExtra::ggMarginal(
          p,
          type = "histogram",
          margins = "both",
          groupColour = FALSE,
          groupFill = FALSE
        )
      )
    }
    
    boot_raw <- .bootstrap_knot_distribution(
      data = as.data.frame(df_est[, c("x_value", "y_value")]),
      f = y_value ~ x_value,
      trim = horm_trim,
      n_draws = horm_boot_draws,
      seed = horm_boot_seed
    )
    
    final_knot <- .select_final_knot(boot_raw$knot_values, method = horm_knot_select)
    
    if (!is.finite(final_knot)) {
      fail_reason <- "No successful bootstrap knot estimates were obtained."
      
      hormetic_rows[[paste(x_name, y_name, sep = "___")]] <<-
        .make_hormesis_summary_row(
          x_name, y_name, NULL,
          boot = c(boot_raw, list(knot_select_method = horm_knot_select)),
          horm_extra = list(horm_ref = NA_real_, horm_ref_x = NA_real_, horm_inf_x = NA_real_),
          failure_reason = fail_reason
        )
      
      p <- ggplot2::ggplot(df_plot, ggplot2::aes(x = x_value, y = y_value)) +
        ggplot2::geom_point(alpha = point_alpha, size = point_size, color = "grey40") +
        ggplot2::theme_bw() +
        ggplot2::theme(
          legend.position = "none",
          plot.title = ggplot2::element_text(hjust = 0.5)
        ) +
        ggplot2::labs(
          title = paste(y_name, "vs", x_name),
          x = x_name,
          y = y_name
        )
      
      return(
        ggExtra::ggMarginal(
          p,
          type = "histogram",
          margins = "both",
          groupColour = FALSE,
          groupFill = FALSE
        )
      )
    }
    
    tl_out <- tryCatch(
      {
        res <- .reg2_internal(
          data = as.data.frame(df_est[, c("x_value", "y_value")]),
          f = y_value ~ x_value,
          horm_vrt = final_knot
        )
        list(result = res, reason = NA_character_)
      },
      error = function(e) {
        list(result = NULL, reason = e$message)
      }
    )
    
    tl_res <- tl_out$result
    fail_reason <- tl_out$reason
    
    if (is.null(tl_res) || any(!is.finite(c(
      tl_res$b1, tl_res$z1, tl_res$p1,
      tl_res$b2, tl_res$z2, tl_res$p2, tl_res$horm_vrt
    )))) {
      if (is.null(fail_reason) || is.na(fail_reason)) {
        fail_reason <- "Final fixed-knot refit returned non-finite estimates."
      }
      
      hormetic_rows[[paste(x_name, y_name, sep = "___")]] <<-
        .make_hormesis_summary_row(
          x_name, y_name, NULL,
          boot = c(boot_raw, list(knot_select_method = horm_knot_select)),
          horm_extra = list(horm_ref = NA_real_, horm_ref_x = NA_real_, horm_inf_x = NA_real_),
          failure_reason = fail_reason
        )
      
      p <- ggplot2::ggplot(df_plot, ggplot2::aes(x = x_value, y = y_value)) +
        ggplot2::geom_point(alpha = point_alpha, size = point_size, color = "grey40") +
        ggplot2::theme_bw() +
        ggplot2::theme(
          legend.position = "none",
          plot.title = ggplot2::element_text(hjust = 0.5)
        ) +
        ggplot2::labs(
          title = paste(y_name, "vs", x_name),
          x = x_name,
          y = y_name
        )
      
      return(
        ggExtra::ggMarginal(
          p,
          type = "histogram",
          margins = "both",
          groupColour = FALSE,
          groupFill = FALSE
        )
      )
    }
    
    x_obs <- df_est$x_value[is.finite(df_est$x_value)]
    horm_extra <- compute_horm_inflection(tl_res, x_obs)
    
    hormetic_rows[[paste(x_name, y_name, sep = "___")]] <<-
      .make_hormesis_summary_row(
        x_name, y_name, tl_res,
        boot = c(boot_raw, list(knot_select_method = horm_knot_select)),
        horm_extra = horm_extra,
        failure_reason = NA_character_
      )
    
    horm_vrt <- as.numeric(tl_res$horm_vrt)[1]
    
    df_plot <- df_plot %>%
      dplyr::mutate(slope_group = ifelse(x_value <= horm_vrt, "Slope 1", "Slope 2"))
    
    line_df <- tryCatch(
      .make_segment_predictions(tl_res, df_est$x_value),
      error = function(e) NULL
    )
    
    x_rng <- range(df_est$x_value, na.rm = TRUE)
    y_rng <- range(df_est$y_value, na.rm = TRUE)
    xr <- diff(x_rng)
    yr <- diff(y_rng)
    
    if (!is.finite(xr) || xr == 0) xr <- 1
    if (!is.finite(yr) || yr == 0) yr <- 1
    
    y_pad_bottom <- 0.24 * yr
    y_pad_top    <- 0.13 * yr
    y_lower_lim  <- y_rng[1] - y_pad_bottom
    y_upper_lim  <- y_rng[2] + y_pad_top
    
    y_horm_vrt_lab <- y_lower_lim + 0.78 * y_pad_bottom
    y_slope1_lab   <- y_lower_lim + 0.46 * y_pad_bottom
    y_slope2_lab   <- y_lower_lim + 0.20 * y_pad_bottom
    
    x_left_default  <- x_rng[1] + 0.03 * xr
    x_right_default <- x_rng[1] + 0.97 * xr
    too_close_left  <- abs(x_left_default  - horm_vrt) < 0.12 * xr
    too_close_right <- abs(x_right_default - horm_vrt) < 0.12 * xr
    x_left_lab  <- if (too_close_left)  x_rng[1] + 0.18 * xr else x_left_default
    x_right_lab <- if (too_close_right) x_rng[1] + 0.82 * xr else x_right_default
    x_horm_vrt_lab <- min(max(horm_vrt, x_rng[1] + 0.10 * xr), x_rng[2] - 0.10 * xr)
    
    point_cols <- c("Slope 1" = "deepskyblue2", "Slope 2" = "tomato2")
    line_cols  <- c("Slope 1" = "blue4", "Slope 2" = "red4")
    combined_cols <- c(
      "Slope 1" = "deepskyblue2",
      "Slope 2" = "tomato2",
      "Slope 1_line" = "blue4",
      "Slope 2_line" = "red4"
    )
    
    density_y_base <- y_rng[2] + 0.22 * y_pad_top
    density_height <- 0.55 * y_pad_top
    density_df <- .make_boot_density_strip(
      knot_vals = boot_raw$knot_values,
      x_rng = x_rng,
      y_base = density_y_base,
      strip_height = density_height
    )
    
    p <- ggplot2::ggplot(df_plot, ggplot2::aes(x = x_value, y = y_value)) +
      ggplot2::geom_point(
        ggplot2::aes(color = slope_group),
        alpha = point_alpha,
        size = point_size
      ) +
      ggplot2::annotate(
        "segment",
        x = horm_vrt, xend = horm_vrt,
        y = y_rng[1], yend = y_rng[2],
        color = "black",
        linewidth = 1.2
      ) +
      ggplot2::coord_cartesian(
        xlim = x_rng,
        ylim = c(y_lower_lim, y_upper_lim)
      ) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        legend.position = "none",
        plot.title = ggplot2::element_text(hjust = 0.5),
        axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 16))
      ) +
      ggplot2::labs(
        title = paste(y_name, "vs", x_name),
        x = x_name,
        y = y_name
      ) +
      ggplot2::annotate(
        "text",
        x = x_horm_vrt_lab,
        y = y_horm_vrt_lab,
        label = paste0("Horm_vrt = ", round(horm_vrt, 2)),
        color = "black",
        fontface = "bold",
        vjust = 0.5,
        hjust = 0.5,
        size = 3.2
      ) +
      ggplot2::annotate(
        "text",
        x = x_left_lab,
        y = y_slope1_lab,
        label = paste0(
          "Slope 1: b = ", round(tl_res$b1, 2),
          ", z = ", round(tl_res$z1, 2),
          ", ", .rp(tl_res$p1)
        ),
        color = line_cols[["Slope 1"]],
        hjust = 0,
        vjust = 0.5,
        size = 3.0
      ) +
      ggplot2::annotate(
        "text",
        x = x_right_lab,
        y = y_slope2_lab,
        label = paste0(
          "Slope 2: b = ", round(tl_res$b2, 2),
          ", z = ", round(tl_res$z2, 2),
          ", ", .rp(tl_res$p2)
        ),
        color = line_cols[["Slope 2"]],
        hjust = 1,
        vjust = 0.5,
        size = 3.0
      )
    
    if (!is.null(line_df)) {
      line_df$plot_col <- ifelse(line_df$segment == "Slope 1", "Slope 1_line", "Slope 2_line")
      
      p <- p +
        ggplot2::geom_line(
          data = line_df,
          ggplot2::aes(x = x, y = y, color = plot_col),
          linewidth = 1.8,
          inherit.aes = FALSE,
          alpha = 1
        )
    }
    
    p <- p +
      ggplot2::scale_color_manual(values = combined_cols, guide = "none")
    
    if (isTRUE(is.finite(horm_extra$horm_inf_x))) {
      p <- p +
        ggplot2::annotate(
          "rect",
          xmin = horm_extra$horm_inf_x,
          xmax = x_rng[2],
          ymin = y_rng[1],
          ymax = y_rng[2],
          fill = "grey70",
          alpha = .20
        )
    }
    
    if (!is.null(density_df)) {
      p <- p +
        ggplot2::geom_area(
          data = density_df,
          ggplot2::aes(x = x, y = y),
          inherit.aes = FALSE,
          fill = "red3",
          alpha = 0.35
        ) +
        ggplot2::annotate(
          "segment",
          x = min(density_df$x, na.rm = TRUE),
          xend = max(density_df$x, na.rm = TRUE),
          y = density_y_base,
          yend = density_y_base,
          color = "red3",
          linewidth = 0.6,
          alpha = 0.8
        )
    }
    
    ggExtra::ggMarginal(
      p,
      type = "histogram",
      margins = "both",
      groupColour = FALSE,
      groupFill = FALSE
    )
  }
  
  combos <- expand.grid(
    x_var = unique(c(x_static, x_long)),
    y_var = unique(c(y_static, y_long)),
    stringsAsFactors = FALSE
  )
  
  plot_list <- vector("list", nrow(combos))
  
  for (i in seq_len(nrow(combos))) {
    plot_list[[i]] <- make_joint_plot(plot_data, combos$x_var[i], combos$y_var[i])
  }
  
  plot_list <- plot_list[!sapply(plot_list, is.null)]
  
  hormetic_thresholds <- NULL
  if (isTRUE(horm_thresh)) {
    hormetic_thresholds <- dplyr::bind_rows(hormetic_rows)
    assign("hormetic_thresholds", hormetic_thresholds, envir = parent.frame())
  }
  
  print_pages <- function(plots) {
    n <- length(plots)
    idx <- seq(1, n, by = plots_per_page)
    for (i in idx) {
      page_plots <- plots[i:min(i + plots_per_page - 1, n)]
      gridExtra::grid.arrange(grobs = page_plots, ncol = 2, nrow = 2)
    }
  }
  
  if (!is.null(save_file)) {
    grDevices::pdf(save_file, width = width, height = height)
    print_pages(plot_list)
    grDevices::dev.off()
  } else if (print_plots) {
    print_pages(plot_list)
  }
  
  if (isTRUE(horm_thresh) && !is.null(hormetic_thresholds)) {
    csv_name <- .derive_csv_name(save_file)
    utils::write.csv(hormetic_thresholds, file = csv_name, row.names = FALSE)
  }
  
  settings_summary <- list(
    id_var = id_var,
    data_format = data_format,
    pair_mode = pair_mode,
    x_static = x_static,
    x_long = x_long,
    y_static = y_static,
    y_long = y_long,
    sample_n = sample_n,
    x_center = x_center,
    y_center = y_center,
    x_center_method = if (x_center) x_center_method else NA_character_,
    y_center_method = if (y_center) y_center_method else NA_character_,
    horm_thresh = horm_thresh,
    horm_trim = if (horm_thresh) horm_trim else NA_real_,
    horm_boot_draws = if (horm_thresh) horm_boot_draws else NA_integer_,
    horm_knot_select = if (horm_thresh) horm_knot_select else NA_character_,
    n_plots = length(plot_list),
    hormetic_table_created = isTRUE(horm_thresh),
    summary_csv = if (isTRUE(horm_thresh)) .derive_csv_name(save_file) else NA_character_
  )
  
  .print_settings_summary()
  
  invisible(list(
    data = plot_data,
    plots = plot_list,
    hormetic_thresholds = hormetic_thresholds,
    settings_summary = settings_summary
  ))
}
#applied

setwd("/home/cjh37695/ABCD_PROJECTS/HORMESIS_METHODS/")
df <- read.csv("ABCD_HORM_METH_SEM_STRUCTRUAL_4.24.26.csv")

names(df)


res <- ScatterSplatter(
  data = df,
  id_var = "subID",
  data_format = "wide",
  x_static = c("WP_fHSES"),
  x_long = c("WP_mSchAdv", "WP_mNBdang", "WP_FAMCONY"),
  y_long = c("mEmoSup", "mReApp"),
  pair_mode = "all",
  sample_n = 2500,
  x_center = FALSE,
  y_center = TRUE,
  y_center_method = "wave_grand_mean",
  horm_thresh = TRUE,
  horm_trim = 0.05,
  horm_boot_draws = 100,
  horm_knot_select = "median",
  save_file = "HORMETIC_INFLECTIONS_BS100_WAVE_GMC.pdf"
)

write.csv(hormetic_thresholds, "HORMETIC_INFLECTIONS_BS100_WAVE_GMC.csv", row.names = F)


res <- ScatterSplatter(
  data = df,
  id_var = "subID",
  data_format = "wide",
  x_static = c("WP_fHSES"),
  x_long = c("WP_mSchAdv", "WP_mNBdang", "WP_FAMCONY"),
  y_long = c("mEmoSup", "mReApp"),
  pair_mode = "same_wave",
  sample_n = 2500,
  x_center = FALSE,
  y_center = TRUE,
  y_center_method = "wave_grand_mean",
  horm_thresh = TRUE,
  horm_trim = 0.05,
  horm_boot_draws = 100,
  horm_knot_select = "median",
  save_file = "HORMETIC_INFLECTIONS_BS100_WAVE_GMC_TEST.pdf"
)

write.csv(hormetic_thresholds, "HORMETIC_INFLECTIONS_BS100_WAVE_GMC_TEST.csv", row.names = F)
