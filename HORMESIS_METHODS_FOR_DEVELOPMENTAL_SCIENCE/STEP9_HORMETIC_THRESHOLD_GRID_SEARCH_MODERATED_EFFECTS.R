################################################################################
# Moderated Hormetic Reference-Return Threshold Grid Search
# Family Conflict x Parental Monitoring
# Mplus LGC spline output -> R grid search
#
# Purpose:
#   Identify whether and where the post-knot family-conflict trajectory returns
#   to, or crosses beyond, the fixed low-adversity reference trajectory at each
#   selected level of parental monitoring.
#
# Probed moderation:
#   Adversity: MS1_FCONemosup / MS2_FCONemosup
#   Moderator: mLMC_mPmon
#   Interactions: S1int = MS1_FCONemosup * mLMC_mPmon
#                 S2int = MS2_FCONemosup * mLMC_mPmon
#
# Non-probed moderator:
#   mLMC_mPNH and S3int/S4int are retained in the Mplus model, but held at
#   zero/mean for this grid search.
#
# Fixed reference:
#   The reference is the low-family-conflict trajectory with all predictors held
#   at their mean/zero. It does NOT vary across parental-monitoring levels.
#
# Primary threshold rule:
#   Area-under-difference / average maladaptive difference.
#   The threshold is the first post-knot adversity value where the expected
#   trajectory is, on average across modeled time, at or beyond the fixed
#   low-adversity reference trajectory in the maladaptive direction.
################################################################################

library(MplusAutomation)
library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(ggrepel)
library(scales)

################################################################################
# 1. User settings
################################################################################

# Path to the Mplus output file produced by:
#   supp_lineargrowth_fcon_pmon_phn_MODTHRESH_PMON.inp
out_file <- "C:\\Users\\cjh37695\\Dropbox\\HORMESIS_METHODS\\ANALYSIS\\STRUCTURAL\\MODERATED\\goals_lineargrowth_schadv_pmon_pnh_threshold_search_int2.out"

# Output prefix for saved files.
out_prefix <- "SUPR_FCONxPMON_THRESH"

# Growth-model time scores. Update if the LGC loadings are not 0, 1, 2.
time_values <- c(0, 1, 2)

# Family-conflict grid resolution.
# Smaller = more precise but larger data frame.
x_step <- 0.001

# Direction of worse functioning for the outcome.
# For emotion suppression, higher scores are usually worse; verify for your scale.
# Use:
#   "higher_is_worse" for symptoms, dysregulation, suppression, impairment, risk.
#   "lower_is_worse"  for competence, adaptive functioning, wellbeing, integrity.
outcome_direction <- "higher_is_worse"

# Optional tolerance around exact equality.
# Keep at 0 for the pure mathematical reference-return threshold.
threshold_margin <- 0

# Plot settings.
n_plot_trajectories <- 10

# Save outputs?
save_outputs <- FALSE

################################################################################
# 2. Helper functions
################################################################################

extract_mplus_param <- function(param_df, param_name) {
  
  possible_name_cols <- c("param", "parameter", "label")
  name_col <- possible_name_cols[possible_name_cols %in% names(param_df)][1]
  
  if (is.na(name_col)) {
    stop(
      "Could not identify a parameter-name column in the MplusAutomation output. ",
      "Expected one of: ", paste(possible_name_cols, collapse = ", ")
    )
  }
  
  possible_est_cols <- c("est", "estimate", "est_std")
  est_col <- possible_est_cols[possible_est_cols %in% names(param_df)][1]
  
  if (is.na(est_col)) {
    stop(
      "Could not identify an estimate column in the MplusAutomation output. ",
      "Expected one of: ", paste(possible_est_cols, collapse = ", ")
    )
  }
  
  param_df_clean <- param_df %>%
    mutate(
      .param_upper = toupper(as.character(.data[[name_col]]))
    )
  
  out <- param_df_clean %>%
    filter(.data$.param_upper == toupper(param_name))
  
  if (nrow(out) == 0) {
    stop(paste0("Could not find parameter: ", param_name))
  }
  
  if (nrow(out) > 1) {
    warning(paste0(
      "More than one row found for parameter: ", param_name,
      ". Using the first row."
    ))
  }
  
  as.numeric(out[[est_col]][1])
}

make_spline_terms <- function(x) {
  tibble(
    x  = x,
    S1 = ifelse(x < 0, x, 0),
    S2 = ifelse(x > 0, x, 0)
  )
}

compute_moderated_trajectory <- function(x, m, time_values, pars) {
  
  spline_df <- make_spline_terms(x)
  
  traj_df <- tidyr::crossing(
    spline_df,
    tibble(m_value = m)
  ) %>%
    mutate(
      I_xm = pars$BMI +
        pars$BIX1*S1 + pars$BIX2*S2 +
        pars$BIM*m_value +
        pars$BIX1M*S1*m_value + pars$BIX2M*S2*m_value,
      
      S_xm = pars$BMS +
        pars$BSX1*S1 + pars$BSX2*S2 +
        pars$BSM*m_value +
        pars$BSX1M*S1*m_value + pars$BSX2M*S2*m_value
    ) %>%
    tidyr::crossing(time = time_values) %>%
    mutate(
      y_hat = I_xm + S_xm*time
    )
  
  traj_df
}

compute_fixed_reference <- function(time_values, pars) {
  tibble(
    time = time_values,
    y_ref = pars$REFI + pars$REFS*time
  )
}

add_directional_difference <- function(grid_df,
                                       outcome_direction = c("higher_is_worse",
                                                             "lower_is_worse"),
                                       threshold_margin = 0) {
  
  outcome_direction <- match.arg(outcome_direction)
  
  grid_df %>%
    mutate(
      # Raw difference from fixed low-adversity reference trajectory.
      # Positive = conditional trajectory is above reference.
      # Negative = conditional trajectory is below reference.
      diff_from_ref = y_hat - y_ref,
      
      # Directional difference where positive always means worse-than-reference.
      maladaptive_diff = case_when(
        outcome_direction == "higher_is_worse" ~ diff_from_ref,
        outcome_direction == "lower_is_worse"  ~ -diff_from_ref,
        TRUE ~ NA_real_
      ),
      
      crossed_reference = maladaptive_diff >= threshold_margin,
      
      direction_label = case_when(
        outcome_direction == "higher_is_worse" ~ "Higher values indicate worse functioning",
        outcome_direction == "lower_is_worse"  ~ "Lower values indicate worse functioning",
        TRUE ~ NA_character_
      )
    )
}

summarise_grid_by_m_x <- function(grid_df, threshold_margin = 0) {
  
  grid_df %>%
    filter(x >= 0) %>%
    group_by(m_label, m_value, x) %>%
    summarise(
      mean_maladaptive_diff = mean(maladaptive_diff, na.rm = TRUE),
      min_maladaptive_diff  = min(maladaptive_diff, na.rm = TRUE),
      max_maladaptive_diff  = max(maladaptive_diff, na.rm = TRUE),
      
      prop_time_crossed = mean(maladaptive_diff >= threshold_margin, na.rm = TRUE),
      n_time_crossed    = sum(maladaptive_diff >= threshold_margin, na.rm = TRUE),
      n_time_total      = sum(!is.na(maladaptive_diff)),
      
      crossed_any_time = any(maladaptive_diff >= threshold_margin, na.rm = TRUE),
      crossed_all_time = all(maladaptive_diff >= threshold_margin, na.rm = TRUE),
      
      .groups = "drop"
    ) %>%
    mutate(
      crossed_average = mean_maladaptive_diff >= threshold_margin
    )
}

find_reference_return_average <- function(grid_summary, threshold_margin = 0) {
  
  grid_summary %>%
    group_by(m_label, m_value) %>%
    group_modify(~ {
      
      threshold_df <- .x %>%
        filter(mean_maladaptive_diff >= threshold_margin) %>%
        arrange(x)
      
      if (nrow(threshold_df) == 0) {
        return(
          tibble(
            threshold_found = FALSE,
            threshold_x = NA_real_,
            threshold_margin = threshold_margin,
            mean_maladaptive_diff_at_threshold = NA_real_,
            prop_time_crossed_at_threshold = NA_real_,
            n_time_crossed_at_threshold = NA_real_,
            n_time_total = NA_real_,
            threshold_rule = "average_maladaptive_difference"
          )
        )
      }
      
      tibble(
        threshold_found = TRUE,
        threshold_x = threshold_df$x[1],
        threshold_margin = threshold_margin,
        mean_maladaptive_diff_at_threshold = threshold_df$mean_maladaptive_diff[1],
        prop_time_crossed_at_threshold = threshold_df$prop_time_crossed[1],
        n_time_crossed_at_threshold = threshold_df$n_time_crossed[1],
        n_time_total = threshold_df$n_time_total[1],
        threshold_rule = "average_maladaptive_difference"
      )
    }) %>%
    ungroup() %>%
    arrange(m_value)
}

find_threshold_by_time <- function(grid_traj, threshold_margin = 0) {
  
  grid_traj %>%
    filter(x >= 0) %>%
    group_by(m_label, m_value, time) %>%
    summarise(
      threshold_found = any(maladaptive_diff >= threshold_margin, na.rm = TRUE),
      threshold_x = ifelse(
        threshold_found,
        min(x[maladaptive_diff >= threshold_margin], na.rm = TRUE),
        NA_real_
      ),
      min_maladaptive_diff = min(maladaptive_diff, na.rm = TRUE),
      max_maladaptive_diff = max(maladaptive_diff, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      threshold_margin = threshold_margin
    )
}

find_alternative_thresholds <- function(grid_summary, threshold_result) {
  
  grid_summary %>%
    group_by(m_label, m_value) %>%
    group_modify(~ {
      this_threshold <- threshold_result %>%
        filter(m_label == .y$m_label, m_value == .y$m_value)
      
      tibble(
        rule = c("any_time", "majority_time", "all_time", "average_maladaptive_difference"),
        threshold_x = c(
          .x %>%
            filter(crossed_any_time) %>%
            summarise(x = ifelse(n() > 0, min(x), NA_real_)) %>%
            pull(x),
          
          .x %>%
            filter(prop_time_crossed >= 0.50) %>%
            summarise(x = ifelse(n() > 0, min(x), NA_real_)) %>%
            pull(x),
          
          .x %>%
            filter(crossed_all_time) %>%
            summarise(x = ifelse(n() > 0, min(x), NA_real_)) %>%
            pull(x),
          
          this_threshold$threshold_x[1]
        )
      )
    }) %>%
    ungroup() %>%
    arrange(m_value, rule)
}

make_pieceoftraj_df <- function(threshold_result, pars, time_values) {
  
  # Fixed reference: one row only. PieceOfTraj() can duplicate it across facets.
  ref_row <- tibble(
    moderator_label = "Overall",
    moderator_value = 0,
    trajectory_type = "Reference",
    adversity_x = pars$XLOW,
    threshold_found = TRUE,
    Intercept = pars$REFI,
    Slope = pars$REFS
  )
  
  # Vertex and toxic rows for every moderator level.
  nonthreshold_rows <- threshold_result %>%
    select(m_label, m_value, threshold_found, threshold_x) %>%
    tidyr::crossing(
      tibble(
        trajectory_type = c("Vertex", "Toxic"),
        adversity_x = c(0, pars$XHIGH)
      )
    ) %>%
    mutate(
      S1 = ifelse(adversity_x < 0, adversity_x, 0),
      S2 = ifelse(adversity_x > 0, adversity_x, 0),
      Intercept = pars$BMI +
        pars$BIX1*S1 + pars$BIX2*S2 +
        pars$BIM*m_value +
        pars$BIX1M*S1*m_value + pars$BIX2M*S2*m_value,
      Slope = pars$BMS +
        pars$BSX1*S1 + pars$BSX2*S2 +
        pars$BSM*m_value +
        pars$BSX1M*S1*m_value + pars$BSX2M*S2*m_value,
      moderator_label = m_label,
      moderator_value = m_value
    ) %>%
    select(
      moderator_label, moderator_value, trajectory_type, adversity_x,
      threshold_found, Intercept, Slope
    )
  
  # Threshold rows only where a reference-return threshold was found.
  threshold_rows <- threshold_result %>%
    filter(threshold_found) %>%
    transmute(
      moderator_label = m_label,
      moderator_value = m_value,
      trajectory_type = "Threshold",
      adversity_x = threshold_x,
      threshold_found = threshold_found
    ) %>%
    mutate(
      S1 = ifelse(adversity_x < 0, adversity_x, 0),
      S2 = ifelse(adversity_x > 0, adversity_x, 0),
      Intercept = pars$BMI +
        pars$BIX1*S1 + pars$BIX2*S2 +
        pars$BIM*moderator_value +
        pars$BIX1M*S1*moderator_value + pars$BIX2M*S2*moderator_value,
      Slope = pars$BMS +
        pars$BSX1*S1 + pars$BSX2*S2 +
        pars$BSM*moderator_value +
        pars$BSX1M*S1*moderator_value + pars$BSX2M*S2*moderator_value
    ) %>%
    select(
      moderator_label, moderator_value, trajectory_type, adversity_x,
      threshold_found, Intercept, Slope
    )
  
  bind_rows(ref_row, nonthreshold_rows, threshold_rows) %>%
    mutate(
      trajectory_type = factor(
        trajectory_type,
        levels = c("Reference", "Vertex", "Threshold", "Toxic")
      ),
      moderator_label = factor(
        moderator_label,
        levels = c("Overall", levels(threshold_result$m_label))
      )
    ) %>%
    arrange(moderator_value, trajectory_type)
}

################################################################################
# 3. Read Mplus output and extract needed parameters
################################################################################

mplus_out <- readModels(out_file, what = "parameters")
param_df <- mplus_out$parameters$unstandardized

needed_params <- c(
  # adversity range and fixed reference
  "XLOW", "XHIGH", "REFI", "REFS", "RT0", "RT1", "RT2",
  
  # moderator probe values
  "M2L", "M15L", "M1L", "M05L", "M0", "M05H", "M1H", "M15H", "M2H",
  
  # model ingredients
  "BMI", "BMS",
  "BIX1", "BIX2", "BSX1", "BSX2",
  "BIM", "BSM", "BIX1M", "BIX2M", "BSX1M", "BSX2M"
)

pars <- map_dbl(
  needed_params,
  ~ extract_mplus_param(param_df, .x)
)

names(pars) <- needed_params
pars <- as.list(pars)

print(pars)

################################################################################
# 4. Build adversity and moderator grids
################################################################################

x_grid <- seq(
  from = pars$XLOW,
  to   = pars$XHIGH,
  by   = x_step
)

# Make sure exact reference, knot, and high-adversity values are included.
x_grid <- sort(unique(c(x_grid, pars$XLOW, 0, pars$XHIGH)))

m_grid <- tibble(
  m_label = c("-2 SD", "-1.5 SD", "-1 SD", "-0.5 SD", "Mean",
              "+0.5 SD", "+1 SD", "+1.5 SD", "+2 SD"),
  m_value = c(pars$M2L, pars$M15L, pars$M1L, pars$M05L, pars$M0,
              pars$M05H, pars$M1H, pars$M15H, pars$M2H),
  m_order = 1:9
) %>%
  mutate(
    m_label = factor(m_label, levels = m_label)
  )

################################################################################
# 5. Compute expected trajectories across adversity x moderator grid
################################################################################

grid_traj <- compute_moderated_trajectory(
  x = x_grid,
  m = m_grid$m_value,
  time_values = time_values,
  pars = pars
) %>%
  left_join(m_grid, by = "m_value") %>%
  mutate(
    m_label = factor(m_label, levels = levels(m_grid$m_label))
  )

################################################################################
# 6. Attach fixed low-adversity reference trajectory
################################################################################

ref_traj <- compute_fixed_reference(
  time_values = time_values,
  pars = pars
)

grid_traj <- grid_traj %>%
  left_join(ref_traj, by = "time") %>%
  add_directional_difference(
    outcome_direction = outcome_direction,
    threshold_margin = threshold_margin
  )

################################################################################
# 7. Identify threshold within each moderator level
################################################################################

grid_summary <- summarise_grid_by_m_x(
  grid_df = grid_traj,
  threshold_margin = threshold_margin
)

threshold_result <- find_reference_return_average(
  grid_summary = grid_summary,
  threshold_margin = threshold_margin
) %>%
  mutate(
    outcome_direction = outcome_direction
  ) %>%
  arrange(m_value)

print(threshold_result)

################################################################################
# 8. Identify the earliest moderator level where crossover occurs
################################################################################

moderator_emergence <- threshold_result %>%
  filter(threshold_found) %>%
  arrange(m_value) %>%
  slice(1) %>%
  mutate(
    emergence_rule = "lowest_tested_moderator_level_with_average_reference_return"
  )

if (nrow(moderator_emergence) == 0) {
  moderator_emergence <- tibble(
    m_label = factor(NA_character_, levels = levels(m_grid$m_label)),
    m_value = NA_real_,
    threshold_found = FALSE,
    threshold_x = NA_real_,
    threshold_margin = threshold_margin,
    mean_maladaptive_diff_at_threshold = NA_real_,
    prop_time_crossed_at_threshold = NA_real_,
    n_time_crossed_at_threshold = NA_real_,
    n_time_total = NA_real_,
    threshold_rule = "average_maladaptive_difference",
    outcome_direction = outcome_direction,
    emergence_rule = "no_tested_moderator_level_with_average_reference_return"
  )
}

print(moderator_emergence)

################################################################################
# 9. Diagnostic summaries by time point and alternative rules
################################################################################

threshold_by_time <- find_threshold_by_time(
  grid_traj = grid_traj,
  threshold_margin = threshold_margin
) %>%
  mutate(outcome_direction = outcome_direction)

print(threshold_by_time)

alternative_thresholds <- find_alternative_thresholds(
  grid_summary = grid_summary,
  threshold_result = threshold_result
)

print(alternative_thresholds)

################################################################################
# 10. Create PieceOfTraj-ready output
################################################################################

pieceoftraj_df <- make_pieceoftraj_df(
  threshold_result = threshold_result,
  pars = pars,
  time_values = time_values
)

print(pieceoftraj_df)

################################################################################
# 11. Create plotting data for diagnostic trajectory plots
################################################################################

plot_x_values_by_m <- threshold_result %>%
  select(m_label, m_value, threshold_found, threshold_x) %>%
  group_by(m_label, m_value) %>%
  group_modify(~ {
    base_x <- seq(
      from = pars$XLOW,
      to   = pars$XHIGH,
      length.out = n_plot_trajectories
    )
    
    tibble(
      x = sort(unique(c(
        pars$XLOW,
        0,
        base_x,
        .x$threshold_x[.x$threshold_found],
        pars$XHIGH
      )))
    )
  }) %>%
  ungroup()

plot_traj <- plot_x_values_by_m %>%
  mutate(
    S1 = ifelse(x < 0, x, 0),
    S2 = ifelse(x > 0, x, 0),
    I_xm = pars$BMI +
      pars$BIX1*S1 + pars$BIX2*S2 +
      pars$BIM*m_value +
      pars$BIX1M*S1*m_value + pars$BIX2M*S2*m_value,
    S_xm = pars$BMS +
      pars$BSX1*S1 + pars$BSX2*S2 +
      pars$BSM*m_value +
      pars$BSX1M*S1*m_value + pars$BSX2M*S2*m_value
  ) %>%
  tidyr::crossing(time = time_values) %>%
  mutate(y_hat = I_xm + S_xm*time) %>%
  left_join(threshold_result %>% select(m_label, m_value, threshold_found, threshold_x),
            by = c("m_label", "m_value")) %>%
  mutate(
    x_label = case_when(
      abs(x - pars$XLOW) < 1e-8 ~ "Low adversity reference",
      abs(x - 0) < 1e-8 ~ "Knot / vertex",
      threshold_found & !is.na(threshold_x) & abs(x - threshold_x) < 1e-8 ~
        "Reference-return threshold",
      abs(x - pars$XHIGH) < 1e-8 ~ "Highest adversity",
      TRUE ~ paste0("x = ", round(x, 2))
    ),
    trajectory_type = case_when(
      x_label == "Low adversity reference" ~ "Reference",
      x_label == "Knot / vertex" ~ "Knot",
      x_label == "Reference-return threshold" ~ "Threshold",
      x_label == "Highest adversity" ~ "Highest adversity",
      TRUE ~ "Grid trajectory"
    ),
    is_key = trajectory_type %in% c(
      "Reference", "Knot", "Threshold", "Highest adversity"
    )
  )

label_df <- plot_traj %>%
  filter(is_key) %>%
  group_by(m_label, m_value, x_label) %>%
  filter(time == max(time)) %>%
  ungroup()

################################################################################
# 12. Plot expected trajectories across adversity, faceted by moderator
################################################################################

color_breaks <- c(
  pars$XLOW,
  pars$XLOW + 0.10 * (0 - pars$XLOW),
  pars$XLOW + 0.30 * (0 - pars$XLOW),
  pars$XLOW + 0.60 * (0 - pars$XLOW),
  0,
  pars$XHIGH * 0.35,
  pars$XHIGH * 0.70,
  pars$XHIGH
)

color_values <- rescale(color_breaks)

p_traj <- ggplot(
  plot_traj,
  aes(x = time, y = y_hat, group = interaction(m_label, x))
) +
  geom_line(
    data = subset(plot_traj, !is_key),
    aes(color = x),
    linewidth = 0.55,
    alpha = 0.40,
    lineend = "round"
  ) +
  geom_line(
    data = subset(plot_traj, is_key),
    aes(color = x),
    linewidth = 1.20,
    alpha = 1,
    lineend = "round"
  ) +
  geom_point(
    data = subset(plot_traj, is_key),
    aes(color = x),
    size = 2.1
  ) +
  geom_label_repel(
    data = label_df,
    aes(label = x_label, color = x),
    nudge_x = 0.30,
    direction = "y",
    hjust = 0,
    box.padding = 0.20,
    point.padding = 0.15,
    segment.alpha = 0.45,
    segment.size = 0.35,
    fill = "white",
    size = 2.8,
    family = "Arial",
    show.legend = FALSE
  ) +
  scale_color_gradientn(
    colors = c(
      "#2B3A8C", "#3B6FB6", "#4FA3C7", "#86CFA5",
      "#C7E77A", "#F2C14E", "#F28E2B", "#D55E00"
    ),
    values = color_values,
    breaks = pretty(c(pars$XLOW, pars$XHIGH), n = 6),
    name = "Family conflict"
  ) +
  scale_x_continuous(
    breaks = time_values,
    expand = expansion(mult = c(0.03, 0.32))
  ) +
  facet_wrap(~ m_label, nrow = 3) +
  theme_classic(base_size = 12) +
  theme(
    text = element_text(family = "Arial"),
    plot.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold"),
    legend.position = "right",
    plot.margin = margin(5.5, 100, 5.5, 5.5)
  ) +
  labs(
    title = "Expected Growth Trajectories Across school adveristy and Parental Monitoring",
    subtitle = paste0(
      "Reference is fixed at low school adveristy with moderators/covariates held at mean/zero. ",
      "Threshold rule: average maladaptive difference."
    ),
    x = "Time",
    y = "DERS Goals"
  ) +
  coord_cartesian(clip = "off")

print(p_traj)

ggsave("GOALS_SCHADVxPMON.png",
       p_traj, width = 11, height = 8.5, dpi = 300, bg = "white")

################################################################################
# 13. Plot average maladaptive difference by moderator
################################################################################

p_mean_diff <- grid_summary %>%
  ggplot(aes(x = x, y = mean_maladaptive_diff)) +
  geom_hline(
    yintercept = threshold_margin,
    linetype = "dashed",
    linewidth = 0.70
  ) +
  geom_line(
    linewidth = 0.95,
    alpha = 0.90
  ) +
  geom_vline(
    data = threshold_result %>% filter(threshold_found),
    aes(xintercept = threshold_x),
    linetype = "dotted",
    linewidth = 0.80
  ) +
  facet_wrap(~ m_label, nrow = 3) +
  theme_classic(base_size = 12) +
  theme(
    text = element_text(family = "Arial"),
    plot.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold")
  ) +
  labs(
    title = "Average Maladaptive Difference from Fixed Low-Adversity Reference",
    subtitle = paste0(
      "Vertical dotted line is shown only when a threshold is found. ",
      "Outcome direction: ", outcome_direction
    ),
    x = "Post-knot family-conflict value",
    y = "Mean maladaptive difference from reference"
  )

print(p_mean_diff)

################################################################################
# 14. Plot time-specific differences by moderator
################################################################################

p_diff <- grid_traj %>%
  filter(x >= 0) %>%
  ggplot(aes(x = x, y = maladaptive_diff, group = factor(time))) +
  geom_hline(
    yintercept = threshold_margin,
    linetype = "dashed",
    linewidth = 0.65
  ) +
  geom_line(
    linewidth = 0.80,
    alpha = 0.88
  ) +
  geom_vline(
    data = threshold_result %>% filter(threshold_found),
    aes(xintercept = threshold_x),
    linetype = "dotted",
    linewidth = 0.75
  ) +
  facet_wrap(~ m_label, nrow = 3) +
  theme_classic(base_size = 12) +
  theme(
    text = element_text(family = "Arial"),
    plot.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold")
  ) +
  labs(
    title = "Time-Specific Maladaptive Difference from Fixed Reference",
    subtitle = "Positive values indicate movement beyond the fixed reference in the maladaptive direction.",
    x = "Post-knot family-conflict value",
    y = "Maladaptive difference from reference",
    group = "Time"
  )

print(p_diff)

################################################################################
# 15. Inspect exact grid rows at identified thresholds
################################################################################

threshold_rows <- grid_traj %>%
  semi_join(
    threshold_result %>% filter(threshold_found) %>% select(m_label, m_value, threshold_x),
    by = c("m_label", "m_value")
  ) %>%
  inner_join(
    threshold_result %>% filter(threshold_found) %>% select(m_label, m_value, threshold_x),
    by = c("m_label", "m_value")
  ) %>%
  filter(abs(x - threshold_x) < 1e-10)

print(threshold_rows)

threshold_summary_rows <- grid_summary %>%
  inner_join(
    threshold_result %>% filter(threshold_found) %>% select(m_label, m_value, threshold_x),
    by = c("m_label", "m_value")
  ) %>%
  filter(abs(x - threshold_x) < 1e-10)

print(threshold_summary_rows)

################################################################################
# 16. Save outputs
################################################################################

if (save_outputs) {
  
  write.csv(grid_traj,
            paste0(out_prefix, "_grid_trajectory_values.csv"),
            row.names = FALSE)
  
  write.csv(grid_summary,
            paste0(out_prefix, "_grid_summary_by_moderator_adversity.csv"),
            row.names = FALSE)
  
  write.csv(threshold_result,
            paste0(out_prefix, "_threshold_by_moderator.csv"),
            row.names = FALSE)
  
  write.csv(moderator_emergence,
            paste0(out_prefix, "_moderator_emergence_summary.csv"),
            row.names = FALSE)
  
  write.csv(threshold_by_time,
            paste0(out_prefix, "_threshold_by_time.csv"),
            row.names = FALSE)
  
  write.csv(alternative_thresholds,
            paste0(out_prefix, "_alternative_thresholds.csv"),
            row.names = FALSE)
  
  write.csv(pieceoftraj_df,
            paste0(out_prefix, "_PieceOfTraj_input.csv"),
            row.names = FALSE)
  
  ggsave(paste0(out_prefix, "_expected_trajectories.png"),
         p_traj, width = 11, height = 8.5, dpi = 300, bg = "white")
  
  ggsave(paste0(out_prefix, "_average_maladaptive_difference.png"),
         p_mean_diff, width = 11, height = 8.5, dpi = 300, bg = "white")
  
  ggsave(paste0(out_prefix, "_time_specific_differences.png"),
         p_diff, width = 11, height = 8.5, dpi = 300, bg = "white")
}
