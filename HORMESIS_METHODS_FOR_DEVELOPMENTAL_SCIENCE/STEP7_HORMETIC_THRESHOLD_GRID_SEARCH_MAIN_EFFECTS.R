################################################################################
# Hormetic Reference-Return Threshold Grid Search
# Mplus LGC spline output -> R grid search
#
# Purpose:
#   Identify the first post-knot adversity value where the model-implied
#   trajectory returns to, or crosses beyond, the low-adversity reference
#   trajectory in the direction of worse functioning.
#
# Primary threshold rule:
#   Area-under-difference / average maladaptive difference.
#
#   The threshold is the first post-knot adversity value where the expected
#   trajectory is, on average across modeled time, at or beyond the
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

# Path to your Mplus output file
out_file <- "C:\\Users\\cjh37695\\Dropbox\\HORMESIS_METHODS\\ANALYSIS\\STRUCTURAL\\reapp_lineargrowth_hses_threshold_search.out"

# Growth-model time scores.
# Update this if your LGC loadings are not 0, 1, 2.
time_values <- c(0, 1, 2)

# Adversity-grid resolution.
# Smaller = more precise but larger data frame.
x_step <- 0.001

# Direction of worse functioning for the outcome.
#
# Use:
#   "higher_is_worse" for symptoms, dysregulation, impairment, risk, etc.
#   "lower_is_worse"  for competence, adaptive functioning, wellbeing,
#                     neural integrity, etc.
#
outcome_direction <- "higher_is_worse" ## flipping this around because it is post-traumatic growth

# Optional tolerance around exact equality.
#
# Keep at 0 for the pure mathematical reference-return threshold.
# Use a positive value only if you want a more conservative threshold,
# e.g., requiring the average maladaptive difference to exceed the reference
# by a substantively meaningful amount.
#
# Example:
#   threshold_margin <- 0.05
#
threshold_margin <- 0

# Number of x-values to show in the trajectory plot.
# This affects plotting only, not the threshold search.
n_plot_trajectories <- 45

################################################################################
# 2. Helper functions
################################################################################

extract_mplus_param <- function(param_df, param_name) {
  
  param_df_clean <- param_df %>%
    mutate(
      param_upper = toupper(.data$param)
    )
  
  out <- param_df_clean %>%
    filter(.data$param_upper == toupper(param_name))
  
  if (nrow(out) == 0) {
    stop(paste0("Could not find parameter: ", param_name))
  }
  
  if (nrow(out) > 1) {
    warning(paste0(
      "More than one row found for parameter: ", param_name,
      ". Using the first row."
    ))
  }
  
  out$est[1]
}


make_spline_terms <- function(x) {
  
  tibble(
    x = x,
    S1 = ifelse(x < 0, x, 0),
    S2 = ifelse(x > 0, x, 0)
  )
}


compute_trajectory <- function(x, time_values, pars) {
  
  spline_df <- make_spline_terms(x)
  
  traj_df <- spline_df %>%
    mutate(
      I_x = pars$B_MUI + pars$B_IX1*S1 + pars$B_IX2*S2,
      S_x = pars$B_MUS + pars$B_SX1*S1 + pars$B_SX2*S2
    ) %>%
    tidyr::crossing(time = time_values) %>%
    mutate(
      y_hat = I_x + S_x*time
    )
  
  traj_df
}


add_directional_difference <- function(grid_df,
                                       outcome_direction = c("higher_is_worse",
                                                             "lower_is_worse")) {
  
  outcome_direction <- match.arg(outcome_direction)
  
  grid_df %>%
    mutate(
      # Raw difference from the low-adversity reference trajectory.
      # Positive means y_hat is above reference.
      # Negative means y_hat is below reference.
      diff_from_ref = y_hat - y_ref,
      
      # Directional difference where positive values always mean
      # worse-than-reference movement.
      maladaptive_diff = case_when(
        outcome_direction == "higher_is_worse" ~ diff_from_ref,
        outcome_direction == "lower_is_worse"  ~ -diff_from_ref,
        TRUE ~ NA_real_
      ),
      
      # Time-specific crossing, retained as diagnostic information.
      crossed_reference = maladaptive_diff >= threshold_margin,
      
      direction_label = case_when(
        outcome_direction == "higher_is_worse" ~ "Higher values indicate worse functioning",
        outcome_direction == "lower_is_worse"  ~ "Lower values indicate worse functioning",
        TRUE ~ NA_character_
      )
    )
}


summarise_grid_by_x <- function(grid_df, threshold_margin = 0) {
  
  grid_df %>%
    filter(x >= 0) %>%
    group_by(x) %>%
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


find_reference_return_average <- function(grid_summary,
                                          threshold_margin = 0) {
  
  threshold_df <- grid_summary %>%
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
}

################################################################################
# 3. Read Mplus output and extract needed parameters
################################################################################

mplus_out <- readModels(out_file, what = "parameters")

param_df <- mplus_out$parameters$unstandardized

needed_params <- c(
  "XLOW", "XHIGH",
  "B_MUI", "B_MUS",
  "B_IX1", "B_IX2",
  "B_SX1", "B_SX2"
)

pars <- map_dbl(
  needed_params,
  ~ extract_mplus_param(param_df, .x)
)

names(pars) <- needed_params

pars <- as.list(pars)

print(pars)

################################################################################
# 4. Build adversity grid
################################################################################

x_grid <- seq(
  from = pars$XLOW,
  to   = pars$XHIGH,
  by   = x_step
)

# Make sure exact reference, knot, and high-adversity values are included.
x_grid <- sort(unique(c(x_grid, pars$XLOW, 0, pars$XHIGH)))

################################################################################
# 5. Compute expected trajectories across adversity grid
################################################################################

grid_traj <- compute_trajectory(
  x = x_grid,
  time_values = time_values,
  pars = pars
)

################################################################################
# 6. Compute low-adversity reference trajectory
################################################################################

ref_traj <- compute_trajectory(
  x = pars$XLOW,
  time_values = time_values,
  pars = pars
) %>%
  select(time, y_ref = y_hat)

grid_traj <- grid_traj %>%
  left_join(ref_traj, by = "time") %>%
  add_directional_difference(
    outcome_direction = outcome_direction
  )

################################################################################
# 7. Identify threshold using average maladaptive difference
################################################################################

grid_summary <- summarise_grid_by_x(
  grid_df = grid_traj,
  threshold_margin = threshold_margin
)

threshold_result <- find_reference_return_average(
  grid_summary = grid_summary,
  threshold_margin = threshold_margin
) %>%
  mutate(
    outcome_direction = outcome_direction
  )

print(threshold_result)

################################################################################
# 8. Diagnostic summaries by time point
################################################################################

# This is no longer the primary threshold rule.
# It is retained to show when each modeled occasion crosses the reference.
threshold_by_time <- grid_traj %>%
  filter(x >= 0) %>%
  group_by(time) %>%
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
    outcome_direction = outcome_direction,
    threshold_margin = threshold_margin
  )

print(threshold_by_time)

################################################################################
# 9. Optional diagnostic: compare alternative rules
################################################################################

alternative_thresholds <- tibble(
  rule = c("any_time", "majority_time", "all_time", "average_maladaptive_difference"),
  threshold_x = c(
    grid_summary %>%
      filter(crossed_any_time) %>%
      summarise(x = ifelse(n() > 0, min(x), NA_real_)) %>%
      pull(x),
    
    grid_summary %>%
      filter(prop_time_crossed >= 0.50) %>%
      summarise(x = ifelse(n() > 0, min(x), NA_real_)) %>%
      pull(x),
    
    grid_summary %>%
      filter(crossed_all_time) %>%
      summarise(x = ifelse(n() > 0, min(x), NA_real_)) %>%
      pull(x),
    
    threshold_result$threshold_x
  )
)

print(alternative_thresholds)

################################################################################
# 10. Create plotting data
################################################################################

plot_x_values <- seq(
  from = pars$XLOW,
  to   = pars$XHIGH,
  length.out = n_plot_trajectories
)

plot_x_values <- sort(unique(c(
  pars$XLOW,
  0,
  plot_x_values,
  threshold_result$threshold_x,
  pars$XHIGH
)))

plot_traj <- compute_trajectory(
  x = plot_x_values,
  time_values = time_values,
  pars = pars
) %>%
  mutate(
    x_label = case_when(
      abs(x - pars$XLOW) < 1e-8 ~ "Low adversity reference",
      abs(x - 0) < 1e-8 ~ "Knot / vertex",
      !is.na(threshold_result$threshold_x) &
        abs(x - threshold_result$threshold_x) < 1e-8 ~ "Reference-return threshold",
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
      "Reference",
      "Knot",
      "Threshold",
      "Highest adversity"
    )
  )

# Label only key trajectories at the final time point.
label_df <- plot_traj %>%
  group_by(x_label) %>%
  filter(time == max(time)) %>%
  ungroup() %>%
  filter(is_key)

################################################################################
# 11. Plot expected trajectories with distinct adversity gradient
################################################################################

# Custom gradient anchors.
# These intentionally allocate more color space around the low-to-knot region
# than a default continuous palette usually does.
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
  aes(x = time, y = y_hat, group = x)
) +
  # background trajectories
  geom_line(
    data = subset(plot_traj, !is_key),
    aes(color = x),
    linewidth = 0.65,
    alpha = 0.45,
    lineend = "round"
  ) +
  
  # key trajectories
  geom_line(
    data = subset(plot_traj, is_key),
    aes(color = x),
    linewidth = 1.35,
    alpha = 1,
    lineend = "round"
  ) +
  
  geom_point(
    data = subset(plot_traj, is_key),
    aes(color = x),
    size = 2.6
  ) +
  
  # repelled labels off to the right
  geom_label_repel(
    data = label_df,
    aes(
      label = x_label,
      color = x
    ),
    nudge_x = 0.35,
    direction = "y",
    hjust = 0,
    box.padding = 0.25,
    point.padding = 0.20,
    segment.alpha = 0.50,
    segment.size = 0.40,
    fill = "white",
    size = 3.4,
    family = "Arial",
    show.legend = FALSE
  ) +
  
  scale_color_gradientn(
    colors = c(
      "#2B3A8C",  # deep indigo-blue
      "#3B6FB6",  # blue
      "#4FA3C7",  # blue-teal
      "#86CFA5",  # green
      "#C7E77A",  # yellow-green around knot
      "#F2C14E",  # gold
      "#F28E2B",  # orange
      "#D55E00"   # red-orange
    ),
    values = color_values,
    breaks = pretty(c(pars$XLOW, pars$XHIGH), n = 6),
    name = "Adversity"
  ) +
  
  scale_x_continuous(
    breaks = time_values,
    expand = expansion(mult = c(0.03, 0.30))
  ) +
  
  theme_classic(base_size = 13) +
  theme(
    text = element_text(family = "Arial"),
    plot.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    legend.position = "right",
    plot.margin = margin(5.5, 100, 5.5, 5.5)
  ) +
  labs(
    title = "Expected Growth Trajectories Across Adversity Levels",
    subtitle = paste0(
      "Highlighted lines indicate the low-adversity reference, knot, reference-return threshold, and highest adversity.\n",
      "Threshold rule: average maladaptive difference; threshold x = ",
      ifelse(
        threshold_result$threshold_found,
        round(threshold_result$threshold_x, 3),
        "not found"
      )
    ),
    x = "Time",
    y = "Expected outcome trajectory"
  ) +
  coord_cartesian(clip = "off")

print(p_traj)

################################################################################
# 12. Plot average maladaptive difference across adversity
################################################################################

p_mean_diff <- grid_summary %>%
  ggplot(aes(x = x, y = mean_maladaptive_diff)) +
  geom_hline(
    yintercept = threshold_margin,
    linetype = "dashed",
    linewidth = 0.75
  ) +
  geom_line(
    linewidth = 1.05,
    alpha = 0.90
  ) +
  geom_vline(
    data = threshold_result %>% filter(threshold_found),
    aes(xintercept = threshold_x),
    linetype = "dotted",
    linewidth = 0.85
  ) +
  theme_classic(base_size = 13) +
  theme(
    text = element_text(family = "Arial"),
    plot.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold")
  ) +
  labs(
    title = "Average Maladaptive Difference from Low-Adversity Reference",
    subtitle = paste0(
      "Threshold is the first post-knot adversity value where the time-averaged difference reaches the reference.\n",
      "Outcome direction: ", outcome_direction,
      "; threshold x = ",
      ifelse(
        threshold_result$threshold_found,
        round(threshold_result$threshold_x, 3),
        "not found"
      )
    ),
    x = "Post-knot adversity value",
    y = "Mean maladaptive difference from reference"
  )

print(p_mean_diff)

################################################################################
# 13. Plot time-specific differences from reference
################################################################################

p_diff <- grid_traj %>%
  filter(x >= 0) %>%
  ggplot(aes(x = x, y = maladaptive_diff, group = factor(time))) +
  geom_hline(
    yintercept = threshold_margin,
    linetype = "dashed",
    linewidth = 0.70
  ) +
  geom_line(
    linewidth = 0.90,
    alpha = 0.90
  ) +
  geom_vline(
    data = threshold_result %>% filter(threshold_found),
    aes(xintercept = threshold_x),
    linetype = "dotted",
    linewidth = 0.80
  ) +
  theme_classic(base_size = 13) +
  theme(
    text = element_text(family = "Arial"),
    plot.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold")
  ) +
  labs(
    title = "Time-Specific Maladaptive Difference from Low-Adversity Reference",
    subtitle = paste0(
      "Positive values indicate movement beyond the reference in the maladaptive direction. ",
      "Outcome direction: ", outcome_direction
    ),
    x = "Post-knot adversity value",
    y = "Maladaptive difference from reference",
    group = "Time"
  )

print(p_diff)

################################################################################
# 14. Inspect exact grid rows at the threshold
################################################################################

threshold_rows <- grid_traj %>%
  filter(
    threshold_result$threshold_found,
    abs(x - threshold_result$threshold_x) < 1e-10
  )

print(threshold_rows)

threshold_summary_row <- grid_summary %>%
  filter(
    threshold_result$threshold_found,
    abs(x - threshold_result$threshold_x) < 1e-10
  )

print(threshold_summary_row)

################################################################################
# 15. Optional: save outputs
################################################################################

# write.csv(grid_traj, "hormetic_grid_trajectory_values.csv", row.names = FALSE)
# write.csv(grid_summary, "hormetic_grid_summary_by_adversity.csv", row.names = FALSE)
# write.csv(threshold_result, "hormetic_reference_return_threshold.csv", row.names = FALSE)
# write.csv(threshold_by_time, "hormetic_threshold_by_time.csv", row.names = FALSE)
# write.csv(alternative_thresholds, "hormetic_alternative_thresholds.csv", row.names = FALSE)

# ggsave("hormetic_expected_trajectories.png", p_traj, width = 7.5, height = 5, dpi = 300)
# ggsave("hormetic_average_maladaptive_difference.png", p_mean_diff, width = 7, height = 5, dpi = 300)
# ggsave("hormetic_time_specific_differences.png", p_diff, width = 7, height = 5, dpi = 300)