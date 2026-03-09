#Library

library(dplyr)
library(MplusAutomation)

# Set working directory 

setwd("/home/cjh37695/ABCD_PROJECTS/HORMESIS_METHODS/")

#Library

df <- read.csv("ABCD_HORM_METH_12.29.25.csv")


################################################################################
##################### PREP ADVERSITY VARIABLES  
names(df)

### EARLY LIFE COGNITIVE RESOURCES (HOUSEHOLD SHORT FORM)
# THESE ARE SCALED WEIRD, SOME ARE 0 - 3 AND SOME 0 - 5. 
# Coded so that bad are high

df$ErlHs1PR <- 
  dplyr::recode(df$ErlHs1P_11,
                `0` = 4,
                `1` = 3,
                `2` = 2,
                `3` = 1)
df$ErlHs2PR <- 
  dplyr::recode(df$ErlHs2P_11,
                `0` = 4,
                `1` = 3,
                `2` = 2,
                `3` = 1)
df$ErlHs3PR <- 
  dplyr::recode(df$ErlHs3P_11,
                `0` = 4,
                `1` = 3,
                `2` = 2,
                `3` = 1)
df$ErlHs4PR <- 
  dplyr::recode(df$ErlHs4P_11,
                `0` = 6,
                `1` = 5,
                `2` = 4,
                `3` = 3,
                `4` = 2,
                `5` = 1)
df$ErlHs5PR <- 
  dplyr::recode(df$ErlHs5P_11,
                `0` = 4,
                `1` = 3,
                `2` = 2,
                `3` = 1)
# ITEM 6
df$ErlHs6PR <- 
  dplyr::recode(df$ErlHs6P_11,
                `0` = 4,
                `1` = 3,
                `2` = 2,
                `3` = 1)

### NEIGHBORHOOD SAFETY
#W1
df$NBHsaf1P_1R <- 
  max(df$NBHsaf1P_1, na.rm = TRUE) +
  min(df$NBHsaf1P_1, na.rm = TRUE) -
  df$NBHsaf1P_1
df$NBHsaf2P_1R <- 
  max(df$NBHsaf2P_1, na.rm = TRUE) +
  min(df$NBHsaf2P_1, na.rm = TRUE) -
  df$NBHsaf2P_1
df$NBHsaf3P_1R <- 
  max(df$NBHsaf3P_1, na.rm = TRUE) +
  min(df$NBHsaf3P_1, na.rm = TRUE) -
  df$NBHsaf3P_1
#W3
df$NBHsaf1P_3R <- 
  max(df$NBHsaf1P_3, na.rm = TRUE) +
  min(df$NBHsaf1P_3, na.rm = TRUE) -
  df$NBHsaf1P_3
df$NBHsaf2P_3R <- 
  max(df$NBHsaf2P_3, na.rm = TRUE) +
  min(df$NBHsaf2P_3, na.rm = TRUE) -
  df$NBHsaf2P_3
df$NBHsaf3P_3R <- 
  max(df$NBHsaf3P_3, na.rm = TRUE) +
  min(df$NBHsaf3P_3, na.rm = TRUE) -
  df$NBHsaf3P_3
#W5
df$NBHsaf1P_5R <- 
  max(df$NBHsaf1P_5, na.rm = TRUE) +
  min(df$NBHsaf1P_5, na.rm = TRUE) -
  df$NBHsaf1P_5
df$NBHsaf2P_5R <- 
  max(df$NBHsaf2P_5, na.rm = TRUE) +
  min(df$NBHsaf2P_5, na.rm = TRUE) -
  df$NBHsaf2P_5
df$NBHsaf3P_5R <- 
  max(df$NBHsaf3P_5, na.rm = TRUE) +
  min(df$NBHsaf3P_5, na.rm = TRUE) -
  df$NBHsaf3P_5
#W7
df$NBHsaf1P_7R <- 
  max(df$NBHsaf1P_7, na.rm = TRUE) +
  min(df$NBHsaf1P_7, na.rm = TRUE) -
  df$NBHsaf1P_7
df$NBHsaf2P_7R <- 
  max(df$NBHsaf2P_7, na.rm = TRUE) +
  min(df$NBHsaf2P_7, na.rm = TRUE) -
  df$NBHsaf2P_7
df$NBHsaf3P_7R <- 
  max(df$NBHsaf3P_7, na.rm = TRUE) +
  min(df$NBHsaf3P_7, na.rm = TRUE) -
  df$NBHsaf3P_7
#W9
df$NBHsaf1P_9R <- 
  max(df$NBHsaf1P_9, na.rm = TRUE) +
  min(df$NBHsaf1P_9, na.rm = TRUE) -
  df$NBHsaf1P_9
df$NBHsaf2P_9R <- 
  max(df$NBHsaf2P_9, na.rm = TRUE) +
  min(df$NBHsaf2P_9, na.rm = TRUE) -
  df$NBHsaf2P_9
df$NBHsaf3P_9R <- 
  max(df$NBHsaf3P_9, na.rm = TRUE) +
  min(df$NBHsaf3P_9, na.rm = TRUE) -
  df$NBHsaf3P_9
#W11
df$NBHsaf1P_11R <- 
  max(df$NBHsaf1P_11, na.rm = TRUE) +
  min(df$NBHsaf1P_11, na.rm = TRUE) -
  df$NBHsaf1P_11
df$NBHsaf2P_11R <- 
  max(df$NBHsaf2P_11, na.rm = TRUE) +
  min(df$NBHsaf2P_11, na.rm = TRUE) -
  df$NBHsaf2P_11
df$NBHsaf3P_11R <- 
  max(df$NBHsaf3P_11, na.rm = TRUE) +
  min(df$NBHsaf3P_11, na.rm = TRUE) -
  df$NBHsaf3P_11

### SCHOOL SAFETY
#W1
df$SCHenv1Y_1R <- 
  max(df$SCHenv1Y_1, na.rm = TRUE) +
  min(df$SCHenv1Y_1, na.rm = TRUE) -
  df$SCHenv1Y_1
df$SCHenv2Y_1R <- 
  max(df$SCHenv2Y_1, na.rm = TRUE) +
  min(df$SCHenv2Y_1, na.rm = TRUE) -
  df$SCHenv2Y_1
df$SCHenv3Y_1R <- 
  max(df$SCHenv3Y_1, na.rm = TRUE) +
  min(df$SCHenv3Y_1, na.rm = TRUE) -
  df$SCHenv3Y_1
df$SCHenv4Y_1R <- 
  max(df$SCHenv4Y_1, na.rm = TRUE) +
  min(df$SCHenv4Y_1, na.rm = TRUE) -
  df$SCHenv4Y_1
df$SCHenv5Y_1R <- 
  max(df$SCHenv5Y_1, na.rm = TRUE) +
  min(df$SCHenv5Y_1, na.rm = TRUE) -
  df$SCHenv5Y_1
df$SCHenv6Y_1R <- 
  max(df$SCHenv6Y_1, na.rm = TRUE) +
  min(df$SCHenv6Y_1, na.rm = TRUE) -
  df$SCHenv6Y_1
#W3
df$SCHenv1Y_3R <- 
  max(df$SCHenv1Y_3, na.rm = TRUE) +
  min(df$SCHenv1Y_3, na.rm = TRUE) -
  df$SCHenv1Y_3
df$SCHenv2Y_3R <- 
  max(df$SCHenv2Y_3, na.rm = TRUE) +
  min(df$SCHenv2Y_3, na.rm = TRUE) -
  df$SCHenv2Y_3
df$SCHenv3Y_3R <- 
  max(df$SCHenv3Y_3, na.rm = TRUE) +
  min(df$SCHenv3Y_3, na.rm = TRUE) -
  df$SCHenv3Y_3
df$SCHenv4Y_3R <- 
  max(df$SCHenv4Y_3, na.rm = TRUE) +
  min(df$SCHenv4Y_3, na.rm = TRUE) -
  df$SCHenv4Y_3
df$SCHenv5Y_3R <- 
  max(df$SCHenv5Y_3, na.rm = TRUE) +
  min(df$SCHenv5Y_3, na.rm = TRUE) -
  df$SCHenv5Y_3
df$SCHenv6Y_3R <- 
  max(df$SCHenv6Y_3, na.rm = TRUE) +
  min(df$SCHenv6Y_3, na.rm = TRUE) -
  df$SCHenv6Y_3
#W5
df$SCHenv1Y_5R <- 
  max(df$SCHenv1Y_5, na.rm = TRUE) +
  min(df$SCHenv1Y_5, na.rm = TRUE) -
  df$SCHenv1Y_5
df$SCHenv2Y_5R <- 
  max(df$SCHenv2Y_5, na.rm = TRUE) +
  min(df$SCHenv2Y_5, na.rm = TRUE) -
  df$SCHenv2Y_5
df$SCHenv3Y_5R <- 
  max(df$SCHenv3Y_5, na.rm = TRUE) +
  min(df$SCHenv3Y_5, na.rm = TRUE) -
  df$SCHenv3Y_5
df$SCHenv4Y_5R <- 
  max(df$SCHenv4Y_5, na.rm = TRUE) +
  min(df$SCHenv4Y_5, na.rm = TRUE) -
  df$SCHenv4Y_5
df$SCHenv5Y_5R <- 
  max(df$SCHenv5Y_5, na.rm = TRUE) +
  min(df$SCHenv5Y_5, na.rm = TRUE) -
  df$SCHenv5Y_5
df$SCHenv6Y_5R <- 
  max(df$SCHenv6Y_5, na.rm = TRUE) +
  min(df$SCHenv6Y_5, na.rm = TRUE) -
  df$SCHenv6Y_5
#W7
df$SCHenv1Y_7R <- 
  max(df$SCHenv1Y_7, na.rm = TRUE) +
  min(df$SCHenv1Y_7, na.rm = TRUE) -
  df$SCHenv1Y_7
df$SCHenv2Y_7R <- 
  max(df$SCHenv2Y_7, na.rm = TRUE) +
  min(df$SCHenv2Y_7, na.rm = TRUE) -
  df$SCHenv2Y_7
df$SCHenv3Y_7R <- 
  max(df$SCHenv3Y_7, na.rm = TRUE) +
  min(df$SCHenv3Y_7, na.rm = TRUE) -
  df$SCHenv3Y_7
df$SCHenv4Y_7R <- 
  max(df$SCHenv4Y_7, na.rm = TRUE) +
  min(df$SCHenv4Y_7, na.rm = TRUE) -
  df$SCHenv4Y_7
df$SCHenv5Y_7R <- 
  max(df$SCHenv5Y_7, na.rm = TRUE) +
  min(df$SCHenv5Y_7, na.rm = TRUE) -
  df$SCHenv5Y_7
df$SCHenv6Y_7R <- 
  max(df$SCHenv6Y_7, na.rm = TRUE) +
  min(df$SCHenv6Y_7, na.rm = TRUE) -
  df$SCHenv6Y_7
#W9
df$SCHenv1Y_9R <- 
  max(df$SCHenv1Y_9, na.rm = TRUE) +
  min(df$SCHenv1Y_9, na.rm = TRUE) -
  df$SCHenv1Y_9
df$SCHenv2Y_9R <- 
  max(df$SCHenv2Y_9, na.rm = TRUE) +
  min(df$SCHenv2Y_9, na.rm = TRUE) -
  df$SCHenv2Y_9
df$SCHenv3Y_9R <- 
  max(df$SCHenv3Y_9, na.rm = TRUE) +
  min(df$SCHenv3Y_9, na.rm = TRUE) -
  df$SCHenv3Y_9
df$SCHenv4Y_9R <- 
  max(df$SCHenv4Y_9, na.rm = TRUE) +
  min(df$SCHenv4Y_9, na.rm = TRUE) -
  df$SCHenv4Y_9
df$SCHenv5Y_9R <- 
  max(df$SCHenv5Y_9, na.rm = TRUE) +
  min(df$SCHenv5Y_9, na.rm = TRUE) -
  df$SCHenv5Y_9
df$SCHenv6Y_9R <- 
  max(df$SCHenv6Y_9, na.rm = TRUE) +
  min(df$SCHenv6Y_9, na.rm = TRUE) -
  df$SCHenv6Y_9

##### HARSH SES 

#Reverse Score SOCIAL MOBILITY 
df$SocMob_1R <- 
  max(df$SocMob_1, na.rm = TRUE) +
  min(df$SocMob_1, na.rm = TRUE) -
  df$SocMob_1



# Reverse score parents highest education
df$HiParEdu_1R <- 
  dplyr::recode(df$HiParEdu_1,
                `1` = 5,
                `2` = 4,
                `3` = 3,
                `4` = 2,
                `5` = 1)



##### COMPUTE FAMILY CONFLICT 
# HELPER FUNCTION FOR CONFLICT SUM 

compute_scale <- function(df, items, max_na = 1) {
  na_count <- rowSums(is.na(df[items]))
  score    <- rowSums(df[items], na.rm = TRUE)
  ifelse(na_count <= max_na, score, NA_real_)
}

# DEFINE WAVES WE WANT 
waves <- c(1, 3, 5, 7, 9, 11)

# PARENT REPORT 
## REMOVING ITEMS 7 & 9 BECAUSE IT IS KNOWN THEY DO NOT FIT WELL
for (w in waves) {
  items <- c(
    paste0("Fcon1P_", w),
    paste0("Fcon2rP_", w),
    paste0("Fcon3P_", w),
    paste0("Fcon4rP_", w),
    paste0("Fcon5P_", w),
    paste0("Fcon6P_", w),
#    paste0("Fcon7rP_", w),
    paste0("Fcon8P_", w)#,
#    paste0("Fcon8P.1_", w)
  )
  
  df[[paste0("FAMCONP", w)]] <- compute_scale(df, items)
}


# YOUTH REPORT 
## REMOVING ITEMS 7 & 9 BECAUSE IT IS KNOWN THEY DO NOT FIT WELL
for (w in waves) {
  items <- c(
    paste0("Fcon1Y_", w),
    paste0("Fcon2rY_", w),
    paste0("Fcon3Y_", w),
    paste0("Fcon4rY_", w),
    paste0("Fcon5Y_", w),
    paste0("Fcon6Y_", w),
#    paste0("Fcon7rY_", w),
    paste0("Fcon8Y_", w)#,
#    paste0("Fcon9Y_", w)
  )
  
  df[[paste0("FAMCONY", w)]] <- compute_scale(df, items)
}

## TODO: NEED TO COMPUTE INTERNAL RELIABILITY


################################################################################
##################### DROP NEIGHBORHOOD VALUES NOT AT EARLY TIMEPOINTS 

df_RED <- df %>%
  select(-c("NBHcoh1P_5",    "NBHcoh1P_9",   
            "NBHcoh2P_5",    "NBHcoh2P_9",    "NBHcoh3P_5",    "NBHcoh3P_9",    "NBHcoh4P_5",    "NBHcoh4P_9",   
            "NBHcoh5P_5",    "NBHcoh5P_9",    "NBHisc1P_5",    "NBHisc1P_9",    "NBHisc2P_5",    "NBHisc2P_9",   
            "NBHisc3P_5",    "NBHisc3P_9",    "NBHisc4P_5",    "NBHisc4P_9",    "NBHisc5P_5",    "NBHisc5P_9"))


################################################################################
#####################  Visualize data distributions 

## ASSESS LONGITUDINAL DISTRIBUTIONS 

plot_longitudinal_change <- function(df, columns, ncol = NULL, max_per_page = 9, 
                                     max_participants = NULL) {
  
  # Load required libraries
  require(ggplot2)
  require(tidyr)
  require(dplyr)
  require(patchwork)
  
  # Validate inputs
  if (!is.numeric(columns)) {
    stop("'columns' must be a numeric vector of column positions")
  }
  
  if (any(columns < 1 | columns > ncol(df))) {
    stop(paste("Column positions must be between 1 and", ncol(df)))
  }
  
  # Get column names
  column_names <- names(df)[columns]
  
  # Check if columns are numeric
  numeric_check <- sapply(df[, columns, drop = FALSE], is.numeric)
  if (!all(numeric_check)) {
    non_numeric <- column_names[!numeric_check]
    warning(sprintf("Excluding %d non-numeric column(s): %s", 
                    length(non_numeric), 
                    paste(head(non_numeric, 3), collapse = ", ")))
    column_names <- column_names[numeric_check]
    if (length(column_names) == 0) {
      stop("No numeric columns found in the specified columns.")
    }
  }
  
  # Extract base names and wave numbers
  base_names <- sub("(\\d{1,2})$", "", column_names)
  
  # Extract wave numbers
  wave_pattern <- ".*?(\\d{1,2})$"
  has_wave <- grepl("\\d{1,2}$", column_names)
  
  if (!any(has_wave)) {
    stop("No variables with trailing wave numbers found. Variables need trailing digits (1-2) indicating wave.")
  }
  
  # Filter to only variables with wave numbers
  if (!all(has_wave)) {
    warning(sprintf("Excluding %d variable(s) without trailing wave numbers", sum(!has_wave)))
    column_names <- column_names[has_wave]
    base_names <- base_names[has_wave]
  }
  
  wave_numbers <- as.numeric(sub(wave_pattern, "\\1", column_names))
  
  # Check for NA wave numbers
  if (any(is.na(wave_numbers))) {
    bad_vars <- column_names[is.na(wave_numbers)]
    warning(sprintf("Could not extract wave numbers from: %s", paste(head(bad_vars, 3), collapse = ", ")))
    valid_idx <- !is.na(wave_numbers)
    column_names <- column_names[valid_idx]
    base_names <- base_names[valid_idx]
    wave_numbers <- wave_numbers[valid_idx]
  }
  
  # Create grouping dataframe
  var_info <- data.frame(
    original_name = column_names,
    base_name = base_names,
    wave = wave_numbers,
    stringsAsFactors = FALSE
  )
  
  # Group by base name
  var_groups <- split(var_info, var_info$base_name)
  
  # Remove groups with only 1 wave
  var_groups <- var_groups[sapply(var_groups, nrow) > 1]
  
  if (length(var_groups) == 0) {
    stop("No longitudinal variable groups found. Variables need at least 2 waves with trailing digits.")
  }
  
  cat(sprintf("Found %d longitudinal variable groups:\n", length(var_groups)))
  for (i in seq_along(var_groups)) {
    cat(sprintf("  - %s: %d waves\n", names(var_groups)[i], nrow(var_groups[[i]])))
  }
  
  # If ncol is NULL, calculate a reasonable default
  if (is.null(ncol)) {
    ncol <- ceiling(sqrt(min(length(var_groups), max_per_page)))
  }
  
  # Determine number of pages
  n_groups <- length(var_groups)
  n_pages <- ceiling(n_groups / max_per_page)
  
  # Add participant ID column if not present
  n_rows <- nrow(df)
  participant_ids <- 1:n_rows
  
  # Randomly sample participants if max_participants is specified
  if (!is.null(max_participants) && max_participants < n_rows) {
    set.seed(123)  # For reproducibility
    sampled_ids <- sample(participant_ids, max_participants)
    cat(sprintf("Randomly sampling %d of %d participants for visualization\n", 
                max_participants, n_rows))
  } else {
    sampled_ids <- participant_ids
  }
  
  # Create plot for each variable group
  all_plots <- list()
  
  for (base_var in names(var_groups)) {
    
    group_info <- var_groups[[base_var]]
    group_info <- group_info[order(group_info$wave), ]
    
    # Prepare long-format data
    long_data_all <- data.frame()
    
    for (i in 1:nrow(group_info)) {
      var_name <- group_info$original_name[i]
      wave_num <- group_info$wave[i]
      
      temp_df <- data.frame(
        participant_id = participant_ids,
        value = as.numeric(df[[var_name]]),
        wave = wave_num,
        stringsAsFactors = FALSE
      )
      
      long_data_all <- rbind(long_data_all, temp_df)
    }
    
    long_data_all <- long_data_all %>% filter(!is.na(value))
    
    if (nrow(long_data_all) == 0) {
      warning(sprintf("Skipping '%s': no valid data", base_var))
      next
    }
    
    summary_data <- long_data_all %>%
      group_by(wave) %>%
      summarise(
        mean_value = mean(value, na.rm = TRUE),
        n_obs = n(),
        .groups = "drop"
      )
    
    if (nrow(summary_data) < 2) {
      warning(sprintf("Skipping '%s': only %d wave(s) with data", base_var, nrow(summary_data)))
      next
    }
    
    long_data <- long_data_all %>% filter(participant_id %in% sampled_ids)
    
    # Create plot
    p <- ggplot() +
      geom_line(data = long_data, 
                aes(x = wave, y = value, group = participant_id),
                color = "grey80", alpha = 0.5, size = 0.3) +
      geom_point(data = long_data,
                 aes(x = wave, y = value),
                 color = "grey70", alpha = 0.3, size = 0.8) +
      geom_line(data = summary_data, 
                aes(x = wave, y = mean_value),
                color = "black", size = 1.5) +
      geom_point(data = summary_data, 
                 aes(x = wave, y = mean_value),
                 color = "black", size = 3) +
      theme_minimal() +
      theme(
        panel.border = element_rect(color = "black", fill = NA, size = 0.5),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
        plot.margin = margin(5, 5, 5, 5)
      ) +
      scale_x_continuous(breaks = sort(unique(summary_data$wave))) +
      labs(
        title = base_var,
        x = "Wave",
        y = "Value"
      )
    
    all_plots[[base_var]] <- p
  }
  
  if (length(all_plots) == 0) {
    stop("No valid longitudinal plots could be created. Check your data and variable naming.")
  }
  
  n_actual_plots <- length(all_plots)
  n_pages <- ceiling(n_actual_plots / max_per_page)
  
  page_plots <- list()
  for (i in 1:n_pages) {
    start_idx <- (i - 1) * max_per_page + 1
    end_idx <- min(i * max_per_page, n_actual_plots)
    
    plots_for_page <- all_plots[start_idx:end_idx]
    
    page_plot <- wrap_plots(plots_for_page, ncol = ncol) +
      plot_annotation(
        title = sprintf("Longitudinal Change Analysis (Page %d of %d)", i, n_pages),
        theme = theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
      )
    
    page_plots[[i]] <- page_plot
  }
  
  if (n_pages == 1) {
    return(page_plots[[1]])
  } else {
    cat(sprintf("\nCreated %d pages of plots. Access individual pages using [[1]], [[2]], etc.\n", n_pages))
    return(page_plots)
  }
}



## USAGE 
names(df_RED)

# Plot ALL participants (small dataset)
long_plots <- plot_longitudinal_change(
  df = df_RED,
  columns = c(211:270, 277:508), #Your longitudinal variables
  ncol = 3,
  max_per_page = 9,
  max_participants = 350  # change to a number to randomly sample a cleaner plot
)


# View and save
long_plots


# Determine if long_plots is a list or a single plot
if (!is.list(long_plots)) {
  # Single page → save directly
  ggsave(
    filename = "longitudinal_plots_page1.png",
    plot = long_plots,
    width = 12, height = 8, dpi = 300
  )
} else {
  # Multiple pages → loop and save each
  for (i in seq_along(long_plots)) {
    ggsave(
      filename = sprintf("longitudinal_ER_plots_page%d.png", i),
      plot = long_plots[[i]],
      width = 12, height = 8, dpi = 300
    )
  }
}

cat("All plots saved to working directory.\n")


################################################################################
##################### SAVE DATA 

# FIX FAMILY ID 0 TO BE NOT ZERO 
df_RED$FamilyID <- ifelse(df_RED$FamilyID == 0, 11884, df_RED$FamilyID)

FILENAME <- "ABCD_HORM_METH_PREP_12.31.25"

write.csv(df_RED, paste0(FILENAME, ".csv"), row.names=FALSE, na="")

prepareMplusData(df_RED,"ABCD_HORM_METH_PREP_12.31.25.dat", inpfile = TRUE)

names(df_RED)

