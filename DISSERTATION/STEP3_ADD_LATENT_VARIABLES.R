
### LIBRARY
library(dplyr)
library(ggplot2)
library(tidyr)

POUNDTOWN <- 0

# set working directory
if (POUNDTOWN == 1) {
  work_dir <- 'D:\\DISSERTATION\\ANALYSIS\\STATISTICAL_MODELS\\DORRY_ONLY\\'
} else {
  work_dir <- 'F:\\DISSERTATION\\ANALYSIS\\STATISTICAL_MODELS\\DORRY_ONLY\\'
}

setwd(work_dir)


################################################################################
####################### MERGE IN SAVED LATENT FACTOR SCORES 


################# P-FACTOR PARENT REPORTS 

# PARENT REPORT - LATENT MEANS OF T1 and T2 P-FACTOR 

PAR_PFAC_T1_T2 <- read.table("F:\\DISSERTATION\\ANALYSIS\\STATISTICAL_MODELS\\DORRY_ONLY\\SEM\\MEASUREMENT_MODELS\\P_FACTOR_DORRY.txt", header = FALSE)

# Subset df to the ID and desired latent scores 

PAR_PFAC_T1_T2_RED <- PAR_PFAC_T1_T2[, c("V11", "V13", "V15")]

# Rename the variables
colnames(PAR_PFAC_T1_T2_RED) <- c("Pfac_P1", "Pfac_P2", "ID")


# PARENT REPORT - LATENT CHANGE IN P-FACTOR 

PAR_CHANGE_P_FACTORS_DF <- read.table("F:\\DISSERTATION\\ANALYSIS\\STATISTICAL_MODELS\\DORRY_ONLY\\SEM\\MEASUREMENT_MODELS\\LATENT_CHANGE_P_FACTOR_DORRY.txt", header = FALSE)

# Subset df to the ID and desired latent scores 
names(PAR_CHANGE_P_FACTORS_DF)

PAR_CHANGE_P_FACTORS_DF_RED <- PAR_CHANGE_P_FACTORS_DF[, c("V11", "V13", "V15")]

# Rename the variables
colnames(PAR_CHANGE_P_FACTORS_DF_RED) <- c("Pfac_PLI", "Pfac_PLC", "ID")


################# P-FACTOR YOUTH REPORTS


# YOUTH REPORT - LATENT MEANS OF T1 and T2 P-FACTOR 

Y_PFAC_T1_T2 <- read.table("F:\\DISSERTATION\\ANALYSIS\\STATISTICAL_MODELS\\DORRY_ONLY\\SEM\\MEASUREMENT_MODELS\\YOUTH_P_FACTOR_DORRY.txt", header = FALSE)

# Subset df to the ID and desired latent scores 

Y_PFAC_T1_T2_RED <- Y_PFAC_T1_T2[, c("V11", "V13", "V15")]

# Rename the variables
colnames(Y_PFAC_T1_T2_RED) <- c("Pfac_Y1", "Pfac_Y2", "ID")


# YOUTH REPORT - LATENT CHANGE IN P-FACTOR 

Y_CHANGE_P_FACTORS_DF <- read.table("F:\\DISSERTATION\\ANALYSIS\\STATISTICAL_MODELS\\DORRY_ONLY\\SEM\\MEASUREMENT_MODELS\\YOUTH_LATENT_CHANGE_P_FACTOR_DORRY.txt", header = FALSE)

# Subset df to the ID and desired latent scores 
names(Y_CHANGE_P_FACTORS_DF)

Y_CHANGE_P_FACTORS_DF_RED <- Y_CHANGE_P_FACTORS_DF[, c("V11", "V13", "V15")]

# Rename the variables
colnames(Y_CHANGE_P_FACTORS_DF_RED) <- c("Pfac_YLI", "Pfac_YLC", "ID")


################# SOCIAL UNACCEPTANCE


# YOUTH REPORT - LATENT MEANS OF T1 and T2 SOCIAL UNACCEPTANCE  

Y_SOCACC_T1_T2 <- read.table("F:\\DISSERTATION\\ANALYSIS\\STATISTICAL_MODELS\\DORRY_ONLY\\SEM\\MEASUREMENT_MODELS\\SOCIAL_ACCEPT_DORRY.txt", header = FALSE)

# Subset df to the ID and desired latent scores 

Y_SOCACC_T1_T2_RED <- Y_SOCACC_T1_T2[, c("V9", "V11", "V13")]

# Rename the variables
colnames(Y_SOCACC_T1_T2_RED) <- c("SOACC_1", "SOACC_2", "ID")

# YOUTH REPORT - LATENT CHANGE IN SOCIAL UNACCEPTANCE

Y_CHANGE_SOCACC_DF <- read.table("F:\\DISSERTATION\\ANALYSIS\\STATISTICAL_MODELS\\DORRY_ONLY\\SEM\\MEASUREMENT_MODELS\\LATENT_CHANGE_SOCIAL_ACCEPTANCE_DORRY.txt", header = FALSE)

# Subset df to the ID and desired latent scores 
names(Y_CHANGE_SOCACC_DF)

Y_CHANGE_SOCACC_DF_RED <- Y_CHANGE_SOCACC_DF[, c("V9", "V11", "V13")]

# Rename the variables
colnames(Y_CHANGE_SOCACC_DF_RED) <- c("SOACC_LI", "SOACC_LC", "ID")

################# PERSPECTIVE TAKING & EMPATHIC CONCERN


# YOUTH REPORT - LATENT MEANS OF T1 and T2 PERSPECTIVE TAKING & EMPATHIC CONCERN  

Y_PTEMC_T1_T2 <- read.table("F:\\DISSERTATION\\ANALYSIS\\STATISTICAL_MODELS\\DORRY_ONLY\\SEM\\MEASUREMENT_MODELS\\PERSPECTIVE_TAKING_EMPATHIC_CONCERN_DORRY.txt", header = FALSE)

# Subset df to the ID and desired latent scores 

Y_PTEMC_T1_T2_RED <- Y_PTEMC_T1_T2[, c("V19", "V21", "V23", "V25", "V27")]

# Rename the variables
colnames(Y_PTEMC_T1_T2_RED) <- c("PERTK_1", "EMPCN_1", "PERTK_2", "EMPCN_2", "ID")

# YOUTH REPORT - LATENT CHANGE IN PERSPECTIVE TAKING & EMPATHIC CONCERN

Y_CHANGE_PTEMC_DF <- read.table("F:\\DISSERTATION\\ANALYSIS\\STATISTICAL_MODELS\\DORRY_ONLY\\SEM\\MEASUREMENT_MODELS\\LATENT_CHANGE_PERSPECTIVE_TAKING_EMP_CONCERN_DORRY.txt", header = FALSE)

# Subset df to the ID and desired latent scores 
names(Y_CHANGE_PTEMC_DF)

Y_CHANGE_PTEMC_DF_RED <- Y_CHANGE_PTEMC_DF[, c("V19", "V21", "V23", "V25", "V27")]

# Rename the variables
colnames(Y_CHANGE_PTEMC_DF_RED) <- c("PERTK_LI", "PERTK_LC", "EMPCN_LI", "EMPCN_LC", "ID")


################################################################################
#################### MERGE DATA TOGETHER  

#### LOAD IN ORIGINAL ADVERISTY SCORE DF 

ADV_DF <- read.csv("DISSERTATION_ADVERSITY_DORRY_7.6.25.csv")


JOINT_LAT_VARS <- full_join(PAR_PFAC_T1_T2_RED, Y_PFAC_T1_T2_RED, by = "ID") %>%
  full_join(PAR_CHANGE_P_FACTORS_DF_RED, by = "ID") %>%
  full_join(Y_CHANGE_P_FACTORS_DF_RED, by = "ID") %>%
  full_join(Y_SOCACC_T1_T2_RED, by = "ID") %>%
  full_join(Y_CHANGE_SOCACC_DF_RED, by = "ID") %>%
  full_join(Y_PTEMC_T1_T2_RED, by = "ID") %>%
  full_join(Y_CHANGE_PTEMC_DF_RED, by = "ID")


##### Merge MASTER DF  

FULL_DF <- full_join(ADV_DF, JOINT_LAT_VARS, by = "ID")

##### REDUCE DF TO WANTED VARIABLES 
names(FULL_DF)
RED_FULL_DF <- FULL_DF %>%
  select(c("ID", "P_SEX", "P_EDU", "C_SEX",  "C_AGE",            
          "C_RACE", "MAR_STAT", "PAR_EDU", "PAR_EMPL",
          "INCOME_W1", "SUBSEScm", "SUBSESus",
          "pTHR_Y", "pDEP_Y", "pUNP_Y", "pCUM_Y",
          "pTHR_P", "pDEP_P", "pCUM_P", 
          "C_pTHR_Y", "Q_C_pTHR_Y", "C_pDEP_Y", "Q_C_pDEP_Y", "C_pUNP_Y", "Q_C_pUNP_Y", "C_pCUM_Y", "Q_C_pCUM_Y",
          "C_pTHR_P", "Q_C_pTHR_P", "C_pDEP_P", "Q_C_pDEP_P", "C_pCUM_P", "Q_C_pCUM_P",
          "Pfac_P1", "Pfac_P2", "Pfac_Y1", "Pfac_Y2", "Pfac_PLI", "Pfac_PLC", "Pfac_YLI",         
          "Pfac_YLC", "SOACC_1", "SOACC_2", "SOACC_LI", "SOACC_LC", "PERTK_1", "EMPCN_1",
          "PERTK_2", "EMPCN_2", "PERTK_LI", "PERTK_LC", "EMPCN_LI", "EMPCN_LC"))

##### SAVE FULL AND REDUCED DATA FRAMES 

## FULL FOR DESCRIPTIVES & HEADACHE 

write.csv(FULL_DF, "DISSERTATION_ALL_SURVEY_DATA_DORRY_7.6.25.csv", row.names=FALSE, na="")

## REDUCED FOR EASE OF ANALYSIS

write.csv(RED_FULL_DF, "DISSERTATION_ANALYSIS_SURVEY_DATA_DORRY_7.6.25.csv", row.names=FALSE, na="")

################################################################################
######################## VISUALIZE DISTRIBUTIONS 

## ADVERSITY

names(RED_FULL_DF)
# Convert to long format
df_long <- RED_FULL_DF %>%
  select(ID, pTHR_Y, pDEP_Y, pUNP_Y, pCUM_Y) %>%
  pivot_longer(cols = c(pTHR_Y, pDEP_Y, pUNP_Y, pCUM_Y), 
               names_to = "variable", 
               values_to = "value")

df_long <- RED_FULL_DF %>%
  select(ID, pTHR_Y, pDEP_Y, pUNP_Y, pCUM_Y) %>%
  pivot_longer(
    cols = c(pTHR_Y, pDEP_Y, pUNP_Y, pCUM_Y),
    names_to = "variable",
    values_to = "value",
    names_transform = list(variable = ~ factor(.,
                                               levels = c("pUNP_Y", "pTHR_Y", "pDEP_Y", "pCUM_Y"), # Specify order: top to bottom
                                               labels = c("Unpredictability", "Threat", "Deprivation", "Cumulative") # Custom labels
    ))
  )
# Determine x-axis limits (optional, based on data range)
x_limits <- range(df_long$value, na.rm = TRUE)

# Create stacked histograms
p <- ggplot(df_long, aes(x = value)) +
  geom_histogram(bins = 30, fill = "#1f77b4", color = "black", alpha = 0.8) +
  facet_wrap(~ variable, ncol = 1, scales = "free_y") +
  scale_x_continuous(limits = x_limits, breaks = scales::pretty_breaks(n = 6)) +
  labs(
    #title = "Distribution of Variables",
    x = "Severity of Adversity (Std. 0 - 100)",
    y = "Count"
  ) +
  theme_minimal(base_size = 14, base_family = "serif") +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey80", linetype = "dashed"),
    strip.text = element_text(face = "bold", size = 12),
    strip.background = element_rect(fill = "grey95", color = "black"),
    plot.margin = margin(10, 10, 10, 10)
  )

# Display the plot
print(p)

### SUMMARIZE DISTRIBUTIONS 

#THREAT
skew(RED_FULL_DF$pTHR_Y) #.63
kurtosi(RED_FULL_DF$pTHR_Y) #.67

#DEPRIVATION
skew(RED_FULL_DF$pDEP_Y) #.69
kurtosi(RED_FULL_DF$pDEP_Y) #-.25

#UNPREDICTABILITY
skew(RED_FULL_DF$pUNP_Y) #1.48
kurtosi(RED_FULL_DF$pUNP_Y) #3.56

#CUMULATIVE
skew(RED_FULL_DF$pCUM_Y) #.47
kurtosi(RED_FULL_DF$pCUM_Y) # -.62

# Save as high-resolution image for publication
ggsave("ADVERSITY_DISTRUBUTIONS.png", plot = p, width = 6, height = 8, dpi = 300)

## P-FACTOR DISTRIBUTIONS 

#### YOUTH REPORT 


# Calculate change and categorize with tolerance
RED_FULL_DF <- RED_FULL_DF %>%
  mutate(
    change_youth = Pfac_Y2 - Pfac_Y1,
    change_category = case_when(
      is.na(change_youth) ~ NA_character_,
      change_youth > 1e-6 ~ "Increase",
      TRUE ~ "Decrease"
    ),
    change_category = factor(change_category, levels = c("Increase", "Decrease"))
  )

# Convert to long format for plotting
df_long <- RED_FULL_DF %>%
  select(ID, Pfac_Y1, Pfac_Y2, change_category) %>%
  pivot_longer(
    cols = c(Pfac_Y1, Pfac_Y2),
    names_to = "Timepoint",
    names_pattern = "Pfac_(.*)",
    values_to = "Pfac",
    names_transform = list(
      Timepoint = ~ factor(.,
                           levels = c("Y1", "Y2"),
                           labels = c("T1", "T2")
      )
    )
  ) %>%
  filter(!is.na(Pfac) & !is.na(change_category))  # Remove NA rows


# Create connected line plot for Youth
C_P_FACTOR <- ggplot(df_long, aes(x = Timepoint, y = Pfac, group = ID, color = change_category, fill = change_category)) +
  geom_line(alpha = 0.3, linewidth = 0.5) +
  geom_point(size = 1, alpha = 0.5) +
  scale_color_manual(values = c("Increase" = "#377eb8", "Decrease" = "#ff7f00")) +
  scale_fill_manual(values = c("Increase" = "#377eb8", "Decrease" = "#ff7f00")) +
  scale_x_discrete(expand = c(0.05, 0.05)) +  # Reduce padding around T1 and T2
  labs(
    title = "Youth-Reported P-Factor Scores",
    x = "Timepoint",
    y = "P-Factor Score"
  ) +
  theme_minimal(base_size = 14, base_family = "serif") +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey80", linetype = "dashed"),
    legend.title = element_blank(),
    legend.position = "top",
    legend.text = element_text(size = 12),
    plot.margin = margin(10, 10, 10, 10)
  )

print(C_P_FACTOR)

# Save as high-resolution image for publication
ggsave("YOUTH_REPORT_P_FACTOR_CHANGE.png", plot = C_P_FACTOR, width = 6, height = 4, dpi = 300)


##### PARENT REPORT 

# Calculate change and categorize with tolerance
RED_FULL_DF <- RED_FULL_DF %>%
  mutate(
    change_youth = Pfac_P2 - Pfac_P1,
    change_category = case_when(
      is.na(change_youth) ~ NA_character_,
      change_youth > 1e-6 ~ "Increase",
      TRUE ~ "Decrease"
    ),
    change_category = factor(change_category, levels = c("Increase", "Decrease"))
  )

# Convert to long format for plotting
df_long_PAR <- RED_FULL_DF %>%
  select(ID, Pfac_P1, Pfac_P2, change_category) %>%
  pivot_longer(
    cols = c(Pfac_P1, Pfac_P2),
    names_to = "Timepoint",
    names_pattern = "Pfac_(.*)",
    values_to = "Pfac",
    names_transform = list(
      Timepoint = ~ factor(.,
                           levels = c("P1", "P2"),
                           labels = c("T1", "T2")
      )
    )
  ) %>%
  filter(!is.na(Pfac) & !is.na(change_category))  # Remove NA rows


# Create connected line plot for Youth
P_P_FACTOR <- ggplot(df_long_PAR, aes(x = Timepoint, y = Pfac, group = ID, color = change_category, fill = change_category)) +
  geom_line(alpha = 0.3, linewidth = 0.5) +
  geom_point(size = 1, alpha = 0.5) +
  scale_color_manual(values = c("Increase" = "#a65628", "Decrease" = "#4daf4a")) +
  scale_fill_manual(values = c("Increase" = "#a65628", "Decrease" = "#4daf4a")) +
  scale_x_discrete(expand = c(0.05, 0.05)) +  # Reduce padding around T1 and T2
  labs(
    title = "Parent-Reported P-Factor Scores",
    x = "Timepoint",
    y = "P-Factor Score"
  ) +
  theme_minimal(base_size = 14, base_family = "serif") +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey80", linetype = "dashed"),
    legend.title = element_blank(),
    legend.position = "top",
    legend.text = element_text(size = 12),
    plot.margin = margin(10, 10, 10, 10)
  )

print(P_P_FACTOR)

# Save as high-resolution image for publication
ggsave("PARENT_REPORT_P_FACTOR_CHANGE.png", plot = P_P_FACTOR, width = 6, height = 4, dpi = 300)

