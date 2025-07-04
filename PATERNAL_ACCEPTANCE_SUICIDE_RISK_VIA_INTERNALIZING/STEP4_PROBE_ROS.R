#Library

library(dplyr)
library(MplusAutomation)
library(ggplot2)
library(tidyr)

# Set working directory 


POUNDTOWN <- 0

# set working directory
if (POUNDTOWN == 1) {
  work_dir <- 'C:\\Users\\cjh37695\\Dropbox\\ABCD_FATHERS\\ABCD_NeuroXInternalizing_Suicide\\ANALYSIS\\'
} else {
  work_dir <- 'D:\\Dropbox\\ABCD_FATHERS\\ABCD_NeuroXInternalizing_Suicide\\ANALYSIS\\FINAL_ANALYSIS\\'
}

setwd(work_dir)

# Load data frame 

df <- read.csv("ABCD_FATHER_X_BRAIN_INT_SUI_FULL_2.10.25.csv")

########################## NEW FROM POBIT MODEL 

################################################################################
########################### PROBE ROS - A-Path 
names(df)
# Your real data
plot_df <- data.frame(
  Father_Accept = df$pFacc_RS,
  Internalizing = df$pPint_RS
)

# X = pFacc_RS
#Z = SalLamy1
# Y = pPint_RS


# Randomly sample 500 rows (with complete data)
sampled_df <- plot_df %>%
  filter(!is.na(Father_Accept) & !is.na(Internalizing)) %>%
  filter(Father_Accept > -20 & Father_Accept < 20) %>%
  sample_n(500)

# Create a sequence of X values (same range as PFACC_RS)
x_vals <- seq(-20, 20, length.out = 100)

# My estimates from Mplus 
int_low <- .415
slope_low <- -.124

int_mean <- -.024
slope_mean <- -.033

int_high <- -0.237
slope_high <- 0.012


# Create data frame with predicted Y for each line
lines_df <- data.frame(
  PFACC_RS = x_vals,
  Low  = int_low  + slope_low  * x_vals,
  Mean = int_mean + slope_mean * x_vals,
  High = int_high + slope_high * x_vals
)

lines_long <- pivot_longer(lines_df, cols = c("Low", "Mean", "High"),
                           names_to = "ModeratorLevel", values_to = "PPINT_RS")

# Plot the lines only (you can add your actual data points later)
ggplot() +
 # geom_point(data = sampled_df, aes(x = Father_Accept, y = Internalizing), alpha = 0.4, size = 1) +
  geom_line(data = lines_long, aes(x = PFACC_RS, y = PPINT_RS, color = ModeratorLevel), size = 1.2) +
  labs(
    title = "Simple Slopes from Mplus with Observed Data",
    x = "Father Acceptance (PFACC_RS)",
    y = "Internalizing (PPINT_RS)",
    color = "Moderator Level"
  ) +
  theme_minimal()


#A Path INT on X and INT on Z harmonized RoS
LOWER_ROS <- subset(df, SalLamy1 < .013 & pFacc_RS < -2)

# Z score of father acceptance for the RoS (observed -2)
((-2) - mean(df$pFacc_RS, na.rm = TRUE))/sd(df$pFacc_RS, na.rm = TRUE) # z = -.19 
1078/8097 # 13%

# Conditional RoS
IND_ROS <- subset(df, SalLamy1 < -.06)
1985/8097 # 25%


################################################################################
##################### ORIGNAL FROM LOGISTIC ANALYSIS DO NOT DELETE 

## SUICIDAL IDEATIONS 
names(df)
# FATHER ACCEPTANCE X SAL - L.AMY rsFC
## INDIRECT EFFECT - only one dimension (the moderator)

LOWER_ROS <- subset(df, pSNLamyRS < -2 & pSNLamyRS > -30)
3131/8097 # 39%
UPPER_ROS <- subset(df, pSNLamyRS > 29)
176/8097 # 2%

LOWER_ROS_DAD <- subset(LOWER_ROS, pFacc_RS < -2)
905/8097 #11%
  
HIGHER_ROS <- subset(df, pSNLamyRS > 33)
83/8097 # 1%



################################################################################
#################### EXPONENTIATE LOGIT VALUES #################################
################################################################################

## FATHER ACC -> SUICIDAL IDEATIONS 

exp(-.02) # 0.9801987
1 - 0.9801987 # ~2%

## FATHER -> INTERN -> SUICIDAL IDEATIONS (UNCONDITIONS INDIRECT)

exp (-.001)

# T3 INTERNALIZING ON SUICIDIAL IDEATIONS 

exp(.022) #2% Increase
exp(0.016) 
exp(.027)

#T3 Internalziing on Self harming
exp(0.019) #2% Increase
exp(0.011) 
exp(0.026)


################################################################################
########################### SIMPLE SLOPES PLOT #################################
################################################################################

#ESTABLISH MODEL WEIGHTS 

# Beta weights from your Mplus output (replace with your actual numbers)
b0 <- -0.015        # Example: intercept
b1 <- -0.032     # Effect of X
b2 <- -0.012     # Effect of moderator M
b3 <- 0.003      # Interaction effect


# Define the range for X (e.g., from min to max)
x_seq <- seq(-30, 30, by = 0.5)

# Define specific moderator values (e.g., low, mean, and high)
# You can choose these based on the mean and ±1 standard deviation, or any other criteria
moderator_levels <- c(LBlowROS = -30, UBlowROS = -2, mean = 0, high = 15)  # Replace with your chosen values

## CREATE PREDICTED DATA 

# Create an empty data frame to store predictions
pred_df <- data.frame()

# Loop through each moderator level and calculate predicted Y
for (mod_val in moderator_levels) {
  y_pred <- b0 + b1 * x_seq + b2 * mod_val + b3 * x_seq * mod_val
  temp_df <- data.frame(
    X = x_seq,
    Moderator = as.factor(mod_val),
    Y = y_pred
  )
  pred_df <- rbind(pred_df, temp_df)
}

## PLOT IT 

library(ggplot2)

ggplot(pred_df, aes(x = X, y = Y, color = Moderator)) +
  geom_line(size = 1.2) +
  labs(
    title = "Simple Slopes for the Interaction",
    x = "Predictor (X)",
    y = "Predicted Outcome (Y)",
    color = "Moderator Level"
  ) +
  theme_minimal()
