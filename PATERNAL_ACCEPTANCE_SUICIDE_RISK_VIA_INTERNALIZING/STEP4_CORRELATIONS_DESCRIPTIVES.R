#Library

library(dplyr)
library(psych)
library(stringr)
library(purrr)
library(MplusAutomation)

# Set working directory 

setwd("/home/cjh37695/ABCD_PROJECTS/SAPPENFIELD_FATHERS/")

#Library

df <- read.csv("SAPPENFIELD_STRUCTURAL_LAST_FULL_V3.csv")


################################################################################
################### CORRELATION TABLE 
names(df)

CORR_TAB <- df %>%
  select(c("INCOME6L_1", 
           "Y_AGE_1", 
           "D_acc_M_1", 
           "M_acc_M_1",
           "CIntSTP_1", 
           "mS_SALNAC1", 
           "mS_SALAMY1"))

# Rename varibales for a better plot 

CORR_TAB_NAMED <- CORR_TAB %>%
  rename(Income = INCOME6L_1)%>%
  rename(Age = Y_AGE_1)%>%
  rename(F_Accept = D_acc_M_1)%>%
  rename(M_Accept = M_acc_M_1)%>%
  rename(Intern = CIntSTP_1)%>%
  rename(Sal_Nacc = mS_SALNAC1)%>%
  rename(Sal_Amy = mS_SALAMY1)

names(CORR_TAB_NAMED)


# MAKE THE TOTAL VARIABLES DUMM

# Compute the correlation matrix with significance
corr_results <- Hmisc::rcorr(as.matrix(CORR_TAB_NAMED))

# Extract correlations, p-values, means, and standard deviations
correlations <- round(corr_results$r, 2)  # Correlation coefficients
p_values <- corr_results$P  # p-values
means <- round(colMeans(CORR_TAB_NAMED, na.rm = TRUE), 2)  # Means
sds <- round(apply(CORR_TAB_NAMED, 2, sd, na.rm = TRUE), 2)  # Standard deviations

# Add significance stars to the correlation matrix
stars <- ifelse(p_values < 0.001, "***",
                ifelse(p_values < 0.01, "**",
                       ifelse(p_values < 0.05, "*", "")))
cor_with_stars <- matrix(paste(correlations, stars, sep = ""), 
                         nrow = nrow(correlations), 
                         dimnames = dimnames(correlations))

# Convert to a data frame
cor_table <- as.data.frame(cor_with_stars)
cor_table <- cbind(Variable = rownames(cor_table), cor_table)

# Create a summary statistics table
summary_stats <- data.frame(
  Variable = colnames(CORR_TAB_NAMED),
  Mean = means,
  SD = sds
)

# Print the summary statistics
cat("Summary Statistics:\n")
print(summary_stats)

# Print the correlation table
cat("\nCorrelation Table:\n")
print(cor_table, row.names = FALSE)

write.csv(cor_table, "RAW_T1_CORRELATIIONS.csv", row.names = F)

################################################################################
################# T-TESTS 

# Acceptance
library(psych)
describeBy(df$CIntSTP_1,df$Y_SEX)
## GENDER 
t.test(df$D_acc_M_1 ~ df$Y_SEX) #t(8470) = -5.29, p < .001
t.test(df$M_acc_M_1 ~ df$Y_SEX) #t(9500) = -4.85, p < .001
## SAMPLE 
t.test(df$D_acc_M_1 ~ df$IMG_PT1T9) #t(7893) = -1.79, p < .07
t.test(df$M_acc_M_1 ~ df$IMG_PT1T9) #t(8965) = -1.25, p < .21

# Internalizing 
names(df)
t.test(df$CIntSTP_1 ~ df$Y_SEX) #t(9991) = 9.06, p < .001

## SAMPLE 
t.test(df$CIntSTP_1 ~ df$IMG_PT1T9) #t(9274) = .70, p < .49

# Suicidal ideations 
# SEX
SI_tab_SX <- table(df$Y_SEX, df$SV_SuIdPasT_1)
SI_tab_SX_chi <- chisq.test(SI_tab_SX)
SI_tab_SX_N <- sum(SI_tab_SX)
SI_tab_SX_phi <- sqrt(SI_tab_SX_chi$statistic / SI_tab_SX_N)
# print results
cat("Chi-square:", round(SI_tab_SX_chi$statistic, 3), "\n")
cat("df:", SI_tab_SX_chi$parameter, "\n")
cat("p:", round(SI_tab_SX_chi$p.value, 4), "\n")
cat("Phi:", round(SI_tab_SX_phi, 3), "\n")

#SAMPLE 
SI_tab_SM <- table(df$IMG_PT1T9, df$SV_SuIdPasT_1)
SI_tab_SM_chi <- chisq.test(SI_tab_SM)
SI_tab_SM_N <- sum(SI_tab_SM)
SI_tab_SM_phi <- sqrt(SI_tab_SM_chi$statistic / SI_tab_SM_N)
# print results
cat("Chi-square:", round(SI_tab_SM_chi$statistic, 3), "\n")
cat("df:", SI_tab_SM_chi$parameter, "\n")
cat("p:", round(SI_tab_SM_chi$p.value, 4), "\n")
cat("Phi:", round(SI_tab_SM_phi, 3), "\n")


# NSSI 
NSSI_tab_SX <- chisq.test(table(df$Y_SEX, df$SV_NSSIT_1))
print(NSSI_tab_SX)
NSSI_tab_SM <- chisq.test(table(df$IMG_PT1T9, df$SV_NSSIT_1))
print(NSSI_tab_SM)

# SEX
NSSI_tab_SX <- table(df$Y_SEX, df$SV_NSSIT_1)
NSSI_tab_SX_chi <- chisq.test(NSSI_tab_SX)
NSSI_tab_SX_N <- sum(NSSI_tab_SX)
NSSI_tab_SX_phi <- sqrt(NSSI_tab_SX_chi$statistic / NSSI_tab_SX_N)
# print results
cat("Chi-square:", round(NSSI_tab_SX_chi$statistic, 3), "\n")
cat("df:", NSSI_tab_SX_chi$parameter, "\n")
cat("p:", round(NSSI_tab_SX_chi$p.value, 4), "\n")
cat("Phi:", round(NSSI_tab_SX_phi, 3), "\n")

#SAMPLE 
NSSI_tab_SM <- table(df$IMG_PT1T9, df$SV_NSSIT_1)
NSSI_tab_SM_chi <- chisq.test(NSSI_tab_SM)
NSSI_tab_SM_N <- sum(NSSI_tab_SM)
NSSI_tab_SM_phi <- sqrt(NSSI_tab_SM_chi$statistic / NSSI_tab_SM_N)
# print results
cat("Chi-square:", round(NSSI_tab_SM_chi$statistic, 3), "\n")
cat("df:", NSSI_tab_SM_chi$parameter, "\n")
cat("p:", round(NSSI_tab_SM_chi$p.value, 4), "\n")
cat("Phi:", round(NSSI_tab_SM_phi, 3), "\n")

