#Library

library(dplyr)
library(MplusAutomation)
library(ggplot2)

# Set working directory 


POUNDTOWN <- 1

# set working directory
if (POUNDTOWN == 1) {
  work_dir <- 'C:\\Users\\cjh37695\\Dropbox\\ABCD_FATHERS\\ABCD_NeuroXInternalizing_Suicide\\ANALYSIS\\'
} else {
  work_dir <- 'C:\\Users\\0910h\\Dropbox\\ABCD_FATHERS\\ABCD_NeuroXInternalizing_Suicide\\ANALYSIS\\'
}

setwd(work_dir)

# Load data frame 

df <- read.csv("ABCD_FATHER_X_BRAIN_INT_SUI_FULL_2.6.25.csv")

################################################################################
#################### NEW CORRELATIONS ##########################################
################################################################################


names(df)

CORR_TAB <- df %>%
  select(c("income", "Par.Edu",
           "yage1", "yage5",
           "D_acc_M1", "M_acc_M1",
           "INTPRT1", "INTPRT3",
           "SalSal1", 
           "SalLamy1", "SalRamy1",
           "SALLnac1",  "SALRnac1",
           "SUI_ID_CHA",
           "SLF_HRM_CHA"))

# Rename varibales for a better plot 

CORR_TAB_NAMED <- CORR_TAB %>%
  rename(Father.Acceptance = D_acc_M1)%>%
  rename(Mother.Acceptance = M_acc_M1)%>%
  rename(Intern.T1 = INTPRT1)%>%
  rename(InternT3 = INTPRT3)%>%
  rename(Age.T1 = yage1)%>%
  rename(Age.T5 = yage5)%>%
  rename(Income = income)%>%
  rename(Sui.Thoughts.T5 = SUI_ID_CHA)%>%
  rename(Self.Harm.T5 = SLF_HRM_CHA)

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

## PROPORTIONS OF 0's IN DUMMY CODES 

table(CORR_TAB_NAMED$Sui.Thoughts.T5)
#1014/8097 = 13%

table(CORR_TAB_NAMED$Self.Harm.T5)
#535/8097 = .07%


################################################
########### NEW T-TESTS AND ANOVAS #############
################################################

#### GENDER DIFFERENCES 
#Ensure 2 levels 
table(df$ysex)

df$ysex <- ifelse(df$ysex == 1, 1,
                  ifelse(df$ysex == 2, 2,
                         ifelse(df$ysex == 3, NA,NA)))

## BIOLOGICAL SEX 
df$ysex <- as.factor(df$ysex)

SEX_INTERN1 <- t.test(INTPRT1 ~ ysex, data = df, var.equal = TRUE)  # Assume equal variances
SEX_INTERN3 <- t.test(INTPRT3 ~ ysex, data = df, var.equal = TRUE)  # Assume equal variances
SEX_FACC <- t.test(D_acc_M1 ~ ysex, data = df, var.equal = TRUE)  # Assume equal variances
SEX_MACC <- t.test(M_acc_M1 ~ ysex, data = df, var.equal = TRUE)  # Assume equal variances


# Print results
print(SEX_INTERN1)
print(SEX_INTERN3)
print(SEX_FACC)
print(SEX_MACC)


### CHI SQUARE TESTS ###
SEX_SUI_ID_CON_TAB <- table(df$ysex, df$SUI_ID_CHA)

# Perform Chi-square test
SEX_SUID_CHI_result <- chisq.test(SEX_SUI_ID_CON_TAB)

# Display the result
SEX_SUID_CHI_result

SEX_SLF_HRM_CON_TAB <- table(df$ysex, df$SLF_HRM_CHA)

# Perform Chi-square test
SEX_SLF_HRM_CHI_result <- chisq.test(SEX_SLF_HRM_CON_TAB)

# Display the result
SEX_SLF_HRM_CHI_result

## RACE/ETHNICITY
df$yrace <- as.factor(df$yrace)
# Perform one-way ANOVA
ANOVA_Facc <- aov(D_acc_M1 ~ yrace, data = df)
ANOVA_Macc <- aov(M_acc_M1 ~ yrace, data = df)
ANOVA_INTERN1 <- aov(INTPRT1 ~ yrace, data = df)
ANOVA_INTERN3 <- aov(INTPRT3 ~ yrace, data = df)

# POST HOC
ANOVA_Macc_PH <- TukeyHSD(ANOVA_Macc)
print(ANOVA_Macc_PH)

summary(ANOVA_INTERN1)
summary(ANOVA_INTERN3)
# POST HOC


### CHI SQUARE TESTS ###
RACE_SUI_ID_CON_TAB <- table(df$yrace, df$SUI_ID_CHA)

# Perform Chi-square test
RACE_SUID_CHI_result <- chisq.test(RACE_SUI_ID_CON_TAB)

# Display the result
RACE_SUID_CHI_result

RACE_SLF_HRM_CON_TAB <- table(df$yrace, df$SLF_HRM_CHA)

# Perform Chi-square test
RACE_SLF_HRM_CHI_result <- chisq.test(RACE_SLF_HRM_CON_TAB)

# Display the result
RACE_SLF_HRM_CHI_result
