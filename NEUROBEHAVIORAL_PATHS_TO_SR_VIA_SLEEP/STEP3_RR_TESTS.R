#Library
library(dplyr)
library(MplusAutomation)

# Set working directory 


POUNDTOWN <- 1

# set working directory
if (POUNDTOWN == 1) {
  work_dir <- 'C:\\Users\\cjh37695\\OneDrive - University of Georgia\\Suicide sleep paper\\ANALYSIS\\ANALYSIS_R&R\\'
} else {
  work_dir <- 'C:\\Users\\0910h\\OneDrive - University of Georgia\\Suicide sleep paper\\ANALYSIS\\ANALYSIS_R&R\\\\'
}

setwd(work_dir)


## LOAD DF 

df <- read.csv("ABCD_SLP_DMN_SUI_RR_6.25.25.csv")

################################################################################
################## TEST SCANNER INFORMATION
hist(df$dmn_dmn_W1)
names(df)
table(df$INC_GRP)

hist(df$LAG_1_5)
## DMN rsFC by scanner model 

SCANmod <- aov(dmn_dmn_W1 ~ SCANman_1, data = df)
summary(SCANmod)

## Testing incidental findings on all realted variables 
names(df)
DMN <- t.test(dmn_dmn_W1 ~ INC_GRP, data = df, var.equal = TRUE)
DMN
SLPprob <- t.test(SLPROrs_L ~ INC_GRP, data = df, var.equal = TRUE)
SLPprob
SLPdur <- t.test(PoSLpDrs_L ~ INC_GRP, data = df, var.equal = TRUE)
SLPdur
SA <- t.test(SA_3_5 ~ INC_GRP, data = df, var.equal = TRUE)
SA
SUD <- t.test(sud_3_5 ~ INC_GRP, data = df, var.equal = TRUE)
SUD

# SA contingency Table
SA_cont_table <- table(df$INC_GRP, df$SA_3_5)
# Run chi-square test
chi_result <- chisq.test(SA_cont_table)
print(chi_result)
# SUD contingency Table
SUD_cont_table <- table(df$INC_GRP, df$sud_3_5)
# Run chi-square test
chi_result <- chisq.test(SUD_cont_table)
print(chi_result)

### TEST RACE/ETHNICITY BIAS 

FULL_DF <- read.csv("C:\\Users\\cjh37695\\OneDrive - University of Georgia\\Suicide sleep paper\\ANALYSIS\\ANALYSIS_R&R\\ABCD_SES_DMN_SUI_6.11.25.csv")

## ID THOSE IN INCLUDED SAMPLE 
id_list <- df$subID
## MAKE VARIABLE IN FULL df indexing those 
FULL_DF$in_df <- as.integer(FULL_DF$subID %in% id_list)
table(FULL_DF$in_df)

# RACE/ETHNICITY contingency Table
RACE_EX_cont_table <- table(FULL_DF$in_df, FULL_DF$Y_RACE)
# Run chi-square test
chi_result <- chisq.test(RACE_EX_cont_table)
print(chi_result)

# Step 1: View the contingency table
print(RACE_EX_cont_table)

# Step 2: View observed vs. expected frequencies
print(chi_result$expected)

# Step 3: Calculate standardized residuals
print(chi_result$stdres)

# Step 6: Visualize with a mosaic plot
mosaicplot(RACE_EX_cont_table, main = "Mosaic Plot of in_df by Y_RACE", 
           color = TRUE, xlab = "in_df", ylab = "Y_RACE")

