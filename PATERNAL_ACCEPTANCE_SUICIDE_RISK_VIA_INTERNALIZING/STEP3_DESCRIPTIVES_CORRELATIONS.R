# Set working directory 

POUNDTOWN <- 1

# set working directory
if (POUNDTOWN == 1) {
  work_dir <- 'C:\\Users\\cjh37695\\Dropbox\\ABCD_FATHERS\\ABCD_NeuroXInternalizing_Suicide\\ANALYSIS\\'
} else {
  work_dir <- 'C:\\Users\\0910h\\Dropbox\\ABCD_FATHERS\\ABCD_NeuroXInternalizing_Suicide\\ANALYSIS\\'
}

setwd(work_dir)

#Library

library(dplyr)

library(apaTables)
# Load data frame 

df <- read.csv("ABCD_FATHER_X_BRAIN_INT_SUI_FULL_2.10.25.csv")


## Sample 

## SEX
table(df$ysex)
3896/8097

## Mean age 

## T1
mean(df$yage1, na.rm = T)/12
# 9.92
sd(df$yage1, na.rm = T)/12
# .63


## T3
mean(df$yage3, na.rm = T)/12
# 10.9
sd(df$yage3, na.rm = T)/12
# .65

## T5
mean(df$yage5, na.rm = T)/12
# 12.04
sd(df$yage5, na.rm = T)/12
# .67

## RACE/ETHNICTY 

#T1
table(df$yrace)

#white (59%)
4786/8097
# Hispanic (9%)
747/8097
# Latino/a(20%)
1594/8097
#Asian (2%)
164/8097
#Other 
805/8097


## WAVE DATA 

# Create empty lists for the wave data frames
WAVE1 <- WAVE3 <- WAVE5 <- list()

# Keep the ID column in all wave data frames
WAVE1[['subid']] <- df[['subid']]
WAVE3[['subid']] <- df[['subid']]
WAVE5[['subid']] <- df[['subid']]

# Loop through column names and assign numbered variables to the appropriate wave data frames
for (col in colnames(df)) {
  if (grepl("1$", col)) {
    WAVE1[[col]] <- df[[col]]
  }
  if (grepl("3$", col)) {
    WAVE3[[col]] <- df[[col]]
  }
  if (grepl("5$", col)) {
    WAVE5[[col]] <- df[[col]]
  }
}

# Convert the lists to data frames
WAVE1 <- as.data.frame(WAVE1)
WAVE3 <- as.data.frame(WAVE3)
WAVE5 <- as.data.frame(WAVE5)
names(WAVE3)

# Define a function to remove rows where all values from column 1 to the second-to-last column are NA
remove_na_rows <- function(df) {
  # Keep rows where not all values from column 1 to the second-to-last column are NA
  df[rowSums(is.na(df[, 1:(ncol(df) - 1)])) != (ncol(df) - 1), ]
}

# Apply the function to each WAVE data frame
WAVE1 <- remove_na_rows(WAVE1) #N = 8097
WAVE3 <- remove_na_rows(WAVE3) # N = 7796
WAVE5 <- remove_na_rows(WAVE5) # N = 7592


names(df)

## Compute reliability 

## Acceptance
Facc_ALPHA <- alpha(df[ , c("D_1acc1","D_2acc1","D_3acc1","D_4acc1","D_5acc1")])
print(Facc_ALPHA)

Macc_ALPHA <- alpha(df[ , c("M_1acc1","M_2acc1","M_3acc1","M_4acc1","M_5acc1")])
print(Macc_ALPHA)


hist(df$D_1acc1)

names(df)

## INTERNALIZING

INT <- read.csv("ABCD_INTERNALIZING_CBCL_RAW.csv")
  
subID_list <- df$subID

# Filter df2 to only include rows where subID is in subID_list
INT_RED <- INT[INT$subID %in% subID_list, ]


INT_RED_COM <- INT_RED %>%
  filter(complete.cases(CAxDpRP_1,CSomSRP_1,CWtDpRP_1,
                        CAxDpRP_3,CSomSRP_3,CWtDpRP_3))

## RELIABILITY
INT_W1 <- alpha(INT_RED_COM[ , c("CAxDpRP_1","CSomSRP_1","CWtDpRP_1")])
print(INT_W1)

INT_W3 <- alpha(INT_RED_COM[ , c("CAxDpRP_3","CSomSRP_3","CWtDpRP_3")])
print(INT_W3)

INT_W5 <- alpha(INT_RED_COM[ , c("CAxDpRP_5","CSomSRP_5","CWtDpRP_5")])
print(INT_W5)

################################################
########### ASSESS MISSINGNESS##################
################################################

names(df)

PRIM_STU_DF <- df %>%
  select(c(
    #PARENTING
    "pMacc_RS", "pFacc_RS",
    #INTERNALIZING
    "pINTPRT1", "pINTPRT3", 
    #MODERATORS
    "pSNLamyRS", "pSNRamyRS", "pSNLnaccRS", "pSNRnaccRS",
    #OUTCOMES
    "SLF_HRM_CHA", "SUI_ID_CHA"))

sum(colSums(is.na(PRIM_STU_DF)))

##TOTAL NUMBER OF CELLS 
TOTAL <- 10*8097
##NUMER OF MISSING CELLS 
MISSING <- sum(colSums(is.na(PRIM_STU_DF)))
##PERCENT OF MISSING
MISSING/TOTAL

################################################
########### CORRELATION TABLE OF STUDY VARIABLES
################################################

names(df)

CORR_TAB <- df %>%
  select(c("D_acc_M1", "M_acc_M1",
           "C_SALLnac1", "C_SALRnac1", "C_SalLamy1", "C_SalRamy1",
           "C_INTPRT1", "C_INTPRT3", "C_INTPRT5",
            "SuiSTAB" , 
            "SHRMSTAB"))

# Rename varibales for a better plot 

CORR_TAB_NAMED <- CORR_TAB %>%
  rename(Father.Acceptance = D_acc_M1)%>%
  rename(Mother.Acceptance = M_acc_M1)%>%
  rename(Intern.T1 = C_INTPRT1)%>%
  rename(Intern.T3 = C_INTPRT3)%>%
  rename(Intern.T5 = C_INTPRT5)%>%
  rename(`SN-L.NaCC.T1` = C_SALLnac1)%>%
  rename(`SN-R.NaCC.T1` = C_SALRnac1)%>%
  rename(`SN-L.Amy.T1` = C_SalLamy1)%>%
  rename(`SN-R.Amy.T1` = C_SalRamy1)%>%
  rename(Sui.Thoughts.Total = SuiSTAB)%>%
  rename(Self.Harm.Total = SHRMSTAB)

names(CORR_TAB_NAMED)
describe(CORR_TAB_NAMED)

# Make correlations 

COR <- cor(CORR_TAB_NAMED, use = "pairwise.complete.obs")


# Make a table

APA.TABLE <- apa.cor.table(CORR_TAB_NAMED, filename="apa2.doc", table.number=2)

APA.TABLE

