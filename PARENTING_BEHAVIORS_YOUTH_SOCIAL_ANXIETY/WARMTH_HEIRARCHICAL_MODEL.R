
################################################################################
############ LOAD PACKAGES

library("metafor")
library("dplyr")
library("robumeta")
library("DescTools")


################################################################################
############ SET WORKING DIRECTORY 

POUNDTOWN <- 0

# set working directory
if (POUNDTOWN == 1) {
  work_dir <- 'C:\\Users\\cjh37695\\Dropbox\\Meta Analysis\\Parenting & Social Anxiety\\Analysis\\PUBLICATION VERSION\\Heirarchical Models\\'
} else {
  work_dir <- 'D:\\Dropbox\\Meta Analysis\\Parenting & Social Anxiety\\Analysis\\PUBLICATION VERSION\\Heirarchical Models\\'
}

setwd(work_dir)


################################################################################
############ WARMTH OVERALL MODELS  

# LOAD DATA FRAME 

df <- read.csv("WARM_HEIR_OVERALL_6.8.24.csv")

# REDUCE TO OUR AGE RANGE

df_red <- subset(df, AGE >= 10 & AGE <20)

# ELIMINATE THE STUDY WITH UNDUE INFLUENCE

WARMTH <- subset(df_red, !(ID %in% c('614_3')))


## INTERCEPT_ONLY MODEL 

WARM_INTERCEPT <- robu(PARENTING ~ 1, data = WARMTH, studynum = ID, #Important: Specify studynum = studyID lets the program 'know' the nesting unit
               var.eff.size = PAR_SE, modelweights = "HIER") #The two two avialable options are "CORR" (correlational) an

print(WARM_INTERCEPT)

# BACK CONVERSION

PARB <- -.195 #main effect size
PARLB <- -.224 #confidence interval lower boundary
PARUB <- -.165 #confidence interval upper boundary 

FisherZInv(PARB) #r = -.19
FisherZInv(PARLB) # CI-LB = -.22  
FisherZInv(PARUB) # CI-UB = -.16

# FOREST PLOT 

forest.robu(WARM_INTERCEPT, es.lab = "PARENTING", study.lab = "ID")

## MODERATION ANALYSES

# PARENT GENDER 

WARM_PARGEN <- robu(PARENTING ~ factor(FATHER), data = WARMTH, studynum = ID,
                          var.eff.size = PAR_SE, modelweights = "HIER") 

print(WARM_PARGEN)

# STUDY DESIGN (current vs retrospective)

WARM_DESIGN <- robu(PARENTING ~ factor(DESIGN), data = WARMTH, studynum = ID, #Important: Specify studynum = studyID lets the program 'know' the nesting unit
                    var.eff.size = PAR_SE, modelweights = "HIER") #The two two avialable options are "CORR" (correlational) an

print(WARM_DESIGN)


# PUBLICATION YEAR 

WARM_YEAR <- robu(PARENTING ~ scale(YEAR, center = 1992), data = WARMTH, studynum = ID, #Important: Specify studynum = studyID lets the program 'know' the nesting unit
                     var.eff.size = PAR_SE, modelweights = "HIER") #The two two avialable options are "CORR" (correlational) an

print(WARM_YEAR)


WARM_YEAR_PAR <- robu(PARENTING ~ scale(YEAR, center = 1992)*factor(FATHER), data = WARMTH, studynum = ID, #Important: Specify studynum = studyID lets the program 'know' the nesting unit
                  var.eff.size = PAR_SE, modelweights = "HIER") #The two two avialable options are "CORR" (correlational) an

print(WARM_YEAR_PAR)


# SAMPLE AGE

WARM_AGE <- robu(PARENTING ~ scale(AGE, center = 10), data = WARMTH, studynum = ID, #Important: Specify studynum = studyID lets the program 'know' the nesting unit
                  var.eff.size = PAR_SE, modelweights = "HIER") #The two two avialable options are "CORR" (correlational) an

print(WARM_AGE)


## COUNTRY REGION 

COUNTRY_MOD <- robu(PARENTING ~ factor(COUNT_REG), data = WARMTH, studynum = ID, #Important: Specify studynum = studyID lets the program 'know' the nesting unit
                    var.eff.size = PAR_SE, modelweights = "HIER") #The two two avialable options are "CORR" (correlational) an

print(COUNTRY_MOD)

################################################################################
############ YOUTH SEX SPLIT MODELS

df2 <- read.csv("WARM_HEIR_SEXSPLIT_6.8.24.csv")

# REDUCE TO OUR SAMPLE AGE 

df_ss_red <- subset(df2, AGE >= 10 & AGE <20)

# REMOVE STUDY WITH UNDUE INFLUENCE

WARMTH_SS <- subset(df_ss_red, !(ID %in% c('614_3')))

## INTERCEPT_ONLY MODEL 

WARM_SS_INTERCEPT <- robu(PARENTING ~ 1, data = WARMTH_SS, studynum = ID, 
                       var.eff.size = PAR_SE, modelweights = "HIER") 

print(WARM_SS_INTERCEPT)

WARM_SS_YOUTHSEX <- robu(PARENTING ~ PRC.MALE, data = WARMTH_SS, studynum = ID,
                          var.eff.size = PAR_SE, modelweights = "HIER") 

print(WARM_SS_YOUTHSEX)

WARM_SS_YOUTHSEX_PAR <- robu(PARENTING ~ PRC.MALE*factor(FATHER), data = WARMTH_SS, studynum = ID, 
                         var.eff.size = PAR_SE, modelweights = "HIER") 

print(WARM_SS_YOUTHSEX_PAR)


