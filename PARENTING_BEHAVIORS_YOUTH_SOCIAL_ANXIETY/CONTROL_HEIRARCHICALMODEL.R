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
############ CONTROL OVERALL MODELS  


# LOAD DATA FRAME 

df <- read.csv("CONTROL_HEIR_OVERALL_6.8.24.csv",fileEncoding = "UTF-8")

# LIMIT STUDIES TO ONLY OUR AGR RANGE 

CONTROL <- subset(df, AGE >= 10 & AGE <20)

################################################################################
############ RUN MODELS 


## INTERCEPT ONLY

CON_INTERCEPT <- robu(PARENTING ~ 1, data = CONTROL, studynum = ID, #Important: Specify studynum = studyID lets the program 'know' the nesting unit
               var.eff.size = PAR_SE, modelweights = "HIER")

print(CON_INTERCEPT)


# BACK CONVERSION 

PARB <- .20 #main effect size
PARLB <- .165 #confidence interval lower boundary
PARUB <- .235 #confidence interval upper boundary 

FisherZInv(PARB) #r = .20
FisherZInv(PARLB) # CI-LB = .16  
FisherZInv(PARUB) # CI-UB = .23

# FOREST PLOT 

forest.robu(CON_INTERCEPT, es.lab = "PARENTING", study.lab = "ID")


## MODERATION ANALYSES  

# PARENT GENDER

CON_PARGEN <- robu(PARENTING ~ factor(FATHER), data = CONTROL, studynum = ID, 
                          var.eff.size = PAR_SE, modelweights = "HIER") 

print(CON_PARGEN)


# STUDY DESIGN (current vs retrospective)

CON_DESIGN <- robu(PARENTING ~ factor(DESIGN), data = CONTROL, studynum = ID, 
                    var.eff.size = PAR_SE, modelweights = "HIER") 

print(CON_DESIGN)
table(CONTROL$DESIGN)


######  PUBLICATION YEAR 

CONTROL$MOM <- ifelse(CONTROL$FATHER == 0, 1, 0)

CON_YEAR <- robu(PARENTING ~ scale(YEAR, center = 1992), data = CONTROL, studynum = ID, 
                     var.eff.size = PAR_SE, modelweights = "HIER") 

print(CON_YEAR)


CON_YEAR_PAR <- robu(PARENTING ~ scale(YEAR, center = 1992)*factor(FATHER), data = CONTROL, studynum = ID, 
                  var.eff.size = PAR_SE, modelweights = "HIER") 

print(CON_YEAR_PAR)


# SAMPLE AGE 

CON_AGE <- robu(PARENTING ~ scale(AGE, center = 10), data = CONTROL, studynum = ID, 
                  var.eff.size = PAR_SE, modelweights = "HIER") 

print(CON_AGE)

CON_AGE_PAR <- robu(PARENTING ~ scale(AGE, center = 10)*factor(FATHER), data = CONTROL, studynum = ID, 
                 var.eff.size = PAR_SE, modelweights = "HIER") 

print(CON_AGE_PAR)


## COUNTRY REGION 

COUNTRY_MOD <- robu(PARENTING ~ factor(COUNT_REG), data = CONTROL, studynum = ID, 
                    var.eff.size = PAR_SE, modelweights = "HIER") 

print(COUNTRY_MOD)


################################################################################
############ YOUTH SEX SPLIT ANALYSES 

df2 <- read.csv("CONTROL_HEIR_SEXSPLIT_6.8.24.csv")

# REDUCE TO OUR AGE RANGE 
CONTROL_SS <- subset(df, AGE >= 10 & AGE <20)

## INTERCEPT_ONLY MODEL 


CON_SS_INTERCEPT <- robu(PARENTING ~ 1, data = CONTROL_SS, studynum = ID, 
                       var.eff.size = PAR_SE, modelweights = "HIER") 

print(CON_SS_INTERCEPT)


## MODERATORS

# CHILD GENDER 

CON_SS_YOUTHSEX <- robu(PARENTING ~ PRC.MALE, data = CONTROL_SS, studynum = ID, 
                          var.eff.size = PAR_SE, modelweights = "HIER") 

print(CON_SS_YOUTHSEX)


