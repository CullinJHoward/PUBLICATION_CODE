################################################################################
######################### LIBRARY

library(metaSEM)
library(lavaan)
library(semPlot)
library(dplyr)


################################################################################
######################### SET WORKING DIRECTORY 

POUNDTOWN <- 1


if (POUNDTOWN == 1) {
  work_dir <- 'C:\\Users\\cjh37695\\Dropbox\\Meta Analysis\\Parenting & Social Anxiety\\Analysis\\PUBLICATION VERSION\\META SEM\\'
} else {
  work_dir <- 'C:\\Users\\0910h\\Dropbox\\Meta Analysis\\Parenting & Social Anxiety\\Analysis\\PUBLICATION VERSION\\META SEM\\'
}

setwd(work_dir)


## LOAD IN DF OF ALL STUDIES 
Jsem_FULL <- read.csv("SEM_WARMandCON_ALLSTUDIES_6.15.24.csv", fileEncoding = "UTF-8" )


## REDUCE TO JUST THE QUALIFYING STUDIES 

incStudies <- c(10, 17, 26, 40, 42, 44, 55, 63, 65, 69, 70, 71, 72, 96, 99, 105,
                130, 148, 150, 152, 194, 233, 442, 509, 514, 532, 551, '614_1', '614_2', 692, 713, 
                771, 773, 775, 777, 781, 784, 8, 12, 24,  36,  37,  43,  88, 511,
                785)

Jsem <- Jsem_FULL %>%
  filter(Jsem_FULL$ID %in% incStudies)

names(Jsem)

################################################################################
######################### DATA ORGANIZATION 

## REORDER THE VARIABLE TO MAKE SURE THEY EFFECT SIZES ALIGN CORRECTLY IN THE 
## CORRELATION MATRIX 

Jsem_REORDER <- Jsem %>%
  select("ID", "AUTHORS", "N", 
         "F_WARM", "M_WARM", "F_CON", "M_CON", "Fwarm_Mwarm", 
         "Fcon_Fwarm", "Mcon_Fwarm", "Fcon_Mwarm", "Mcon_Mwarm", "Fcon_Mcon" )


## MAKE LIST OF CORMATRICES (T.cordat), WITH NA ON THE DIAGONAL

T.nvar <- 5
T.varnames <- c("SocAnx", "FatherWarmth", "MotherWarmth", "FatherControl", "MotherControl")
T.labels <- list(T.varnames,T.varnames)

T.cordat <- list()

for (i in 1:nrow(Jsem_REORDER)){	
  T.cordat[[i]] <- vec2symMat(as.matrix(Jsem_REORDER[i,4:13]),diag = FALSE)
  dimnames(T.cordat[[i]]) <- T.labels
}

# put NA on diagonal if variable is missing

for (i in 1:length(T.cordat)){
  for (j in 1:nrow(T.cordat[[i]])){	
    if (sum(is.na(T.cordat[[i]][j,]))==T.nvar-1) 
    {T.cordat[[i]][j,j] <- NA}
  }}

# put NA on diagonal for variable with least present correlations

for (i in 1:length(T.cordat)){
  for (j in 1:nrow(T.cordat[[i]])){
    for (k in 1:T.nvar){				
      if (is.na(T.cordat[[i]][j,k])==TRUE
          &is.na(T.cordat[[i]][j,j])!=TRUE
          &is.na(T.cordat[[i]][k,k])!=TRUE){
        
        if(sum(is.na(T.cordat[[i]])[j,])>sum(is.na(T.cordat[[i]])[k,]))
        {T.cordat[[i]][k,k] = NA}
        if(sum(is.na(T.cordat[[i]])[j,])<=sum(is.na(T.cordat[[i]])[k,]))
        {T.cordat[[i]][j,j] = NA}	
      }}}}

# TEST FOR NON-POSITIVE DEFINITENESS 
is.pd(T.cordat)

## VIEW IT 
T.cordat


################################################################################
######################### STAGE 1 TSSEM 

#Estimating only the study-level variance (type = Diag)

JointParT1 <- tssem1(Cov=T.cordat, n= Jsem_REORDER$N, method="REM", RE.type="Diag")
summary(JointParT1)


################################################################################
######################### Stage 2 SEM 

#### MAKE THE STUCTURAL PATH MODEL 

model1 <- 
"
    ## Youth social anxiety symptoms (SocAnx) is predicted by Father Warmth (F_WARMTH), Mother Warmth (M_WARMTH)
    ## Father Control (F_CONTROL) and Mother Control (M_CONTROL)
    
    SocAnx ~  FatherWarmth + MotherWarmth + FatherControl + MotherControl

    ## Variances of predictors are fixed at 1
  
    FatherWarmth ~~  1* FatherWarmth
    MotherWarmth ~~  1* MotherWarmth
    FatherControl ~~ 1* FatherControl
    MotherControl ~~ 1* MotherControl

    ## Correlation between all of the parenting predictors 
    FatherWarmth ~~ MotherWarmth
    FatherWarmth ~~ FatherControl 
    FatherWarmth ~~ MotherControl
    MotherWarmth ~~ FatherControl
    MotherWarmth ~~ MotherControl
    FatherControl ~~ MotherControl
    
    ## Error variance of Social Anxiety
    SocAnx ~~ ErrorVar_SOCANX * SocAnx
"

# ENSURE THE MODEL IS CORRECT 

plot(model1)

# CONVERTY MODEL TO RAM FORMAT

RAM1 <- lavaan2RAM(model1)

# RUN STRUCTURAL MODEL WITH DATA 

JointParT2 <- tssem2(JointParT1, RAM=RAM1, diag = TRUE) 

summary(JointParT2)

## VIEW THE PLOT 
plot(JointParT2)
