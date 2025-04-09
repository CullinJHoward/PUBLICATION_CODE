################################################################################
######################### Library

library(metaSEM)
library(lavaan)
library(semPlot)
library(dplyr)

################################################################################
######################### Set working Directory 

POUNDTOWN <- 0

# set working directory
if (POUNDTOWN == 1) {
  work_dir <- 'C:\\Users\\cjh37695\\Dropbox\\Meta Analysis\\Parenting & Social Anxiety\\Analysis\\PUBLICATION VERSION\\META SEM\\'
} else {
  work_dir <- 'D:\\Dropbox\\Meta Analysis\\Parenting & Social Anxiety\\Analysis\\PUBLICATION VERSION\\META SEM\\'
}

setwd(work_dir)


# LOAD MY FULL DF 

df <- read.csv("SEM_WARMandCON_ALLSTUDIES_6.15.24.csv", fileEncoding = "UTF-8" )

# REDUCE TO ONLY THE INCLUDED STUDIES (BY ID)

incStudies <- c(10, 17, 26, 40, 42, 44, 55, 63, 65, 69, 70, 71, 72, 96, 99, 105,
130, 148, 150, 152, 194, 233, 442, 509, 514, 532, 551, '614_1', '614_2', 692, 713, 
771, 773, 775, 777, 781, 784)

WARMTH <- df %>%
  filter(df$ID %in% incStudies)

################################################################################
######################### DATA ORGANIZATION 

## MAKE LIST OF CORMATRICES (T.cordat), WITH NA ON THE DIAGONAL

T.nvar <- 3
T.varnames <- c("FATH_CR","MOTH_CR","FA_MA_CR")
T.labels <- list(T.varnames,T.varnames)

T.cordat <- list()

for (i in 1:nrow(WARMTH)){	
  T.cordat[[i]] <- vec2symMat(as.matrix(WARMTH[i,13:15]),diag = FALSE)
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

T.stage1random <- tssem1(Cov=T.cordat, n= WARMTH$N, method="REM", RE.type="Diag")
summary(T.stage1random)


################################################################################
######################### Stage 2 SEM 

#### MAKE THE STUCTURAL PATH MODEL 

model1 <- "
    ## Social Anxiety (SOCANX) is modeled by Mother Control (Mother_CON) and Father Control (Father_CON)
    SOCANX ~ MOTH_CR * Mother_CON + FATH_CR * Father_CON

    ## Variances of predictors are fixed at 1
    Mother_CON ~~ 1 * Mother_CON
    Father_CON ~~ 1 * Father_CON

    ## Correlation between the predictors (Mother Control and Father Control)
    Mother_CON ~~ FA_MA_CR * Father_CON

    ## Error variance of Social Anxiety
    SOCANX ~~ ErrorVar_SOCANX * SOCANX
"

# ENSURE THE MODEL IS CORRECT 

plot(model1)

# CONVERTY MODEL TO RAM FORMAT

RAM1 <- lavaan2RAM(model1)

# RUN STRUCTURAL MODEL WITH DATA 

ParCON2 <- tssem2(T.stage1random, RAM=RAM1) 
summary(ParCON2)

## VIEW THE PLOT 
plot(ParCONT2)