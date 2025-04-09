################################################################################
######################### LIBRARY

library(metaSEM)
library(lavaan)
library(semPlot)
library(dplyr)


################################################################################
######################### SET WORKING DIRECTORY

POUNDTOWN <- 0

if (POUNDTOWN == 1) {
  work_dir <- 'C:\\Users\\cjh37695\\Dropbox\\Meta Analysis\\Parenting & Social Anxiety\\Analysis\\PUBLICATION VERSION\\META SEM\\'
} else {
  work_dir <- 'D:\\Dropbox\\Meta Analysis\\Parenting & Social Anxiety\\Analysis\\PUBLICATION VERSION\\META SEM\\'
}

setwd(work_dir)

# LOAD FULL DF 

df <- read.csv("SEM_WARMandCON_ALLSTUDIES_6.15.24.csv", fileEncoding = "UTF-8" )


# REDUCE TO ONLY THE INCLUDED STUDIES (BY ID)

incStudies <- c(8,  12,  17,  24,  26,  36,  37,  43,  55,  69,
                71,  88,  96,  99, 105, 148, 150, 442, 511, 532,
                551, 771, 773, 777, 781, 784, 785)

CONTROL <- df %>%
  filter(df$ID %in% incStudies)

## REORDER THE VARIABLES TO MATCH MY SEM FORMAT  

names(CONTROL_SEM)

CONTROL_SEM <- CONTROL %>%
  select("ID", "AUTHORS", "N", 
         "F_CON",  "M_CON", "Fcon_Mcon" )


################################################################################
######################### Data Organization 


## MAKE LIST OF CORMATRICES (T.cordat), WITH NA ON THE DIAGONAL

T.nvar <- 3
T.varnames <- c("FATHER", "MOTHER", "SOCANX")
T.labels <- list(T.varnames,T.varnames)

T.cordat <- list()

for (i in 1:nrow(CONTROL_SEM)){	
  T.cordat[[i]] <- vec2symMat(as.matrix(CONTROL_SEM[i,4:6]),diag = FALSE)
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

head(T.cordat)

# TEST FOR NON-POSITIVE DEFINITENESS  
is.pd(T.cordat)

## VIEW IT 
T.cordat


################################################################################
######################### STAGE 1 TSSEM 

ParCONT1 <- tssem1(Cov=T.cordat, n= CONTROL_SEM$N, method="REM", RE.type="Diag")
summary(ParCONT1)


################################################################################
######################### Stage 2 SEM 

#### MAKE THE STUCTURAL PATH MODEL 


model1 <- "
    ## Social Anxiety (SOCANX) is modeled by Mother Control (Mother_CON) and Father Control (Father_CON)
    SOCANX ~  FATHER + MOTHER
  
    
    ## Variances of predictors are fixed at 1
    FATHER ~~ 1 * FATHER
    MOTHER ~~ 1 * MOTHER
    

    ## Correlation between the predictors (Mother Control and Father Control)
    MOTHER ~~ FATHER

    ## Error variance of Social Anxiety
    SOCANX ~~  SOCANX
"

# ENSURE THE MODEL IS CORRECT 

plot(model1)

# CONVERTY MODEL TO RAM FORMAT

RAM1 <- lavaan2RAM(model1)

# RUN STRUCTURAL MODEL WITH DATA 

ParCONT2 <- tssem2(ParCONT1, RAM=RAM1) 
summary(ParCONT2)

## VIEW THE PLOT 
plot(ParCONT2)

