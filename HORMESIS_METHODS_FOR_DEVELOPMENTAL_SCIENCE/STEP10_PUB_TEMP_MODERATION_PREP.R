## liBRARY 
library(dplyr)

## set environment 
setwd("C:\\Users\\cjh37695\\Dropbox\\HORMESIS_METHODS\\ANALYSIS\\") #work


# participant data 

df <- read.csv("ABCD_HORM_METH_SEM_STRUCTRUAL_4.27.26.csv")


## compute pubertal timing 

df$Ptemp113 <- (df$LMC_PUBlev_11 - df$LMC_PUBlev_3)/(df$LMC_Y_AGE_11-df$LMC_Y_AGE_3)

df$GMC_Ptemp113 <- scale(df$Ptemp113, center = T, scale = F)

## save it 

prepareMplusData(df,"ABCD_HORM_METH_SEM_STRUCTRUAL_MOD_5.25.26.dat", inpfile =T)
write.csv(df, "ABCD_HORM_METH_SEM_STRUCTRUAL_MOD_5.25.26.csv", row.names = F)
