## IDENTIFY ROS FOR MODERATION MODEL 


## IDENTIFY PEOPLE IN ROS 
df <- read.csv("C:\\Users\\cjh37695\\Dropbox\\ABCD_FATHERS\\ABCD_NeuroXInternalizing_Suicide\\ANALYSIS\\RR1\\FINAL_MODELS\\SAPPENFIELD_STRUCTURAL_LAST_T1T9_V3.csv")
names(df)

ROS <- subset(df, W_fSALNACt1t9_C > -.80)

X <- -.79
mean <- mean(ROS$W_fSALNACt1t9_C, na.rm = T)
sd <- sd(ROS$W_fSALNACt1t9_C, na.rm = T)

z <- (X - mean) / sd
z
3292/4246