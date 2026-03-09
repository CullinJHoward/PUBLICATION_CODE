#Library
library(dplyr)
library(ggplot2)
library("psych")
library("MplusAutomation")
# Set working directory 


POUNDTOWN <- 0

# set working directory
if (POUNDTOWN == 1) {
  work_dir <- 'C:\\Users\\cjh37695\\OneDrive - University of Georgia\\Suicide sleep paper\\ANALYSIS\\ANALYSIS_R&R\\'
} else {
  work_dir <- 'C:\\Users\\0910h\\OneDrive - University of Georgia\\Suicide sleep paper\\ANALYSIS\\ANALYSIS_R&R\\'
}

setwd(work_dir)

df <- read.csv("C:\\Users\\cjh37695\\OneDrive - University of Georgia\\Suicide sleep paper\\ANALYSIS\\ANALYSIS_R&R\\ABCD_SLP_DMN_SUI_RR_7.1.25.csv")
names(df)
table(df$TRA_HISb_1)
## LOAD DF 

df <- read.csv("PICKAPOINT_PLOT.csv")


# Convert DMN to a factor with specified order
df$DMN <- factor(df$DMN, levels = c("-2 SD", "Mean", "+2 SD"))

ggplot(df, aes(x = SES, y = EST, color = DMN, group = DMN)) +
  geom_line(linewidth = 1) +              # Thicker lines for visibility
  geom_point(size = 3.5, shape = 19) +    # Larger, filled circular points
  scale_x_continuous(breaks = c(-2, -1, 0, 1, 2), 
                     labels = c("-2", "-1", "0", "1", "2")) +  # Explicit x-axis labels
  scale_color_brewer(palette = "Set2",    # Professional color palette
                     name = "DMN rsFC") + # Legend title
  labs(x = "SES-H (SD)",  # Custom x-axis label
       y = "Mediated Probability of Suicidal Ideations 
       via Poor Sleep Duration") +     # Custom y-axis label
  theme_minimal(base_size = 14) +         # Clean theme with larger base font size
  theme(
    text = element_text(family = "serif"), # Use serif font (similar to Times New Roman)
    axis.title.x = element_text(face = "bold", size = 16, margin = margin(t = 10)), # Bold x-axis label
    axis.title.y = element_text(face = "bold", size = 16, margin = margin(r = 10)), # Bold y-axis label
    axis.text = element_text(size = 12, color = "black"), # Clear axis text
    legend.title = element_text(face = "bold", size = 14), # Bold legend title
    legend.text = element_text(size = 12),                # Clear legend text
    legend.position = "bottom",                            # Legend on right
    panel.grid.major = element_line(color = "grey80", linewidth = 0.5), # Subtle grid
    panel.grid.minor = element_blank(),                   # Remove minor grid
    plot.title = element_blank()                          # No overall title
  )

# Save as PNG (300 DPI for publication)
ggsave("CONDITION_SIMPLE_SLOPES.png", width = 7, height = 5, dpi = 300)


################################################################################
###################### JN PLOT REGION OF SIGNIFICANCE 


df <- read.csv("ABCD_SLP_DMN_SUI_RR_7.1.25.csv")

hist(df$C_pC_dmn_dmn_W1)

RoS <- subset(df, C_pC_dmn_dmn_W1 < .18 & C_pC_dmn_dmn_W1 > - .18)

7636/8061 # 94%, n = 7636

## Descriptives 

names(df)

mean(df$Y_AGE_1, na.rm=TRUE)/12 #9.94 (.63)
sd(df$Y_AGE_1, na.rm=TRUE)/12
mean(df$Y_AGE_3, na.rm=TRUE)/12 #10.95 (.65)
sd(df$Y_AGE_3, na.rm=TRUE)/12
mean(df$Y_AGE_5, na.rm=TRUE)/12 #12.05 (.67)
sd(df$Y_AGE_5, na.rm=TRUE)/12

################################################################################
#################### SUBSAMPLE TO NONCLINICAL REFERRALS 
table(df$INC_GRP)
NONCLIN <- subset(df, INC_FND_1 != 4)
table(df$INC_FND_5)

write.csv(NONCLIN, "ABCD_SLP_DMN_SUI_RR_7.14.25_NONCLIN.csv", row.names=FALSE, na="")

prepareMplusData(NONCLIN,"ABCD_SLP_DMN_SUI_RR_7.14.25_NONCLIN.dat")
