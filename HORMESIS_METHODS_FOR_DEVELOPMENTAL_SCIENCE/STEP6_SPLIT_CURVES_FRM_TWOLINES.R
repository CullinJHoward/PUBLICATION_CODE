#Library

library(dplyr)
library(psych)
library(stringr)
library(tidyr)
library(purrr)
library(MplusAutomation)
library(tidyverse)
library(rlang)

# Set working directory 

setwd("/home/cjh37695/ABCD_PROJECTS/HORMESIS_METHODS/")

# LOAD DF 

df <- read.csv("ABCD_HORM_METH_SEM_STRUCTRUAL_4.24.26.csv")

################################################################################
################## CREATE SLOPE VARIABLES FOR TESTING HORMESIS

## create a fucntion for splitting adversity variables

bendItLikeHorm <- function(data,
                             base_name,
                             adv_var,
                             horm_vrt,
                             horm_vrt_sd = NULL,
                             sd_offsets = c(-1, 1)) {
  
  adv_var_quo <- enquo(adv_var)
  adv_vals    <- eval_tidy(adv_var_quo, data)
  
  make_pair <- function(vrt_value, suffix = "") {
    xc <- adv_vals - vrt_value
    
    s1_name <- paste0("S1_", suffix, base_name)
    s2_name <- paste0("S2_", suffix, base_name)
    
    data[[s1_name]] <<- ifelse(xc < 0, xc, 0)
    data[[s2_name]] <<- ifelse(xc > 0, xc, 0)
    
    attr(data[[s1_name]], "horm_vrt") <<- vrt_value
    attr(data[[s2_name]], "horm_vrt") <<- vrt_value
    attr(data[[s1_name]], "base_name") <<- base_name
    attr(data[[s2_name]], "base_name") <<- base_name
  }
  
  # Always create base version (empirical horm_vrt)
  make_pair(horm_vrt, suffix = "")
  
  # Create SD-shifted versions if SD is provided
  if (!is.null(horm_vrt_sd)) {
    for (off in sd_offsets) {
      
      suffix <- if (off < 0) {
        paste0("m", abs(off), "SD_")
      } else {
        paste0("p", off, "SD_")
      }
      
      make_pair(horm_vrt + off * horm_vrt_sd, suffix)
    }
  }
  
  data
}

# applied 
names(df)

df_VRT <- df %>%
#------------------------
####### Negative urgency
#------------------------
  bendItLikeHorm(
    base_name   = "HSESnURG",
    adv_var     = WP_fHSES,
    horm_vrt    = 6.14,
#horm_vrt_sd = 1.30
  ) %>%
  bendItLikeHorm(
    base_name   = "CogEnnURG",
    adv_var     = WP_mLkCogEn,
    horm_vrt    = 4.73,
    #horm_vrt_sd = 1.30
  ) %>%
  bendItLikeHorm(
    base_name   = "SCHnURG_1",
    adv_var     = WP_mSchAdv_1,
    horm_vrt    = 4.44,
    #horm_vrt_sd = 1.30
  ) %>%
  bendItLikeHorm(
    base_name   = "SCHnURG_3",
    adv_var     = WP_mSchAdv_3,
    horm_vrt    = 4.44,
    #horm_vrt_sd = 1.30
  ) %>%
  bendItLikeHorm(
    base_name   = "SCHnURG_5",
    adv_var     = WP_mSchAdv_5,
    horm_vrt    = 4.44,
    #horm_vrt_sd = 1.30
  ) %>%
  bendItLikeHorm(
    base_name   = "SCHnURG_7",
    adv_var     = WP_mSchAdv_7,
    horm_vrt    = 4.44,
    #horm_vrt_sd = 1.30
  ) %>%
  bendItLikeHorm(
    base_name   = "SCHnURG_9",
    adv_var     = WP_mSchAdv_9,
    horm_vrt    = 4.44,
    #horm_vrt_sd = 1.30
  ) %>%
  bendItLikeHorm(
    base_name   = "FCONnURG_1",
    adv_var     = WP_FAMCONY_1,
    horm_vrt    = 5.56,
    #horm_vrt_sd = 1.30
  ) %>%
  bendItLikeHorm(
    base_name   = "FCONnURG_3",
    adv_var     = WP_FAMCONY_3,
    horm_vrt    = 5.56,
    #horm_vrt_sd = 1.30
  )%>%
  bendItLikeHorm(
    base_name   = "FCONnURG_5",
    adv_var     = WP_FAMCONY_5,
    horm_vrt    = 5.56,
    #horm_vrt_sd = 1.30
  )%>%
  bendItLikeHorm(
    base_name   = "FCONnURG_7",
    adv_var     = WP_FAMCONY_7,
    horm_vrt    = 5.56,
    #horm_vrt_sd = 1.30
  )%>%
  bendItLikeHorm(
    base_name   = "FCONnURG_9",
    adv_var     = WP_FAMCONY_9,
    horm_vrt    = 5.56,
    #horm_vrt_sd = 1.30
  ) %>%
  bendItLikeHorm(
    base_name   = "FCONnURG_11",
    adv_var     = WP_FAMCONY_11,
    horm_vrt    = 5.56,
    #horm_vrt_sd = 1.30
  ) %>%
  bendItLikeHorm(
    base_name   = "NBSFnURG_1",
    adv_var     = WP_mNBdang_1,
    horm_vrt    = 1.51,
    #horm_vrt_sd = 1.30
  ) %>%
  bendItLikeHorm(
    base_name   = "NBSFnURG_3",
    adv_var     = WP_mNBdang_3,
    horm_vrt    = 1.51,
    #horm_vrt_sd = 1.30
  ) %>%
  bendItLikeHorm(
    base_name   = "NBSFnURG_5",
    adv_var     = WP_mNBdang_5,
    horm_vrt    = 1.51,
    #horm_vrt_sd = 1.30
  ) %>%
  bendItLikeHorm(
    base_name   = "NBSFnURG_7",
    adv_var     = WP_mNBdang_7,
    horm_vrt    = 1.51,
    #horm_vrt_sd = 1.30
  ) %>%
  bendItLikeHorm(
    base_name   = "NBSFnURG_9",
    adv_var     = WP_mNBdang_9,
    horm_vrt    = 1.51,
    #horm_vrt_sd = 1.30
  ) %>%
  bendItLikeHorm(
    base_name   = "NBSFnURG_11",
    adv_var     = WP_mNBdang_11,
    horm_vrt    = 1.51,
    #horm_vrt_sd = 1.30
  ) %>%
#------------------------
####### Positive urgency
#------------------------
bendItLikeHorm(
    base_name   = "HSESpURG",
    adv_var     = WP_fHSES,
    horm_vrt    = 6.28,
#horm_vrt_sd = 1.30
  ) %>%
  bendItLikeHorm(
    base_name   = "CogEnpURG",
    adv_var     = WP_mLkCogEn,
    horm_vrt    = 4.80,
    #horm_vrt_sd = 1.30
  ) %>%
  bendItLikeHorm(
    base_name   = "SCHpURG_1",
    adv_var     = WP_mSchAdv_1,
    horm_vrt    = 4.44,
    #horm_vrt_sd = 1.30
  ) %>%
  bendItLikeHorm(
    base_name   = "SCHpURG_3",
    adv_var     = WP_mSchAdv_3,
    horm_vrt    = 4.44,
    #horm_vrt_sd = 1.30
  ) %>%
  bendItLikeHorm(
    base_name   = "SCHpURG_5",
    adv_var     = WP_mSchAdv_5,
    horm_vrt    = 4.44,
    #horm_vrt_sd = 1.30
  ) %>%
  bendItLikeHorm(
    base_name   = "SCHpURG_7",
    adv_var     = WP_mSchAdv_7,
    horm_vrt    = 4.44,
    #horm_vrt_sd = 1.30
  ) %>%
  bendItLikeHorm(
    base_name   = "SCHpURG_9",
    adv_var     = WP_mSchAdv_9,
    horm_vrt    = 4.44,
    #horm_vrt_sd = 1.30
  ) %>%
  bendItLikeHorm(
    base_name   = "FCONpURG_1",
    adv_var     = WP_FAMCONY_1,
    horm_vrt    = 5.56,
    #horm_vrt_sd = 1.30
  ) %>%
  bendItLikeHorm(
    base_name   = "FCONpURG_3",
    adv_var     = WP_FAMCONY_3,
    horm_vrt    = 5.56,
    #horm_vrt_sd = 1.30
  )%>%
  bendItLikeHorm(
    base_name   = "FCONpURG_5",
    adv_var     = WP_FAMCONY_5,
    horm_vrt    = 5.56,
    #horm_vrt_sd = 1.30
  )%>%
  bendItLikeHorm(
    base_name   = "FCONpURG_7",
    adv_var     = WP_FAMCONY_7,
    horm_vrt    = 5.56,
    #horm_vrt_sd = 1.30
  )%>%
  bendItLikeHorm(
    base_name   = "FCONpURG_9",
    adv_var     = WP_FAMCONY_9,
    horm_vrt    = 5.56,
    #horm_vrt_sd = 1.30
  ) %>%
  bendItLikeHorm(
    base_name   = "FCONpURG_11",
    adv_var     = WP_FAMCONY_11,
    horm_vrt    = 5.56,
    #horm_vrt_sd = 1.30
  ) %>%
  bendItLikeHorm(
    base_name   = "NBSFpURG_1",
    adv_var     = WP_mNBdang_1,
    horm_vrt    = .833,
    #horm_vrt_sd = 1.30
  ) %>%
  bendItLikeHorm(
    base_name   = "NBSFpURG_3",
    adv_var     = WP_mNBdang_3,
    horm_vrt    = .833,
    #horm_vrt_sd = 1.30
  ) %>%
  bendItLikeHorm(
    base_name   = "NBSFpURG_5",
    adv_var     = WP_mNBdang_5,
    horm_vrt    = .833,
    #horm_vrt_sd = 1.30
  ) %>%
  bendItLikeHorm(
    base_name   = "NBSFpURG_7",
    adv_var     = WP_mNBdang_7,
    horm_vrt    = .833,
    #horm_vrt_sd = 1.30
  ) %>%
  bendItLikeHorm(
    base_name   = "NBSFpURG_9",
    adv_var     = WP_mNBdang_9,
    horm_vrt    = .833,
    #horm_vrt_sd = 1.30
  ) %>%
  bendItLikeHorm(
    base_name   = "NBSFpURG_11",
    adv_var     = WP_mNBdang_11,
    horm_vrt    = .833,
    #horm_vrt_sd = 1.30
  ) %>%
  #------------------------
####### Attune
#------------------------
bendItLikeHorm(
  base_name   = "HSESattn",
  adv_var     = WP_fHSES,
  horm_vrt    = 3.15,
  #horm_vrt_sd = 1.30
) %>%
  bendItLikeHorm(
    base_name   = "CogEnattn",
    adv_var     = WP_mLkCogEn,
    horm_vrt    = 5.35,
    #horm_vrt_sd = 1.30
  ) %>%
  bendItLikeHorm(
    base_name   = "SCHattn_1",
    adv_var     = WP_mSchAdv_1,
    horm_vrt    = 4.44,
    #horm_vrt_sd = 1.30
  ) %>%
  bendItLikeHorm(
    base_name   = "SCHattn_3",
    adv_var     = WP_mSchAdv_3,
    horm_vrt    = 4.44,
    #horm_vrt_sd = 1.30
  ) %>%
  bendItLikeHorm(
    base_name   = "SCHattn_5",
    adv_var     = WP_mSchAdv_5,
    horm_vrt    = 4.44,
    #horm_vrt_sd = 1.30
  ) %>%
  bendItLikeHorm(
    base_name   = "SCHattn_7",
    adv_var     = WP_mSchAdv_7,
    horm_vrt    = 4.44,
    #horm_vrt_sd = 1.30
  ) %>%
  bendItLikeHorm(
    base_name   = "SCHattn_9",
    adv_var     = WP_mSchAdv_9,
    horm_vrt    = 4.44,
    #horm_vrt_sd = 1.30
  ) %>%
  bendItLikeHorm(
    base_name   = "FCONattn_1",
    adv_var     = WP_FAMCONY_1,
    horm_vrt    = 4.73,
    #horm_vrt_sd = 1.30
  ) %>%
  bendItLikeHorm(
    base_name   = "FCONattn_3",
    adv_var     = WP_FAMCONY_3,
    horm_vrt    = 4.73,
    #horm_vrt_sd = 1.30
  )%>%
  bendItLikeHorm(
    base_name   = "FCONattn_5",
    adv_var     = WP_FAMCONY_5,
    horm_vrt    = 4.73,
    #horm_vrt_sd = 1.30
  )%>%
  bendItLikeHorm(
    base_name   = "FCONattn_7",
    adv_var     = WP_FAMCONY_7,
    horm_vrt    = 4.73,
    #horm_vrt_sd = 1.30
  )%>%
  bendItLikeHorm(
    base_name   = "FCONattn_9",
    adv_var     = WP_FAMCONY_9,
    horm_vrt    = 4.73,
    #horm_vrt_sd = 1.30
  ) %>%
  bendItLikeHorm(
    base_name   = "FCONattn_11",
    adv_var     = WP_FAMCONY_11,
    horm_vrt    = 4.73,
    #horm_vrt_sd = 1.30
  ) %>%
  bendItLikeHorm(
    base_name   = "NBSFattn_1",
    adv_var     = WP_mNBdang_1,
    horm_vrt    = 3.33,
    #horm_vrt_sd = 1.30
  ) %>%
  bendItLikeHorm(
    base_name   = "NBSFattn_3",
    adv_var     = WP_mNBdang_3,
    horm_vrt    = 3.33,
    #horm_vrt_sd = 1.30
  ) %>%
  bendItLikeHorm(
    base_name   = "NBSFattn_5",
    adv_var     = WP_mNBdang_5,
    horm_vrt    = 3.33,
    #horm_vrt_sd = 1.30
  ) %>%
  bendItLikeHorm(
    base_name   = "NBSFattn_7",
    adv_var     = WP_mNBdang_7,
    horm_vrt    = 3.33,
    #horm_vrt_sd = 1.30
  ) %>%
  bendItLikeHorm(
    base_name   = "NBSFattn_9",
    adv_var     = WP_mNBdang_9,
    horm_vrt    = 3.33,
    #horm_vrt_sd = 1.30
  ) %>%
  bendItLikeHorm(
    base_name   = "NBSFattn_11",
    adv_var     = WP_mNBdang_11,
    horm_vrt    = 3.33,
    #horm_vrt_sd = 1.30
  ) %>%
#------------------------
####### Distract
#------------------------
bendItLikeHorm(
  base_name   = "HSESdist",
  adv_var     = WP_fHSES,
  horm_vrt    = 2.30,
  #horm_vrt_sd = 1.30
) %>%
  bendItLikeHorm(
    base_name   = "CogEndist",
    adv_var     = WP_mLkCogEn,
    horm_vrt    = 2.54,
    #horm_vrt_sd = 1.30
  ) %>%
  bendItLikeHorm(
    base_name   = "SCHdist_1",
    adv_var     = WP_mSchAdv_1,
    horm_vrt    = 1.06,
    #horm_vrt_sd = 1.30
  ) %>%
  bendItLikeHorm(
    base_name   = "SCHdist_3",
    adv_var     = WP_mSchAdv_3,
    horm_vrt    = 1.06,
    #horm_vrt_sd = 1.30
  ) %>%
  bendItLikeHorm(
    base_name   = "SCHdist_5",
    adv_var     = WP_mSchAdv_5,
    horm_vrt    = 1.06,
    #horm_vrt_sd = 1.30
  ) %>%
  bendItLikeHorm(
    base_name   = "SCHdist_7",
    adv_var     = WP_mSchAdv_7,
    horm_vrt    = 1.06,
    #horm_vrt_sd = 1.30
  ) %>%
  bendItLikeHorm(
    base_name   = "SCHdist_9",
    adv_var     = WP_mSchAdv_9,
    horm_vrt    = 1.06
    #horm_vrt_sd = 1.30
  ) %>%
  bendItLikeHorm(
    base_name   = "FCONdist_1",
    adv_var     = WP_FAMCONY_1,
    horm_vrt    = 1.11,
    #horm_vrt_sd = 1.30
  ) %>%
  bendItLikeHorm(
    base_name   = "FCONdist_3",
    adv_var     = WP_FAMCONY_3,
    horm_vrt    = 1.11,
    #horm_vrt_sd = 1.30
  )%>%
  bendItLikeHorm(
    base_name   = "FCONdist_5",
    adv_var     = WP_FAMCONY_5,
    horm_vrt    = 1.11,
    #horm_vrt_sd = 1.30
  )%>%
  bendItLikeHorm(
    base_name   = "FCONdist_7",
    adv_var     = WP_FAMCONY_7,
    horm_vrt    = 1.11,
    #horm_vrt_sd = 1.30
  )%>%
  bendItLikeHorm(
    base_name   = "FCONdist_9",
    adv_var     = WP_FAMCONY_9,
    horm_vrt    = 1.11,
    #horm_vrt_sd = 1.30
  ) %>%
  bendItLikeHorm(
    base_name   = "FCONdist_11",
    adv_var     = WP_FAMCONY_11,
    horm_vrt    = 1.11,
    #horm_vrt_sd = 1.30
  ) %>%
  bendItLikeHorm(
    base_name   = "NBSFdist_1",
    adv_var     = WP_mNBdang_1,
    horm_vrt    = 1.09,
    #horm_vrt_sd = 1.30
  ) %>%
  bendItLikeHorm(
    base_name   = "NBSFdist_3",
    adv_var     = WP_mNBdang_3,
    horm_vrt    = 1.09,
    #horm_vrt_sd = 1.30
  ) %>%
  bendItLikeHorm(
    base_name   = "NBSFdist_5",
    adv_var     = WP_mNBdang_5,
    horm_vrt    = 1.09,
    #horm_vrt_sd = 1.30
  ) %>%
  bendItLikeHorm(
    base_name   = "NBSFdist_7",
    adv_var     = WP_mNBdang_7,
    horm_vrt    = 1.09,
    #horm_vrt_sd = 1.30
  ) %>%
  bendItLikeHorm(
    base_name   = "NBSFdist_9",
    adv_var     = WP_mNBdang_9,
    horm_vrt    = 1.09,
    #horm_vrt_sd = 1.30
  ) %>%
  bendItLikeHorm(
    base_name   = "NBSFdist_11",
    adv_var     = WP_mNBdang_11,
    horm_vrt    = 1.09,
    #horm_vrt_sd = 1.30
  ) %>%
#------------------------
####### Emotional Suppression
#------------------------
bendItLikeHorm(
  base_name   = "HSESemosup",
  adv_var     = WP_fHSES,
  horm_vrt    = 3.60,
  #horm_vrt_sd = 1.30
) %>%
  bendItLikeHorm(
    base_name   = "CogEnemosup",
    adv_var     = WP_mLkCogEn,
    horm_vrt    = 4.70,
    #horm_vrt_sd = 1.30
  ) %>%
  bendItLikeHorm(
    base_name   = "SCHemosup_1",
    adv_var     = WP_mSchAdv_1,
    horm_vrt    = .56,
    #horm_vrt_sd = 1.30
  ) %>%
  bendItLikeHorm(
    base_name   = "SCHemosup_3",
    adv_var     = WP_mSchAdv_3,
    horm_vrt    = .56,
    #horm_vrt_sd = 1.30
  ) %>%
  bendItLikeHorm(
    base_name   = "SCHemosup_5",
    adv_var     = WP_mSchAdv_5,
    horm_vrt    = .56,
    #horm_vrt_sd = 1.30
  ) %>%
  bendItLikeHorm(
    base_name   = "SCHemosup_7",
    adv_var     = WP_mSchAdv_7,
    horm_vrt    = .56,
    #horm_vrt_sd = 1.30
  ) %>%
  bendItLikeHorm(
    base_name   = "SCHemosup_9",
    adv_var     = WP_mSchAdv_9,
    horm_vrt    = .56,
    #horm_vrt_sd = 1.30
  ) %>%
  bendItLikeHorm(
    base_name   = "FCONemosup_1",
    adv_var     = WP_FAMCONY_1,
    horm_vrt    = 5.56,
    #horm_vrt_sd = 1.30
  ) %>%
  bendItLikeHorm(
    base_name   = "FCONemosup_3",
    adv_var     = WP_FAMCONY_3,
    horm_vrt    = 5.56,
    #horm_vrt_sd = 1.30
  )%>%
  bendItLikeHorm(
    base_name   = "FCONemosup_5",
    adv_var     = WP_FAMCONY_5,
    horm_vrt    = 5.56,
    #horm_vrt_sd = 1.30
  )%>%
  bendItLikeHorm(
    base_name   = "FCONemosup_7",
    adv_var     = WP_FAMCONY_7,
    horm_vrt    = 5.56,
    #horm_vrt_sd = 1.30
  )%>%
  bendItLikeHorm(
    base_name   = "FCONemosup_9",
    adv_var     = WP_FAMCONY_9,
    horm_vrt    = 5.56,
    #horm_vrt_sd = 1.30
  ) %>%
  bendItLikeHorm(
    base_name   = "FCONemosup_11",
    adv_var     = WP_FAMCONY_11,
    horm_vrt    = 5.56,
    #horm_vrt_sd = 1.30
  ) %>%
  bendItLikeHorm(
    base_name   = "NBSFemosup_1",
    adv_var     = WP_mNBdang_1,
    horm_vrt    = .83,
    #horm_vrt_sd = 1.30
  ) %>%
  bendItLikeHorm(
    base_name   = "NBSFemosup_3",
    adv_var     = WP_mNBdang_3,
    horm_vrt    = .83,
    #horm_vrt_sd = 1.30
  ) %>%
  bendItLikeHorm(
    base_name   = "NBSFemosup_5",
    adv_var     = WP_mNBdang_5,
    horm_vrt    = .83,
    #horm_vrt_sd = 1.30
  ) %>%
  bendItLikeHorm(
    base_name   = "NBSFemosup_7",
    adv_var     = WP_mNBdang_7,
    horm_vrt    = .83,
    #horm_vrt_sd = 1.30
  ) %>%
  bendItLikeHorm(
    base_name   = "NBSFemosup_9",
    adv_var     = WP_mNBdang_9,
    horm_vrt    = .83,
    #horm_vrt_sd = 1.30
  ) %>%
  bendItLikeHorm(
    base_name   = "NBSFemosup_11",
    adv_var     = WP_mNBdang_11,
    horm_vrt    = .83,
    #horm_vrt_sd = 1.30
  ) %>%
#------------------------
####### Reappraisal
#------------------------
bendItLikeHorm(
  base_name   = "HSESreapp",
  adv_var     = WP_fHSES,
  horm_vrt    = 7.42,
  #horm_vrt_sd = 1.30
) %>%
  bendItLikeHorm(
    base_name   = "CogEnreapp",
    adv_var     = WP_mLkCogEn,
    horm_vrt    = 4.41,
    #horm_vrt_sd = 1.30
  ) %>%
  bendItLikeHorm(
    base_name   = "SCHreapp_1",
    adv_var     = WP_mSchAdv_1,
    horm_vrt    = .56,
    #horm_vrt_sd = 1.30
  ) %>%
  bendItLikeHorm(
    base_name   = "SCHreapp_3",
    adv_var     = WP_mSchAdv_3,
    horm_vrt    = .56,
    #horm_vrt_sd = 1.30
  ) %>%
  bendItLikeHorm(
    base_name   = "SCHreapp_5",
    adv_var     = WP_mSchAdv_5,
    horm_vrt    = .56,
    #horm_vrt_sd = 1.30
  ) %>%
  bendItLikeHorm(
    base_name   = "SCHreapp_7",
    adv_var     = WP_mSchAdv_7,
    horm_vrt    = .56,
    #horm_vrt_sd = 1.30
  ) %>%
  bendItLikeHorm(
    base_name   = "SCHreapp_9",
    adv_var     = WP_mSchAdv_9,
    horm_vrt    = .56,
    #horm_vrt_sd = 1.30
  ) %>%
  bendItLikeHorm(
    base_name   = "FCONreapp_1",
    adv_var     = WP_FAMCONY_1,
    horm_vrt    = 5.56,
    #horm_vrt_sd = 1.30
  ) %>%
  bendItLikeHorm(
    base_name   = "FCONreapp_3",
    adv_var     = WP_FAMCONY_3,
    horm_vrt    = 5.56,
    #horm_vrt_sd = 1.30
  )%>%
  bendItLikeHorm(
    base_name   = "FCONreapp_5",
    adv_var     = WP_FAMCONY_5,
    horm_vrt    = 5.56,
    #horm_vrt_sd = 1.30
  )%>%
  bendItLikeHorm(
    base_name   = "FCONreapp_7",
    adv_var     = WP_FAMCONY_7,
    horm_vrt    = 5.56,
    #horm_vrt_sd = 1.30
  )%>%
  bendItLikeHorm(
    base_name   = "FCONreapp_9",
    adv_var     = WP_FAMCONY_9,
    horm_vrt    = 5.56,
    #horm_vrt_sd = 1.30
  ) %>%
  bendItLikeHorm(
    base_name   = "FCONreapp_11",
    adv_var     = WP_FAMCONY_11,
    horm_vrt    = 5.56,
    #horm_vrt_sd = 1.30
  ) %>%
  bendItLikeHorm(
    base_name   = "NBSFreapp_1",
    adv_var     = WP_mNBdang_1,
    horm_vrt    = 3.29,
    #horm_vrt_sd = 1.30
  ) %>%
  bendItLikeHorm(
    base_name   = "NBSFreapp_3",
    adv_var     = WP_mNBdang_3,
    horm_vrt    = 3.29,
    #horm_vrt_sd = 1.30
  ) %>%
  bendItLikeHorm(
    base_name   = "NBSFreapp_5",
    adv_var     = WP_mNBdang_5,
    horm_vrt    = 3.29,
    #horm_vrt_sd = 1.30
  ) %>%
  bendItLikeHorm(
    base_name   = "NBSFreapp_7",
    adv_var     = WP_mNBdang_7,
    horm_vrt    = 3.29,
    #horm_vrt_sd = 1.30
  ) %>%
  bendItLikeHorm(
    base_name   = "NBSFreapp_9",
    adv_var     = WP_mNBdang_9,
    horm_vrt    = 3.29,
    #horm_vrt_sd = 1.30
  ) %>%
  bendItLikeHorm(
    base_name   = "NBSFreapp_11",
    adv_var     = WP_mNBdang_11,
    horm_vrt    = 3.29,
    #horm_vrt_sd = 1.30
  ) 

###############################################################################
############### add back in parental monitoring - fixed factor construction 

# remove the wromng (5-item) parental monitoring 

names(df_VRT)

df_VRT <- df_VRT %>%
  select(-c("LMC_WP_mPmon_1", "LMC_WP_mPmon_3", "LMC_WP_mPmon_5", 
            "LMC_WP_mPmon_7",  "LMC_WP_mPmon_9", "LMC_WP_mPmon_11"))

# loadi n the correct (4-item) parental monitoring 
Pmon_df <- read.csv("PARENTAL_MONITORING_4_ITEM_SUB_4.27.26.csv")

# merge it in 

df__VRT_PMON <- left_join(df_VRT, Pmon_df, by = "subID")


###############################################################################
############### Get S1, S2, and Mod means + wave specific deviations

## use CenterForge()
CenterForge <- function(data, vars,
                        types = c("GMC", "LMC", "PMC"),
                        save_means = TRUE,
                        set_name = NULL,
                        check_numeric = TRUE,
                        poly = c("none", "quadratic", "cubic")) {
  
  if (!is.data.frame(data)) {
    stop("data must be a data frame.")
  }
  
  if (missing(vars) || length(vars) == 0) {
    stop("You must provide at least one variable name in 'vars'.")
  }
  
  if (!all(vars %in% names(data))) {
    missing_vars <- vars[!vars %in% names(data)]
    stop("These variables are not in the data: ",
         paste(missing_vars, collapse = ", "))
  }
  
  types <- unique(toupper(types))
  allowed_types <- c("GMC", "LMC", "PMC")
  
  if (!all(types %in% allowed_types)) {
    bad_types <- types[!types %in% allowed_types]
    stop("Invalid type(s): ", paste(bad_types, collapse = ", "),
         ". Allowed types are: ", paste(allowed_types, collapse = ", "))
  }
  
  poly <- match.arg(poly)
  
  df_sub <- data[, vars, drop = FALSE]
  
  if (check_numeric) {
    non_numeric <- vars[!vapply(df_sub, is.numeric, logical(1))]
    if (length(non_numeric) > 0) {
      stop("These variables are not numeric: ",
           paste(non_numeric, collapse = ", "))
    }
  }
  
  out <- data
  
  if (is.null(set_name)) {
    set_name <- "SET"
  }
  
  # helper to add polynomial terms
  add_poly_terms <- function(df, centered_names, poly_type) {
    if (poly_type %in% c("quadratic", "cubic")) {
      for (nm in centered_names) {
        df[[paste0("Q", nm)]] <- df[[nm]]^2
      }
    }
    
    if (poly_type == "cubic") {
      for (nm in centered_names) {
        df[[paste0("C", nm)]] <- df[[nm]]^3
      }
    }
    
    df
  }
  
  # GMC
  if ("GMC" %in% types) {
    var_means <- vapply(df_sub, mean, numeric(1), na.rm = TRUE)
    
    if (save_means) {
      for (j in seq_along(var_means)) {
        out[[paste0("MEAN_GMC_", vars[j])]] <- var_means[j]
      }
    }
    
    gmc_df <- as.data.frame(
      Map(function(x, m) x - m, df_sub, var_means)
    )
    gmc_names <- paste0("GMC_", vars)
    names(gmc_df) <- gmc_names
    out <- cbind(out, gmc_df)
    
    out <- add_poly_terms(out, gmc_names, poly)
  }
  
  # LMC
  if ("LMC" %in% types) {
    pooled_mean <- mean(as.matrix(df_sub), na.rm = TRUE)
    
    if (save_means) {
      out[[paste0("MEAN_LMC_", set_name)]] <- pooled_mean
    }
    
    lmc_df <- as.data.frame(
      lapply(df_sub, function(x) x - pooled_mean)
    )
    lmc_names <- paste0("LMC_", vars)
    names(lmc_df) <- lmc_names
    out <- cbind(out, lmc_df)
    
    out <- add_poly_terms(out, lmc_names, poly)
  }
  
  # PMC
  if ("PMC" %in% types) {
    person_means <- rowMeans(df_sub, na.rm = TRUE)
    person_means[is.nan(person_means)] <- NA
    
    if (save_means) {
      out[[paste0("MEAN_PMC_", set_name)]] <- person_means
    }
    
    pmc_df <- as.data.frame(
      lapply(df_sub, function(x) x - person_means)
    )
    pmc_names <- paste0("PMC_", vars)
    names(pmc_df) <- pmc_names
    out <- cbind(out, pmc_df)
    
    out <- add_poly_terms(out, pmc_names, poly)
  }
  
  return(out)
}


## applied 

# -------------------------
# LMC moderators / covariates
# -------------------------
names(df__VRT_PMON)
df__VRT_PMON <- CenterForge(
  data = df__VRT_PMON,
  vars = c("LMC_DysBpM_5", "LMC_DysBpM_7", "LMC_DysBpM_9", "LMC_DysBpM_11"),
  types = "PMC",
  save_means = TRUE,
  set_name = "mLMC_DysBpM",
  poly = "none"
)

df__VRT_PMON <- CenterForge(
  data = df__VRT_PMON,
  vars = c("LMC_SysBpM_5", "LMC_SysBpM_7", "LMC_SysBpM_9", "LMC_SysBpM_11"),
  types = "PMC",
  save_means = TRUE,
  set_name = "mLMC_SysBpM",
  poly = "none"
)

df__VRT_PMON <- CenterForge(
  data = df__VRT_PMON,
  vars = c("LMC_WP_mPmon_1", "LMC_WP_mPmon_3", "LMC_WP_mPmon_5", "LMC_WP_mPmon_7", "LMC_WP_mPmon_9", "LMC_WP_mPmon_11"),
  types = "PMC",
  save_means = TRUE,
  set_name = "mLMC_mPmon",
  poly = "none"
)

df__VRT_PMON <- CenterForge(
  data = df__VRT_PMON,
  vars = c("LMC_WP_mPNH_5", "LMC_WP_mPNH_7", "LMC_WP_mPNH_9"),
  types = "PMC",
  save_means = TRUE,
  set_name = "mLMC_mPNH",
  poly = "none"
)


# -------------------------
# nURG
# -------------------------

df__VRT_PMON <- CenterForge(
  data = df__VRT_PMON,
  vars = c("S1_SCHnURG_1", "S1_SCHnURG_3", "S1_SCHnURG_5", "S1_SCHnURG_7", "S1_SCHnURG_9"),
  types = "PMC",
  save_means = TRUE,
  set_name = "MS1_SCHnURG",
  poly = "none"
)

df__VRT_PMON <- CenterForge(
  data = df__VRT_PMON,
  vars = c("S2_SCHnURG_1", "S2_SCHnURG_3", "S2_SCHnURG_5", "S2_SCHnURG_7", "S2_SCHnURG_9"),
  types = "PMC",
  save_means = TRUE,
  set_name = "MS2_SCHnURG",
  poly = "none"
)

df__VRT_PMON <- CenterForge(
  data = df__VRT_PMON,
  vars = c("S1_FCONnURG_1", "S1_FCONnURG_3", "S1_FCONnURG_5", "S1_FCONnURG_7", "S1_FCONnURG_9", "S1_FCONnURG_11"),
  types = "PMC",
  save_means = TRUE,
  set_name = "MS1_FCONnURG",
  poly = "none"
)

df__VRT_PMON <- CenterForge(
  data = df__VRT_PMON,
  vars = c("S2_FCONnURG_1", "S2_FCONnURG_3", "S2_FCONnURG_5", "S2_FCONnURG_7", "S2_FCONnURG_9", "S2_FCONnURG_11"),
  types = "PMC",
  save_means = TRUE,
  set_name = "MS2_FCONnURG",
  poly = "none"
)

df__VRT_PMON <- CenterForge(
  data = df__VRT_PMON,
  vars = c("S1_NBSFnURG_1", "S1_NBSFnURG_3", "S1_NBSFnURG_5", "S1_NBSFnURG_7", "S1_NBSFnURG_9", "S1_NBSFnURG_11"),
  types = "PMC",
  save_means = TRUE,
  set_name = "MS1_NBSFnURG",
  poly = "none"
)

df__VRT_PMON <- CenterForge(
  data = df__VRT_PMON,
  vars = c("S2_NBSFnURG_1", "S2_NBSFnURG_3", "S2_NBSFnURG_5", "S2_NBSFnURG_7", "S2_NBSFnURG_9", "S2_NBSFnURG_11"),
  types = "PMC",
  save_means = TRUE,
  set_name = "MS2_NBSFnURG",
  poly = "none"
)

# -------------------------
# pURG
# -------------------------

df__VRT_PMON <- CenterForge(
  data = df__VRT_PMON,
  vars = c("S1_SCHpURG_1", "S1_SCHpURG_3", "S1_SCHpURG_5", "S1_SCHpURG_7", "S1_SCHpURG_9"),
  types = "PMC",
  save_means = TRUE,
  set_name = "MS1_SCHpURG",
  poly = "none"
)

df__VRT_PMON <- CenterForge(
  data = df__VRT_PMON,
  vars = c("S2_SCHpURG_1", "S2_SCHpURG_3", "S2_SCHpURG_5", "S2_SCHpURG_7", "S2_SCHpURG_9"),
  types = "PMC",
  save_means = TRUE,
  set_name = "MS2_SCHpURG",
  poly = "none"
)

df__VRT_PMON <- CenterForge(
  data = df__VRT_PMON,
  vars = c("S1_FCONpURG_1", "S1_FCONpURG_3", "S1_FCONpURG_5", "S1_FCONpURG_7", "S1_FCONpURG_9", "S1_FCONpURG_11"),
  types = "PMC",
  save_means = TRUE,
  set_name = "MS1_FCONpURG",
  poly = "none"
)

df__VRT_PMON <- CenterForge(
  data = df__VRT_PMON,
  vars = c("S2_FCONpURG_1", "S2_FCONpURG_3", "S2_FCONpURG_5", "S2_FCONpURG_7", "S2_FCONpURG_9", "S2_FCONpURG_11"),
  types = "PMC",
  save_means = TRUE,
  set_name = "MS2_FCONpURG",
  poly = "none"
)

df__VRT_PMON <- CenterForge(
  data = df__VRT_PMON,
  vars = c("S1_NBSFpURG_1", "S1_NBSFpURG_3", "S1_NBSFpURG_5", "S1_NBSFpURG_7", "S1_NBSFpURG_9", "S1_NBSFpURG_11"),
  types = "PMC",
  save_means = TRUE,
  set_name = "MS1_NBSFpURG",
  poly = "none"
)

df__VRT_PMON <- CenterForge(
  data = df__VRT_PMON,
  vars = c("S2_NBSFpURG_1", "S2_NBSFpURG_3", "S2_NBSFpURG_5", "S2_NBSFpURG_7", "S2_NBSFpURG_9", "S2_NBSFpURG_11"),
  types = "PMC",
  save_means = TRUE,
  set_name = "MS2_NBSFpURG",
  poly = "none"
)

# -------------------------
# attn
# -------------------------

df__VRT_PMON <- CenterForge(
  data = df__VRT_PMON,
  vars = c("S1_SCHattn_1", "S1_SCHattn_3", "S1_SCHattn_5", "S1_SCHattn_7", "S1_SCHattn_9"),
  types = "PMC",
  save_means = TRUE,
  set_name = "MS1_SCHattn",
  poly = "none"
)

df__VRT_PMON <- CenterForge(
  data = df__VRT_PMON,
  vars = c("S2_SCHattn_1", "S2_SCHattn_3", "S2_SCHattn_5", "S2_SCHattn_7", "S2_SCHattn_9"),
  types = "PMC",
  save_means = TRUE,
  set_name = "MS2_SCHattn",
  poly = "none"
)

df__VRT_PMON <- CenterForge(
  data = df__VRT_PMON,
  vars = c("S1_FCONattn_1", "S1_FCONattn_3", "S1_FCONattn_5", "S1_FCONattn_7", "S1_FCONattn_9", "S1_FCONattn_11"),
  types = "PMC",
  save_means = TRUE,
  set_name = "MS1_FCONattn",
  poly = "none"
)

df__VRT_PMON <- CenterForge(
  data = df__VRT_PMON,
  vars = c("S2_FCONattn_1", "S2_FCONattn_3", "S2_FCONattn_5", "S2_FCONattn_7", "S2_FCONattn_9", "S2_FCONattn_11"),
  types = "PMC",
  save_means = TRUE,
  set_name = "MS2_FCONattn",
  poly = "none"
)

df__VRT_PMON <- CenterForge(
  data = df__VRT_PMON,
  vars = c("S1_NBSFattn_1", "S1_NBSFattn_3", "S1_NBSFattn_5", "S1_NBSFattn_7", "S1_NBSFattn_9", "S1_NBSFattn_11"),
  types = "PMC",
  save_means = TRUE,
  set_name = "MS1_NBSFattn",
  poly = "none"
)

df__VRT_PMON <- CenterForge(
  data = df__VRT_PMON,
  vars = c("S2_NBSFattn_1", "S2_NBSFattn_3", "S2_NBSFattn_5", "S2_NBSFattn_7", "S2_NBSFattn_9", "S2_NBSFattn_11"),
  types = "PMC",
  save_means = TRUE,
  set_name = "MS2_NBSFattn",
  poly = "none"
)

# -------------------------
# dist (goals)
# -------------------------

df__VRT_PMON <- CenterForge(
  data = df__VRT_PMON,
  vars = c("S1_SCHdist_1", "S1_SCHdist_3", "S1_SCHdist_5", "S1_SCHdist_7", "S1_SCHdist_9"),
  types = "PMC",
  save_means = TRUE,
  set_name = "MS1_SCHdist",
  poly = "none"
)

df__VRT_PMON <- CenterForge(
  data = df__VRT_PMON,
  vars = c("S2_SCHdist_1", "S2_SCHdist_3", "S2_SCHdist_5", "S2_SCHdist_7", "S2_SCHdist_9"),
  types = "PMC",
  save_means = TRUE,
  set_name = "MS2_SCHdist",
  poly = "none"
)

df__VRT_PMON <- CenterForge(
  data = df__VRT_PMON,
  vars = c("S1_FCONdist_1", "S1_FCONdist_3", "S1_FCONdist_5", "S1_FCONdist_7", "S1_FCONdist_9", "S1_FCONdist_11"),
  types = "PMC",
  save_means = TRUE,
  set_name = "MS1_FCONdist",
  poly = "none"
)

df__VRT_PMON <- CenterForge(
  data = df__VRT_PMON,
  vars = c("S2_FCONdist_1", "S2_FCONdist_3", "S2_FCONdist_5", "S2_FCONdist_7", "S2_FCONdist_9", "S2_FCONdist_11"),
  types = "PMC",
  save_means = TRUE,
  set_name = "MS2_FCONdist",
  poly = "none"
)

df__VRT_PMON <- CenterForge(
  data = df__VRT_PMON,
  vars = c("S1_NBSFdist_1", "S1_NBSFdist_3", "S1_NBSFdist_5", "S1_NBSFdist_7", "S1_NBSFdist_9", "S1_NBSFdist_11"),
  types = "PMC",
  save_means = TRUE,
  set_name = "MS1_NBSFdist",
  poly = "none"
)

df__VRT_PMON <- CenterForge(
  data = df__VRT_PMON,
  vars = c("S2_NBSFdist_1", "S2_NBSFdist_3", "S2_NBSFdist_5", "S2_NBSFdist_7", "S2_NBSFdist_9", "S2_NBSFdist_11"),
  types = "PMC",
  save_means = TRUE,
  set_name = "MS2_NBSFdist",
  poly = "none"
)

# -------------------------
# emosup
# -------------------------

df__VRT_PMON <- CenterForge(
  data = df__VRT_PMON,
  vars = c("S1_SCHemosup_1", "S1_SCHemosup_3", "S1_SCHemosup_5", "S1_SCHemosup_7", "S1_SCHemosup_9"),
  types = "PMC",
  save_means = TRUE,
  set_name = "MS1_SCHemosup",
  poly = "none"
)

df__VRT_PMON <- CenterForge(
  data = df__VRT_PMON,
  vars = c("S2_SCHemosup_1", "S2_SCHemosup_3", "S2_SCHemosup_5", "S2_SCHemosup_7", "S2_SCHemosup_9"),
  types = "PMC",
  save_means = TRUE,
  set_name = "MS2_SCHemosup",
  poly = "none"
)

df__VRT_PMON <- CenterForge(
  data = df__VRT_PMON,
  vars = c("S1_FCONemosup_1", "S1_FCONemosup_3", "S1_FCONemosup_5", "S1_FCONemosup_7", "S1_FCONemosup_9", "S1_FCONemosup_11"),
  types = "PMC",
  save_means = TRUE,
  set_name = "MS1_FCONemosup",
  poly = "none"
)

df__VRT_PMON <- CenterForge(
  data = df__VRT_PMON,
  vars = c("S2_FCONemosup_1", "S2_FCONemosup_3", "S2_FCONemosup_5", "S2_FCONemosup_7", "S2_FCONemosup_9", "S2_FCONemosup_11"),
  types = "PMC",
  save_means = TRUE,
  set_name = "MS2_FCONemosup",
  poly = "none"
)

df__VRT_PMON <- CenterForge(
  data = df__VRT_PMON,
  vars = c("S1_NBSFemosup_1", "S1_NBSFemosup_3", "S1_NBSFemosup_5", "S1_NBSFemosup_7", "S1_NBSFemosup_9", "S1_NBSFemosup_11"),
  types = "PMC",
  save_means = TRUE,
  set_name = "MS1_NBSFemosup",
  poly = "none"
)

df__VRT_PMON <- CenterForge(
  data = df__VRT_PMON,
  vars = c("S2_NBSFemosup_1", "S2_NBSFemosup_3", "S2_NBSFemosup_5", "S2_NBSFemosup_7", "S2_NBSFemosup_9", "S2_NBSFemosup_11"),
  types = "PMC",
  save_means = TRUE,
  set_name = "MS2_NBSFemosup",
  poly = "none"
)

# -------------------------
# reapp
# -------------------------

df__VRT_PMON <- CenterForge(
  data = df__VRT_PMON,
  vars = c("S1_SCHreapp_1", "S1_SCHreapp_3", "S1_SCHreapp_5", "S1_SCHreapp_7", "S1_SCHreapp_9"),
  types = "PMC",
  save_means = TRUE,
  set_name = "MS1_SCHreapp",
  poly = "none"
)

df__VRT_PMON <- CenterForge(
  data = df__VRT_PMON,
  vars = c("S2_SCHreapp_1", "S2_SCHreapp_3", "S2_SCHreapp_5", "S2_SCHreapp_7", "S2_SCHreapp_9"),
  types = "PMC",
  save_means = TRUE,
  set_name = "MS2_SCHreapp",
  poly = "none"
)

df__VRT_PMON <- CenterForge(
  data = df__VRT_PMON,
  vars = c("S1_FCONreapp_1", "S1_FCONreapp_3", "S1_FCONreapp_5", "S1_FCONreapp_7", "S1_FCONreapp_9", "S1_FCONreapp_11"),
  types = "PMC",
  save_means = TRUE,
  set_name = "MS1_FCONreapp",
  poly = "none"
)

df__VRT_PMON <- CenterForge(
  data = df__VRT_PMON,
  vars = c("S2_FCONreapp_1", "S2_FCONreapp_3", "S2_FCONreapp_5", "S2_FCONreapp_7", "S2_FCONreapp_9", "S2_FCONreapp_11"),
  types = "PMC",
  save_means = TRUE,
  set_name = "MS2_FCONreapp",
  poly = "none"
)

df__VRT_PMON <- CenterForge(
  data = df__VRT_PMON,
  vars = c("S1_NBSFreapp_1", "S1_NBSFreapp_3", "S1_NBSFreapp_5", "S1_NBSFreapp_7", "S1_NBSFreapp_9", "S1_NBSFreapp_11"),
  types = "PMC",
  save_means = TRUE,
  set_name = "MS1_NBSFreapp",
  poly = "none"
)

df__VRT_PMON <- CenterForge(
  data = df__VRT_PMON,
  vars = c("S2_NBSFreapp_1", "S2_NBSFreapp_3", "S2_NBSFreapp_5", "S2_NBSFreapp_7", "S2_NBSFreapp_9", "S2_NBSFreapp_11"),
  types = "PMC",
  save_means = TRUE,
  set_name = "MS2_NBSFreapp",
  poly = "none"
)

# -------------------------
# non-knotted adversity means 
# -------------------------

df__VRT_PMON <- CenterForge(
  data = df__VRT_PMON,
  vars = c("WP_FAMCONY_1", "WP_FAMCONY_3", "WP_FAMCONY_5", "WP_FAMCONY_7", "WP_FAMCONY_9", "WP_FAMCONY_11"),
  types = "PMC",
  save_means = TRUE,
  set_name = "M_FCON",
  poly = "none"
)

df__VRT_PMON <- CenterForge(
  data = df__VRT_PMON,
  vars = c("WP_mSchAdv_1", "WP_mSchAdv_3", "WP_mSchAdv_5", "WP_mSchAdv_7", "WP_mSchAdv_9"),
  types = "PMC",
  save_means = TRUE,
  set_name = "M_SchAdv",
  poly = "none"
)

df__VRT_PMON <- CenterForge(
  data = df__VRT_PMON,
  vars = c("WP_mNBdang_1", "WP_mNBdang_3", "WP_mNBdang_5", "WP_mNBdang_7", "WP_mNBdang_9", "WP_mNBdang_11"),
  types = "PMC",
  save_means = TRUE,
  set_name = "M_NBdang",
  poly = "none"
)

###############################################################################
############### REDUCE DFs
names(df__VRT_PMON)
vars_to_remove <- c(
  # LMC moderators / covariates
  "LMC_DysBpM_5", "LMC_DysBpM_7", "LMC_DysBpM_9", "LMC_DysBpM_11",
  "LMC_SysBpM_5", "LMC_SysBpM_7", "LMC_SysBpM_9", "LMC_SysBpM_11",
  "LMC_WP_mPmon_1", "LMC_WP_mPmon_3", "LMC_WP_mPmon_5", "LMC_WP_mPmon_7", "LMC_WP_mPmon_9", "LMC_WP_mPmon_11",
  "LMC_WP_mPNH_5", "LMC_WP_mPNH_7", "LMC_WP_mPNH_9",
  
  # nURG
  "S1_SCHnURG_1", "S1_SCHnURG_3", "S1_SCHnURG_5", "S1_SCHnURG_7", "S1_SCHnURG_9",
  "S2_SCHnURG_1", "S2_SCHnURG_3", "S2_SCHnURG_5", "S2_SCHnURG_7", "S2_SCHnURG_9",
  "S1_FCONnURG_1", "S1_FCONnURG_3", "S1_FCONnURG_5", "S1_FCONnURG_7", "S1_FCONnURG_9", "S1_FCONnURG_11",
  "S2_FCONnURG_1", "S2_FCONnURG_3", "S2_FCONnURG_5", "S2_FCONnURG_7", "S2_FCONnURG_9", "S2_FCONnURG_11",
  "S1_NBSFnURG_1", "S1_NBSFnURG_3", "S1_NBSFnURG_5", "S1_NBSFnURG_7", "S1_NBSFnURG_9", "S1_NBSFnURG_11",
  "S2_NBSFnURG_1", "S2_NBSFnURG_3", "S2_NBSFnURG_5", "S2_NBSFnURG_7", "S2_NBSFnURG_9", "S2_NBSFnURG_11",
  
  # pURG
  "S1_SCHpURG_1", "S1_SCHpURG_3", "S1_SCHpURG_5", "S1_SCHpURG_7", "S1_SCHpURG_9",
  "S2_SCHpURG_1", "S2_SCHpURG_3", "S2_SCHpURG_5", "S2_SCHpURG_7", "S2_SCHpURG_9",
  "S1_FCONpURG_1", "S1_FCONpURG_3", "S1_FCONpURG_5", "S1_FCONpURG_7", "S1_FCONpURG_9", "S1_FCONpURG_11",
  "S2_FCONpURG_1", "S2_FCONpURG_3", "S2_FCONpURG_5", "S2_FCONpURG_7", "S2_FCONpURG_9", "S2_FCONpURG_11",
  "S1_NBSFpURG_1", "S1_NBSFpURG_3", "S1_NBSFpURG_5", "S1_NBSFpURG_7", "S1_NBSFpURG_9", "S1_NBSFpURG_11",
  "S2_NBSFpURG_1", "S2_NBSFpURG_3", "S2_NBSFpURG_5", "S2_NBSFpURG_7", "S2_NBSFpURG_9", "S2_NBSFpURG_11",
  
  # attn
  "S1_SCHattn_1", "S1_SCHattn_3", "S1_SCHattn_5", "S1_SCHattn_7", "S1_SCHattn_9",
  "S2_SCHattn_1", "S2_SCHattn_3", "S2_SCHattn_5", "S2_SCHattn_7", "S2_SCHattn_9",
  "S1_FCONattn_1", "S1_FCONattn_3", "S1_FCONattn_5", "S1_FCONattn_7", "S1_FCONattn_9", "S1_FCONattn_11",
  "S2_FCONattn_1", "S2_FCONattn_3", "S2_FCONattn_5", "S2_FCONattn_7", "S2_FCONattn_9", "S2_FCONattn_11",
  "S1_NBSFattn_1", "S1_NBSFattn_3", "S1_NBSFattn_5", "S1_NBSFattn_7", "S1_NBSFattn_9", "S1_NBSFattn_11",
  "S2_NBSFattn_1", "S2_NBSFattn_3", "S2_NBSFattn_5", "S2_NBSFattn_7", "S2_NBSFattn_9", "S2_NBSFattn_11",
  
  # dist
  "S1_SCHdist_1", "S1_SCHdist_3", "S1_SCHdist_5", "S1_SCHdist_7", "S1_SCHdist_9",
  "S2_SCHdist_1", "S2_SCHdist_3", "S2_SCHdist_5", "S2_SCHdist_7", "S2_SCHdist_9",
  "S1_FCONdist_1", "S1_FCONdist_3", "S1_FCONdist_5", "S1_FCONdist_7", "S1_FCONdist_9", "S1_FCONdist_11",
  "S2_FCONdist_1", "S2_FCONdist_3", "S2_FCONdist_5", "S2_FCONdist_7", "S2_FCONdist_9", "S2_FCONdist_11",
  "S1_NBSFdist_1", "S1_NBSFdist_3", "S1_NBSFdist_5", "S1_NBSFdist_7", "S1_NBSFdist_9", "S1_NBSFdist_11",
  "S2_NBSFdist_1", "S2_NBSFdist_3", "S2_NBSFdist_5", "S2_NBSFdist_7", "S2_NBSFdist_9", "S2_NBSFdist_11",
  
  # emosup
  "S1_SCHemosup_1", "S1_SCHemosup_3", "S1_SCHemosup_5", "S1_SCHemosup_7", "S1_SCHemosup_9",
  "S2_SCHemosup_1", "S2_SCHemosup_3", "S2_SCHemosup_5", "S2_SCHemosup_7", "S2_SCHemosup_9",
  "S1_FCONemosup_1", "S1_FCONemosup_3", "S1_FCONemosup_5", "S1_FCONemosup_7", "S1_FCONemosup_9", "S1_FCONemosup_11",
  "S2_FCONemosup_1", "S2_FCONemosup_3", "S2_FCONemosup_5", "S2_FCONemosup_7", "S2_FCONemosup_9", "S2_FCONemosup_11",
  "S1_NBSFemosup_1", "S1_NBSFemosup_3", "S1_NBSFemosup_5", "S1_NBSFemosup_7", "S1_NBSFemosup_9", "S1_NBSFemosup_11",
  "S2_NBSFemosup_1", "S2_NBSFemosup_3", "S2_NBSFemosup_5", "S2_NBSFemosup_7", "S2_NBSFemosup_9", "S2_NBSFemosup_11",
  
  # reapp
  "S1_SCHreapp_1", "S1_SCHreapp_3", "S1_SCHreapp_5", "S1_SCHreapp_7", "S1_SCHreapp_9",
  "S2_SCHreapp_1", "S2_SCHreapp_3", "S2_SCHreapp_5", "S2_SCHreapp_7", "S2_SCHreapp_9",
  "S1_FCONreapp_1", "S1_FCONreapp_3", "S1_FCONreapp_5", "S1_FCONreapp_7", "S1_FCONreapp_9", "S1_FCONreapp_11",
  "S2_FCONreapp_1", "S2_FCONreapp_3", "S2_FCONreapp_5", "S2_FCONreapp_7", "S2_FCONreapp_9", "S2_FCONreapp_11",
  "S1_NBSFreapp_1", "S1_NBSFreapp_3", "S1_NBSFreapp_5", "S1_NBSFreapp_7", "S1_NBSFreapp_9", "S1_NBSFreapp_11",
  "S2_NBSFreapp_1", "S2_NBSFreapp_3", "S2_NBSFreapp_5", "S2_NBSFreapp_7", "S2_NBSFreapp_9", "S2_NBSFreapp_11",
  
  # attprb
  "S1_SCHattprb_1", "S1_SCHattprb_3", "S1_SCHattprb_5", "S1_SCHattprb_7", "S1_SCHattprb_9",
  "S2_SCHattprb_1", "S2_SCHattprb_3", "S2_SCHattprb_5", "S2_SCHattprb_7", "S2_SCHattprb_9",
  "S1_FCONattprb_1", "S1_FCONattprb_3", "S1_FCONattprb_5", "S1_FCONattprb_7", "S1_FCONattprb_9", "S1_FCONattprb_11",
  "S2_FCONattprb_1", "S2_FCONattprb_3", "S2_FCONattprb_5", "S2_FCONattprb_7", "S2_FCONattprb_9", "S2_FCONattprb_11",
  "S1_NBSFattprb_1", "S1_NBSFattprb_3", "S1_NBSFattprb_5", "S1_NBSFattprb_7", "S1_NBSFattprb_9", "S1_NBSFattprb_11",
  "S2_NBSFattprb_1", "S2_NBSFattprb_3", "S2_NBSFattprb_5", "S2_NBSFattprb_7", "S2_NBSFattprb_9", "S2_NBSFattprb_11"
  )

df_VRT_RED <- df__VRT_PMON[, !(names(df__VRT_PMON) %in% vars_to_remove), drop = FALSE]

## MOVE VARIABLE ORDER 

# helper 
MoveMeansToEnd <- function(data, extra_vars = NULL) {
  
  if (!is.data.frame(data)) {
    stop("data must be a data frame.")
  }
  
  mean_vars <- names(data)[grepl("^MEAN_", names(data))]
  
  if (is.null(extra_vars)) {
    extra_vars <- character(0)
  }
  
  if (!is.character(extra_vars)) {
    stop("'extra_vars' must be a character vector of variable names.")
  }
  
  bad_vars <- extra_vars[!extra_vars %in% names(data)]
  if (length(bad_vars) > 0) {
    stop("These extra_vars are not in data: ", paste(bad_vars, collapse = ", "))
  }
  
  end_vars <- unique(c(extra_vars, mean_vars))
  front_vars <- names(data)[!names(data) %in% end_vars]
  
  data <- data[, c(front_vars, end_vars), drop = FALSE]
  
  return(data)
}

# applied 
names(df_VRT_RED)

df_VRT_RED_MOVE <- MoveMeansToEnd(
  df_VRT_RED,
  extra_vars = c("S1_HSESnURG",            
                 "S2_HSESnURG",    "S1_CogEnnURG",    "S2_CogEnnURG",   "S1_HSESpURG",            
                 "S2_HSESpURG",    "S1_CogEnpURG",    "S2_CogEnpURG",   "S1_HSESattn",            
                 "S2_HSESattn",    "S1_CogEnattn",    "S2_CogEnattn",   "S1_HSESdist",            
                 "S2_HSESdist",    "S1_CogEndist",    "S2_CogEndist",   "S1_HSESemosup",          
                 "S2_HSESemosup",  "S1_CogEnemosup",  "S2_CogEnemosup",   "S1_HSESreapp",           
                 "S2_HSESreapp",   "S1_CogEnreapp",   "S2_CogEnreapp"))
names(df_VRT_RED_MOVE)
## clean names 
mean_vars <- names(df_VRT_RED_MOVE)[grepl("^MEAN_PMC_", names(df_VRT_RED_MOVE))]

names(df_VRT_RED_MOVE)[match(mean_vars, names(df_VRT_RED_MOVE))] <- 
  sub("^MEAN_PMC_", "", mean_vars)

names(df_VRT_RED_MOVE)

###############################################################################
############### SAVE DATA

## ALL VERSIONS TOGETHER 

write.csv(df_VRT_RED_MOVE, "ABCD_HORM_METH_SEM_STRUCTRUAL_4.27.26.csv", row.names = F)
prepareMplusData(df_VRT_RED_MOVE,"ABCD_HORM_METH_SEM_STRUCTRUAL_4.27.26.dat", inpfile =T)




