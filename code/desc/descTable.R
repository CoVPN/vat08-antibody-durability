rm(list=ls(all=TRUE))
here::i_am("README.md")
repoDir <- here::here()
datDir <- file.path(repoDir,"data")
codeDir <- file.path(repoDir, "code")
outputDir <- file.path(repoDir, "output")
figDir <- file.path(repoDir, "figures")
tabDir <- file.path(repoDir, "tables")
  
library(tidyverse)
library(plyr)
source(file.path(codeDir, "common.R"))

#Number of participants with antibody data in the longitudinal analysis cohorts
nAbdata <- read.csv(file.path(datDir, "vat08_combined_data_processed_longitudinal_nAb.csv"))
bAbdata <- read.csv(file.path(datDir, "vat08_combined_data_processed_longitudinal_bAb.csv"))

numberOfObsTable_nAb <- tibble("Stage" = character(), "Group" = character(), "Treatment" = character(),
                           "D43" = numeric(), "D78" = numeric(), "D134" = numeric(), "D202" = numeric(),
                           "D292" = numeric(), "D387" = numeric())
numberOfObsTable_bAb <- tibble("Stage" = character(), "Group" = character(), "Treatment" = character(),
                               "D43" = numeric(), "D78" = numeric(), "D134" = numeric(), "D202" = numeric(),
                               "D292" = numeric(), "D387" = numeric())
for(s in c(1, 2)){
  for(g in c("Naive", "Non-naive", "Non-naive Group A", "Non-naive Group B")){
    for(treatment in c(1, 0)){
      if(g == "Naive"){
        nAbdata_subset <- filter(nAbdata, Trialstage == s & Bserostatus == 0 & Trt == treatment)
        bAbdata_subset <- filter(bAbdata, Trialstage == s & Bserostatus == 0 & Trt == treatment)
      }
      if(g == "Non-naive"){
        nAbdata_subset <- filter(nAbdata, Trialstage == s & Bserostatus == 1 & Trt == treatment)
        bAbdata_subset <- filter(bAbdata, Trialstage == s & Bserostatus == 1 & Trt == treatment)
      }
      
      if(g == "Non-naive Group A"){
        nAbdata_subset <- filter(nAbdata, Trialstage == s & Bserostatus == 1 & D01_S_pos_only_in_non_naive_group == 1 & Trt == treatment)
        bAbdata_subset <- filter(bAbdata, Trialstage == s & Bserostatus == 1 & D01_S_pos_only_in_non_naive_group == 1 & Trt == treatment)
      }
      
      if(g == "Non-naive Group B"){
        nAbdata_subset <- filter(nAbdata, Trialstage == s & Bserostatus == 1 & D01_S_pos_only_in_non_naive_group == 0 & Trt == treatment)
        bAbdata_subset <- filter(bAbdata, Trialstage == s & Bserostatus == 1 & D01_S_pos_only_in_non_naive_group == 0 & Trt == treatment)
      }
      
      time <- c(43, 78, 134, 202, 292, 387)
      if(dim(nAbdata_subset)[1] >0 ){
        nobs_nAb <- laply(time, function(x){
          tmp <- filter(nAbdata_subset, time2 == x)
          return(ifelse(dim(tmp)[1] >0, dim(tmp)[1], 0))
        })
      }else{
        nobs_nAb <- rep(0, 6)
      }
      
      if(dim(bAbdata_subset)[1] >0 ){
        nobs_bAb <- laply(time, function(x){
          tmp <- filter(bAbdata_subset, time2 == x)
          return(ifelse(dim(tmp)[1] >0, dim(tmp)[1], 0))
        })
      }else{
        nobs_bAb <- rep(0, 6)
      }
      
      numberOfObsTable_nAb <- add_row(numberOfObsTable_nAb , "Stage" = ifelse(s==1, "Stage 1", "Stage 2"), "Group" = g, 
                                      "Treatment" = ifelse(treatment == 1, "Vaccine", "Placebo"),
                                     "D43" = nobs_nAb[1], "D78" = nobs_nAb[2], "D134" = nobs_nAb[3], "D202" = nobs_nAb[4],
                                     "D292" = nobs_nAb[5], "D387" = nobs_nAb[6])
      numberOfObsTable_bAb <- add_row(numberOfObsTable_bAb , "Stage" = ifelse(s==1, "Stage 1", "Stage 2"), "Group" = g, 
                                      "Treatment" = ifelse(treatment == 1, "Vaccine", "Placebo"),
                                      "D43" = nobs_bAb[1], "D78" = nobs_bAb[2], "D134" = nobs_bAb[3], "D202" = nobs_bAb[4],
                                      "D292" = nobs_bAb[5], "D387" = nobs_bAb[6])
    }
   
  }
}

write.csv(numberOfObsTable_nAb, file.path(tabDir, "numberOfObsTable_nAb.csv"))
write.csv(numberOfObsTable_bAb, file.path(tabDir, "numberOfObsTable_bAb.csv"))


#Number of participants in the longitudinal sub-cohort with antibody data after imputation
#longitudinal sub-cohort refers to participants with antibody marker data post D78
df1 <- read.csv(file.path(datDir, "antibodyDurabilityWide_imputed.csv"))
longdf <- read.csv(file.path(datDir, "antibodyDurabilityLong_imputed.csv"))
fulllongdf <- left_join(longdf, df1, by = "Ptid")

postD78df <- filter(longdf, !visitn%in%c("Baseline", "Day22", "Day43", "Day78"))
postD78df  <- postD78df %>%filter(rowSums(is.na(.)) != 13)

fulllongdf <- filter(fulllongdf, Ptid %in% unique(postD78df$Ptid))
write.csv(unique(postD78df$Ptid), file.path(datDir, "longitudinalSubcohortPtid.csv"))


fulllongdf$time2 <- as.numeric(ifelse(fulllongdf$visitn == "Baseline", 0, gsub("Day","",fulllongdf$visitn)))
numberOfObsTable_nAb <- tibble("Stage" = character(), "Group" = character(), "Treatment" = character(),
                               "D1" = numeric(), "D22" = numeric(),
                               "D43" = numeric(), "D78" = numeric(), "D134" = numeric(), "D202" = numeric(),
                               "D292" = numeric(), "D387" = numeric())
numberOfObsTable_bAb <- tibble("Stage" = character(), "Group" = character(), "Treatment" = character(),
                               "D1" = numeric(), "D22" = numeric(),
                               "D43" = numeric(), "D78" = numeric(), "D134" = numeric(), "D202" = numeric(),
                               "D292" = numeric(), "D387" = numeric())
for(s in c(1, 2)){
  for(g in c("Naive", "Non-naive", "Non-naive Group A", "Non-naive Group B")){
    for(treatment in c(1, 0)){
      if(g == "Naive"){
        nAbdata_subset <- filter(fulllongdf, Trialstage == s & Bserostatus == 0 & Trt == treatment)
        bAbdata_subset <- filter(fulllongdf, Trialstage == s & Bserostatus == 0 & Trt == treatment)
      }
      if(g == "Non-naive"){
        nAbdata_subset <- filter(fulllongdf, Trialstage == s & Bserostatus == 1 & Trt == treatment)
        bAbdata_subset <- filter(fulllongdf, Trialstage == s & Bserostatus == 1 & Trt == treatment)
      }
      
      if(g == "Non-naive Group A"){
        nAbdata_subset <- filter(fulllongdf, Trialstage == s & Bserostatus == 1 & D01_S_pos_only_in_non_naive_group == 1 & Trt == treatment)
        bAbdata_subset <- filter(fulllongdf, Trialstage == s & Bserostatus == 1 & D01_S_pos_only_in_non_naive_group == 1 & Trt == treatment)
      }
      
      if(g == "Non-naive Group B"){
        nAbdata_subset <- filter(fulllongdf, Trialstage == s & Bserostatus == 1 & D01_S_pos_only_in_non_naive_group == 0 & Trt == treatment)
        bAbdata_subset <- filter(fulllongdf, Trialstage == s & Bserostatus == 1 & D01_S_pos_only_in_non_naive_group == 0 & Trt == treatment)
      }
      
      time <- c(0, 22, 43, 78, 134, 202, 292, 387)
      if(dim(nAbdata_subset)[1] >0 ){
        nobs_nAb <- laply(time, function(x){
          tmp <- filter(nAbdata_subset, time2 == x &!is.na(pseudoneutid50))
          return(ifelse(dim(tmp)[1] >0, dim(tmp)[1], 0))
        })
      }else{
        nobs_nAb <- rep(0, 8)
      }
      
      if(dim(bAbdata_subset)[1] >0 ){
        nobs_bAb <- laply(time, function(x){
          tmp <- filter(bAbdata_subset, time2 == x &!is.na(bindSpike))
          return(ifelse(dim(tmp)[1] >0, dim(tmp)[1], 0))
        })
      }else{
        nobs_bAb <- rep(0, 8)
      }
      
      numberOfObsTable_nAb <- add_row(numberOfObsTable_nAb , "Stage" = ifelse(s==1, "Stage 1", "Stage 2"), "Group" = g, 
                                      "Treatment" = ifelse(treatment == 1, "Vaccine", "Placebo"),
                                      "D1" = nobs_nAb[1], "D22" = nobs_nAb[2],
                                      "D43" = nobs_nAb[3], "D78" = nobs_nAb[4], "D134" = nobs_nAb[5], "D202" = nobs_nAb[6],
                                      "D292" = nobs_nAb[7], "D387" = nobs_nAb[8])
      numberOfObsTable_bAb <- add_row(numberOfObsTable_bAb , "Stage" = ifelse(s==1, "Stage 1", "Stage 2"), "Group" = g, 
                                      "Treatment" = ifelse(treatment == 1, "Vaccine", "Placebo"),
                                      "D1" = nobs_bAb[1], "D22" = nobs_bAb[2],
                                      "D43" = nobs_bAb[3], "D78" = nobs_bAb[4], "D134" = nobs_bAb[5], "D202" = nobs_bAb[6],
                                      "D292" = nobs_bAb[7], "D387" = nobs_bAb[8])
    }
    
  }
}

write.csv(numberOfObsTable_nAb, file.path(tabDir, "numberOfObsTable_nAb_imputed_longitudinalSubcohort.csv"))
write.csv(numberOfObsTable_bAb, file.path(tabDir, "numberOfObsTable_bAb_imputed_longitudinalSubcohort.csv"))

#percentage of nAb had occasional missing
missingTab <- read.csv(file.path(datDir, "missingTab.csv"))

nAbmissingTab <- filter(missingTab, assay == "nAb")
bAbmissingTab <- filter(missingTab, assay == "bAb")


#Total number of serum samples for participants in the per-protocol analysis set (PPAS) with measured 
#nAb-ID50 titers against at least one antigen, number of samples (percentages) with occasional missingness 
#(missing titers against 1–3 antigens), and number of samples (percentages) with non-occasional missingness 
#(missing titers against 4 antigens) at each scheduled visit.
nAbmissingTab <- ddply(nAbmissingTab, .(time), function(df){
  total = sum(df$total)
  occasionalMiss = sum(df$occasionalMiss)
  oneMarkerObserved = sum(df$oneMarkerObserved)
  return(c("total" = total, "occasionalMiss" = occasionalMiss, "oneMarkerObserved" = oneMarkerObserved))
})
nAbmissingTab <- data.frame(nAbmissingTab)
nAbmissingTab$ocMisP <- round(nAbmissingTab$occasionalMiss/nAbmissingTab$total*100, 1)
nAbmissingTab$oneP <- round(nAbmissingTab$oneMarkerObserved/nAbmissingTab$total*100, 1)
nAbmissingTab$occasionalMiss <- paste0(nAbmissingTab$occasionalMiss, " (", nAbmissingTab$ocMisP, "%",")")
nAbmissingTab$oneMarkerObserved <- paste0(nAbmissingTab$oneMarkerObserved, " (", nAbmissingTab$oneP, "%",")")
nAbmissingTab$ocMisP <- NULL
nAbmissingTab$oneP <- NULL
#Total number of samples (blood samples) for participants in the per-protocol analysis set (PPAS) with measured
#bAb-IgG concentration against at least one antigen, number of samples (percentages) with occasional missingness
#(missing concentrations against 1–6 antigens), and number of samples (percentages) with non-occasional missingness
#(missing concentrations against 7 antigens) at each scheduled visit.
bAbmissingTab <- ddply(bAbmissingTab, .(time), function(df){
  total = sum(df$total)
  occasionalMiss = sum(df$occasionalMiss)
  oneMarkerObserved = sum(df$oneMarkerObserved)
  return(c("total" = total, "occasionalMiss" = occasionalMiss, "oneMarkerObserved" = oneMarkerObserved))
})

bAbmissingTab <- data.frame(bAbmissingTab)
bAbmissingTab$ocMisP <- round(bAbmissingTab$occasionalMiss/bAbmissingTab$total*100, 1)
bAbmissingTab$oneP <- round(bAbmissingTab$oneMarkerObserved/bAbmissingTab$total*100, 1)
bAbmissingTab$occasionalMiss <- paste0(bAbmissingTab$occasionalMiss, " (", bAbmissingTab$ocMisP, "%",")")
bAbmissingTab$oneMarkerObserved <- paste0(bAbmissingTab$oneMarkerObserved, " (", bAbmissingTab$oneP, "%",")")
bAbmissingTab$ocMisP <- NULL
bAbmissingTab$oneP <- NULL

write.csv(nAbmissingTab, file.path(tabDir, "nAbmissing.csv"))
write.csv(bAbmissingTab, file.path(tabDir, "bAbmissing.csv"))

#D1 antibody marker level that were conditioned on
nAbdata0 <- filter(nAbdata, Bserostatus == 1)
bAbdata0 <- filter(bAbdata, Bserostatus == 1) 

vars1 <- paste("B", c("pseudoneutid50", "pseudoneutid50_BA.1", "pseudoneutid50_BA.2", "pseudoneutid50_BA.4.5",
                      "pseudoneutid50_B.1.351"), "_0", sep = "")

vars2 <- paste("B", c("bindSpike", "bindSpike_beta", "bindSpike_alpha", "bindSpike_gamma", "bindSpike_delta1", "bindSpike_delta2",
                      "bindSpike_delta3", "bindSpike_omicron"), "_0", sep = "")

D1mean <- tibble("marker" = character(), "meanlevel" = numeric())
for(var in c(vars1, vars2)){
  if(var %in% vars1){
    df <- distinct(nAbdata0[, c("Ptid", var)])
    markermean <- mean(df[,var], na.rm = TRUE)
  }else{
    df <- distinct(bAbdata0[, c("Ptid", var)])
    markermean <- mean(df[,var], na.rm = TRUE)
  }
  D1mean <- add_row(.data = D1mean, "marker" = var, "meanlevel" = round(markermean, 2))
}
write.csv(D1mean, file.path(tabDir, "meanlog10D1antibodylevel.csv"))
