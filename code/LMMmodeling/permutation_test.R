rm(list=ls(all=TRUE))
here::i_am("Readme.txt")
repoDir <- here::here()
datDir <- file.path(repoDir,"data")
codeDir <- file.path(repoDir, "code")
outputDir <- file.path(repoDir, "output")
figDir <- file.path(repoDir, "figures")
tabDir <- file.path(repoDir, "tables")

library(tidyverse)
library(plyr)
library(lme4)
library(doParallel)
library(doRNG)
registerDoParallel(35)


source(file.path(codeDir, "LMMmodeling/utils.R"))
source(file.path(codeDir, "/common.R"))

nAbdata <- read.csv(file.path(datDir, "vat08_combined_data_processed_longitudinal_nAb.csv"))
bAbdata <- read.csv(file.path(datDir, "vat08_combined_data_processed_longitudinal_bAb.csv"))

bAbdata$BserostatusD01_S <- ifelse(bAbdata$Bserostatus==0, "naive", 
                                   ifelse(bAbdata$D01_S_pos_only_in_non_naive_group == 1, 
                                          "nonnaiveD01SposOnly", "nonnaiveElse"))
nAbdata$BserostatusD01_S <- ifelse(nAbdata$Bserostatus==0, "naive", 
                                   ifelse(nAbdata$D01_S_pos_only_in_non_naive_group == 1, 
                                          "nonnaiveD01SposOnly", "nonnaiveElse"))


markers <- c("bindSpike", "bindSpike_beta", "bindSpike_alpha", "bindSpike_gamma", "bindSpike_delta1", "bindSpike_delta2", "bindSpike_delta3", "bindSpike_omicron",
             "pseudoneutid50", "pseudoneutid50_B.1.351", "pseudoneutid50_BA.1", "pseudoneutid50_BA.2", "pseudoneutid50_BA.4.5")


#permutation test
nk = 6
p_value_Female <- matrix(NA, nrow = dim(contrasts)[1], ncol = length(markers)*3)
p_value_Male <- matrix(NA, nrow = dim(contrasts)[1], ncol = length(markers)*3)

D43_p_value_Female <- matrix(NA, nrow = dim(contrasts)[1], ncol = length(markers)*3)
D43_p_value_Male <- matrix(NA, nrow = dim(contrasts)[1], ncol = length(markers)*3)
D202_p_value_Female <- matrix(NA, nrow = dim(contrasts)[1], ncol = length(markers)*3)
D202_p_value_Male <- matrix(NA, nrow = dim(contrasts)[1], ncol = length(markers)*3)

#D43 to D202 rate of decrease
rate_p_value_Female <- matrix(NA, nrow = dim(contrasts)[1], ncol = length(markers)*3)
rate_p_value_Male <- matrix(NA, nrow = dim(contrasts)[1], ncol = length(markers)*3)


l = 1
load = FALSE
AUC_diff_permu_par <- list()
for(c in 1:dim(contrasts)[1]){
  GA = as.numeric(contrasts[c, 1])
  GB = as.numeric(contrasts[c, 2])
  GAinfo = filter(groups, Group==GA)
  GBinfo = filter(groups, Group==GB)
  if(GAinfo$Baseline == "nonnaive"){
    nAbdataA = filter(nAbdata, Trialstage == GAinfo$Stage & Trt == GAinfo$Treatment & Bserostatus == 1)
    bAbdataA = filter(bAbdata, Trialstage == GAinfo$Stage & Trt == GAinfo$Treatment & Bserostatus == 1)
  }else{
    nAbdataA = filter(nAbdata, Trialstage == GAinfo$Stage & Trt == GAinfo$Treatment & BserostatusD01_S == GAinfo$Baseline)
    bAbdataA = filter(bAbdata, Trialstage == GAinfo$Stage & Trt == GAinfo$Treatment & BserostatusD01_S == GAinfo$Baseline)
  }
  
  if(GBinfo$Baseline == "nonnaive"){
    nAbdataB = filter(nAbdata, Trialstage == GBinfo$Stage & Trt == GBinfo$Treatment & Bserostatus == 1)
    bAbdataB = filter(bAbdata, Trialstage == GBinfo$Stage & Trt == GBinfo$Treatment & Bserostatus == 1)
    
  }else{
    nAbdataB = filter(nAbdata, Trialstage == GBinfo$Stage & Trt == GBinfo$Treatment & BserostatusD01_S == GBinfo$Baseline)
    bAbdataB = filter(bAbdata, Trialstage == GBinfo$Stage & Trt == GBinfo$Treatment & BserostatusD01_S == GBinfo$Baseline)
  }
  
  if(GA %in% c(1, 4) | GB %in% c(1, 4)){ #naive groups
    baselineAdj = FALSE
  }else{
    baselineAdj = TRUE
  }
  j = 1
  for(marki in markers){
    if(marki %in% c("pseudoneutid50", "pseudoneutid50_B.1.351", "pseudoneutid50_BA.1", "pseudoneutid50_BA.2", "pseudoneutid50_BA.4.5")){
      dataA <- nAbdataA
      dataB <- nAbdataB
    }else{
      dataA <- bAbdataA
      dataB <- bAbdataB
    }
    dataA$mark <- dataA[, marki]
    dataB$mark <- dataB[, marki]
    dataA$Bmark <- dataA[, paste0("B", marki)]
    dataB$Bmark <- dataB[, paste0("B", marki)]
    
    #remove observations with missing mark
    dataAi <- filter(dataA, !is.na(mark))
    dataBi <- filter(dataB, !is.na(mark))
    #remove participants with less than 3 obs
    nobsA <- table(dataAi$Ptid)
    ptidstoIncludeA <- names(nobsA)[nobsA >= 3]
    nobsB <- table(dataBi$Ptid)
    ptidstoIncludeB <- names(nobsB)[nobsB >= 3]
    dataAi<- filter(dataAi, Ptid%in%ptidstoIncludeA)
    dataBi<- filter(dataBi, Ptid%in%ptidstoIncludeB)
    
    dataPooli <- rbind(dataAi, dataBi)
    
    
    fitA <- BFpolynomialFit(dataAi, baselineAdj)
    AUCcA <- AUCsummaryConditional(fitA$fmsummary, fitA$p)
    
    fitB <- BFpolynomialFit(dataBi, baselineAdj)
    AUCcB <- AUCsummaryConditional(fitB$fmsummary, fitB$p)
    
    AUC_diff_Female = AUCcA$est[1] - AUCcB$est[1]
    AUC_diff_Male = AUCcA$est[3] - AUCcB$est[3]
    
    D43fitA <- AntibodyLevelsummaryConditional (fitA$fmsummary, fitA$p, 45)
    D202fitA <- AntibodyLevelsummaryConditional (fitA$fmsummary, fitA$p, 204)
    D43fitB <- AntibodyLevelsummaryConditional (fitB$fmsummary, fitB$p, 45)
    D202fitB <- AntibodyLevelsummaryConditional (fitB$fmsummary, fitB$p, 204)
    
    D43_diff_main <- D43fitA[c(1,3)] - D43fitB[c(1,3)]
    D202_diff_main <- D202fitA[c(1,3)] - D202fitB[c(1,3)]
    
    rate_diff_main_A <- D202fitA[c(1,3)]-D43fitA[c(1,3)]
    rate_diff_main_B <- D202fitB[c(1,3)]-D43fitB[c(1,3)]
    rate_diff_main <- rate_diff_main_A - rate_diff_main_B
    
    uniqueptids <- unique(dataPooli$Ptid)
    if(!load){
      AUC_diff_permu_all <- foreach(k = 1:nk, .options.RNG = 1234)%dorng%{
        ptidsk <- sample(uniqueptids)
        nptidA <- length(unique(dataAi$Ptid))
        nptidB <- length(unique(dataBi$Ptid))
        dataAk <- filter(dataPooli, Ptid %in% ptidsk[1:nptidA])
        dataBk <- filter(dataPooli, !Ptid %in% ptidsk[1:nptidA])
        fitA <- BFpolynomialFit(dataAk, baselineAdj)
        AUCcA <- AUCsummaryConditional(fitA$fmsummary, fitA$p)
        
        fitB <- BFpolynomialFit(dataBk, baselineAdj)
        AUCcB <- AUCsummaryConditional(fitB$fmsummary, fitB$p)
        
        D43fitA <- AntibodyLevelsummaryConditional (fitA$fmsummary, fitA$p, 45)
        D202fitA <- AntibodyLevelsummaryConditional (fitA$fmsummary, fitA$p, 204)
        D43fitB <- AntibodyLevelsummaryConditional (fitB$fmsummary, fitB$p, 45)
        D202fitB <- AntibodyLevelsummaryConditional (fitB$fmsummary, fitB$p, 204)
        
        D43diff <- NULL
        D43diff[1:2] <- D43fitA[c(1,3)] - D43fitB[c(1,3)]
        D202diff <- NULL
        D202diff[1:2] <- D202fitA[c(1,3)] - D202fitB[c(1,3)]
        
        rate_A <- D202fitA[c(1,3)]-D43fitA[c(1,3)]
        rate_B <- D202fitB[c(1,3)]-D43fitB[c(1,3)]
        
        ans <- c(AUCcA$est[1] - AUCcB$est[1], AUCcA$est[3] - AUCcB$est[3], fitA$p, fitB$p,
                 D43diff, D202diff, rate_A-rate_B)
        names(ans) <- c("AUC_diff_Female", "AUC_diff_Male", "pA", "pB", "D43_diff_Female", "D43_diff_Male", 
                        "D202_diff_Female", "D202_diff_Male", "rate_diff_Female", "rate_diff_Male")
        return(ans)
      }
      AUC_diff_permu_par[[l]] <- laply(AUC_diff_permu_all, function(x)c(x[1:10]))
      AUC_diff_permu_Female <- laply(AUC_diff_permu_all, function(x)x[1])
      AUC_diff_permu_Male <- laply(AUC_diff_permu_all, function(x)x[2])
      
      D43_diff_permu_Female <- laply(AUC_diff_permu_all, function(x)x[5])
      D43_diff_permu_Male <- laply(AUC_diff_permu_all, function(x)x[6])
      D202_diff_permu_Female <- laply(AUC_diff_permu_all, function(x)x[7])
      D202_diff_permu_Male <- laply(AUC_diff_permu_all, function(x)x[8])
      rate_diff_permu_Female <- laply(AUC_diff_permu_all, function(x)x[9])
      rate_diff_permu_Male <- laply(AUC_diff_permu_all, function(x)x[10])
    }else{
      AUC_diff_permu_all <- readRDS(file.path(outputDir, "permu_par.rds"))[[l]]
      AUC_diff_permu_Female <- AUC_diff_permu_all[, 1]
      AUC_diff_permu_Male <- AUC_diff_permu_all[, 2]
      
      D43_diff_permu_Female <- AUC_diff_permu_all[, 5]
      D43_diff_permu_Male <- AUC_diff_permu_all[, 5]
      D202_diff_permu_Female <- AUC_diff_permu_all[, 7]
      D202_diff_permu_Male <- AUC_diff_permu_all[, 8]
      rate_diff_permu_Female <- AUC_diff_permu_all[, 9]
      rate_diff_permu_Male <- AUC_diff_permu_all[, 10]
    }
    
    
    
    l = l + 1
    p_value_Female[c, j] <- sum(abs(AUC_diff_permu_Female) > abs(AUC_diff_Female))/nk
    p_value_Male[c, j] <- sum(abs(AUC_diff_permu_Male) > abs(AUC_diff_Male))/nk
    p_value_Female[c, j+13] <- sum(AUC_diff_permu_Female > AUC_diff_Female)/nk
    p_value_Male[c, j+13] <- sum(AUC_diff_permu_Male > AUC_diff_Male)/nk
    p_value_Female[c, j+26] <- sum(AUC_diff_permu_Female < AUC_diff_Female)/nk
    p_value_Male[c, j+26] <- sum(AUC_diff_permu_Male < AUC_diff_Male)/nk
    
    D43_p_value_Female[c, j] <- sum(abs(D43_diff_permu_Female) > abs(D43_diff_main[1]))/nk
    D43_p_value_Male[c, j] <- sum(abs(D43_diff_permu_Male) > abs(D43_diff_main[2]))/nk
    D43_p_value_Female[c, j+13] <- sum(D43_diff_permu_Female > D43_diff_main[1])/nk
    D43_p_value_Male[c, j+13] <- sum(D43_diff_permu_Male > D43_diff_main[2])/nk
    D43_p_value_Female[c, j+26] <- sum(D43_diff_permu_Female < D43_diff_main[1])/nk
    D43_p_value_Male[c, j+26] <- sum(D43_diff_permu_Male < D43_diff_main[2])/nk
    
    
    D202_p_value_Female[c, j] <- sum(abs(D202_diff_permu_Female) > abs(D202_diff_main[1]))/nk
    D202_p_value_Male[c, j] <- sum(abs(D202_diff_permu_Male) > abs(D202_diff_main[2]))/nk
    D202_p_value_Female[c, j+13] <- sum(D202_diff_permu_Female > D202_diff_main[1])/nk
    D202_p_value_Male[c, j+13] <- sum(D202_diff_permu_Male > D202_diff_main[2])/nk
    D202_p_value_Female[c, j+26] <- sum(D202_diff_permu_Female < D202_diff_main[1])/nk
    D202_p_value_Male[c, j+26] <- sum(D202_diff_permu_Male < D202_diff_main[2])/nk
    
    rate_p_value_Female[c, j] <- sum(abs(rate_diff_permu_Female) > abs(rate_diff_main[1]))/nk
    rate_p_value_Male[c, j] <- sum(abs(rate_diff_permu_Male) > abs(rate_diff_main[2]))/nk
    rate_p_value_Female[c, j+13] <- sum(rate_diff_permu_Female > rate_diff_main[1])/nk
    rate_p_value_Male[c, j+13] <- sum(rate_diff_permu_Male > rate_diff_main[2])/nk
    rate_p_value_Female[c, j+26] <- sum(rate_diff_permu_Female < rate_diff_main[1])/nk
    rate_p_value_Male[c, j+26] <- sum(rate_diff_permu_Male < rate_diff_main[2])/nk
    
    j = j + 1
  }
  
}
colnames(p_value_Female)  <- c(markers, paste(markers, "AgtB", sep = "_"), paste(markers, "AltB", sep = "_"))
colnames(p_value_Male)  <- c(markers, paste(markers, "AgtB", sep = "_"), paste(markers, "AltB", sep = "_"))

colnames(D43_p_value_Female)  <- c(markers, paste(markers, "AgtB", sep = "_"), paste(markers, "AltB", sep = "_"))
colnames(D43_p_value_Male)  <- c(markers, paste(markers, "AgtB", sep = "_"), paste(markers, "AltB", sep = "_"))
colnames(D202_p_value_Female) <- c(markers, paste(markers, "AgtB", sep = "_"), paste(markers, "AltB", sep = "_"))
colnames(D202_p_value_Male)  <- c(markers, paste(markers, "AgtB", sep = "_"), paste(markers, "AltB", sep = "_"))

colnames(rate_p_value_Female)  <- c(markers, paste(markers, "AgtB", sep = "_"), paste(markers, "AltB", sep = "_"))
colnames(rate_p_value_Male)  <- c(markers, paste(markers, "AgtB", sep = "_"), paste(markers, "AltB", sep = "_"))

write.csv(p_value_Female, file.path(tabDir, "AUC_p_value_Female.csv"))
write.csv(p_value_Male, file.path(tabDir, "AUC_p_value_Male.csv"))
write.csv(D43_p_value_Female, file.path(tabDir, "D43_p_value_Female.csv"))
write.csv(D43_p_value_Male, file.path(tabDir, "D43_p_value_Male.csv"))
write.csv(D202_p_value_Female, file.path(tabDir, "D202_p_value_Female.csv"))
write.csv(D202_p_value_Male, file.path(tabDir, "D202_p_value_Male.csv"))
write.csv(rate_p_value_Female, file.path(tabDir, "D202toD43GMR_p_value_Female.csv"))
write.csv(rate_p_value_Male, file.path(tabDir, "D202toD43GMR_p_value_Male.csv"))

if(!load){
  saveRDS(AUC_diff_permu_par, file.path(outputDir, "permu_par.rds"))
}








