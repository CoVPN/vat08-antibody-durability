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
source(file.path(codeDir, "common.R"))

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


modelsFitted <- data.frame()

nk = 6
set.seed(2)
runBootPar = TRUE

if(runBootPar){
  AUC_boot_par <- list()
  l = 1
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
      
     
      AUC_boot_all <- foreach(k = 1:nk, .options.RNG = 1234)%dorng%{
        uniqueptidsA <- unique(dataAi$Ptid)
        sampleptidA <- sample(uniqueptidsA, replace = TRUE)
        uniqueptidsB <- unique(dataBi$Ptid)
        sampleptidB <- sample(uniqueptidsB, replace = TRUE)
        dataAk <- data.frame()
        #create bootstrap samples
        Ptid <- 1
        l_ply(sampleptidA, function(x){
          datai <- filter(dataAi, Ptid%in%x)
          datai$Ptid <- Ptid
          dataAk <<- rbind(dataAk, datai)
          Ptid <<- Ptid + 1
        })
        dataBk <- data.frame()
        Ptid <- 1
        l_ply(sampleptidB, function(x){
          datai <- filter(dataBi, Ptid%in%x)
          datai$Ptid <- Ptid
          dataBk <<- rbind(dataBk, datai)
          Ptid <<- Ptid + 1
        })
        
        
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
        D43diff[3:4] <- D43fitA[c(2,4)] + D43fitB[c(2,4)]
        
        D202diff <- NULL
        D202diff[1:2] <- D202fitA[c(1,3)] - D202fitB[c(1,3)]
        D202diff[3:4] <- D202fitA[c(2,4)] + D202fitB[c(2,4)]
        
        rate_A <- D202fitA[c(1,3)]-D43fitA[c(1,3)]
        rate_B <- D202fitB[c(1,3)]-D43fitB[c(1,3)]
        
        ans <- c(est_Female = AUCcA$est[1] - AUCcB$est[1], var_Female = AUCcA$est[2] + AUCcB$est[2], 
                 est_Male = AUCcA$est[3] - AUCcB$est[3], var_Male = AUCcA$est[4] + AUCcB$est[4], 
                 pA = fitA$p, pB = fitB$p, 
                 nptidA = fitA$nptid, nptidB = fitB$nptid,
                 AUCcA$est[1:4], 
                 AUCcB$est[1:4],
                 D43fitA[1:4],
                 D202fitA[1:4], 
                 D43fitB[1:4],
                 D202fitB[1:4],
                 D43diff[1:4],
                 D202diff[1:4],
                 rate_A,
                 rate_B,
                 rate_A - rate_B
        )
        names(ans) <- c("diff_Female_est", "diff_Female_var", "diff_Male_est", "diff_Male_var",
                        "pA", "pB", "nptidA", "nptidB",
                        "AUC_Female_est_GA", "AUC_Female_var_GA", "AUC_Male_est_GA", "AUC_Male_var_GA",
                        "AUC_Female_est_GB", "AUC_Female_var_GB", "AUC_Male_est_GB", "AUC_Male_var_GB",
                        "D43_Female_est_GA", "D43_Female_var_GA", "D43_Male_est_GA", "D43_Male_var_GA", 
                        "D202_Female_est_GA", "D202_Female_var_GA", "D202_Male_est_GA", "D202_Male_var_GA", 
                        "D43_Female_est_GB", "D43_Female_var_GB", "D43_Male_est_GB", "D43_Male_var_GB",
                        "D202_Female_est_GB", "D202_Female_var_GB", "D202_Male_est_GB", "D202_Male_var_GB",
                        "D43_Female_est_diff", "D43_Male_est_diff", "D43_Female_var_diff", "D43_Male_var_diff",
                        "D202_Female_est_diff", "D202_Male_est_diff", "D202_Female_var_diff", "D202_Male_var_diff",
                        "rate_Female_GA",  "rate_Male_GA",
                        "rate_Female_GB", "rate_Male_GB",
                        "rate_diff_Female", "rate_diff_Male"
                        )
        return(ans)
      }
      AUC_boot_par[[l]] <- laply(AUC_boot_all, function(x)x[1:46])
      l = l + 1
    }}
  
  saveRDS(AUC_boot_par, file.path(outputDir, "boot_par.rds"))
  
}else{
  AUC_boot_par <- readRDS(file.path(outputDir, "boot_par.rds"))
}



AUC_Female_CL <- matrix(NA, nrow = dim(contrasts)[1], ncol = length(markers)*12)
AUC_Female_CU <- matrix(NA, nrow = dim(contrasts)[1], ncol = length(markers)*12)
AUC_Female <- matrix(NA, nrow = dim(contrasts)[1], ncol = length(markers)*12)
n_Female <- matrix(NA, nrow = dim(contrasts)[1], ncol = length(markers)*2)

AUC_Male_CL <- matrix(NA, nrow = dim(contrasts)[1], ncol = length(markers)*12)
AUC_Male_CU <- matrix(NA, nrow = dim(contrasts)[1], ncol = length(markers)*12)
AUC_Male <- matrix(NA, nrow = dim(contrasts)[1], ncol = length(markers)*12)
#confidence intervals
alpha = 0.05

l = 1
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
    
    
    
    fitA <- BFpolynomialFit(dataAi, baselineAdj)
    AUCcA <- AUCsummaryConditional(fitA$fmsummary, fitA$p)
   
    fitB <- BFpolynomialFit(dataBi, baselineAdj)
    AUCcB <- AUCsummaryConditional(fitB$fmsummary, fitB$p)
   
    D43fitA <- AntibodyLevelsummaryConditional (fitA$fmsummary, fitA$p, 45)
    D202fitA <- AntibodyLevelsummaryConditional (fitA$fmsummary, fitA$p, 204)
    D43fitB <- AntibodyLevelsummaryConditional (fitB$fmsummary, fitB$p, 45)
    D202fitB <- AntibodyLevelsummaryConditional (fitB$fmsummary, fitB$p, 204)
    
    
    
    AUC_Female_i = AUCcA$est[1] - AUCcB$est[1]
    AUC_Male_i = AUCcA$est[3] - AUCcB$est[3]
   
    AUC_Female_var = AUCcA$est[2] + AUCcB$est[2]
    AUC_Male_var = AUCcA$est[4] + AUCcB$est[4]
    
    D43diff <- NULL
    D43diff[1:2] <- D43fitA[c(1,3)] - D43fitB[c(1,3)]
    D43diff[3:4] <- D43fitA[c(2,4)] + D43fitB[c(2,4)]
    
    D202diff <- NULL
    D202diff[1:2] <- D202fitA[c(1,3)] - D202fitB[c(1,3)]
    D202diff[3:4] <- D202fitA[c(2,4)] + D202fitB[c(2,4)]
    
    rate_A <- D202fitA[c(1,3)]-D43fitA[c(1,3)]
    rate_B <- D202fitB[c(1,3)]-D43fitB[c(1,3)]
    
    AUC_boot_all <- AUC_boot_par[[l]]
    
    #female diff, female group A, female group B 
    AUC_boot_Female <- AUC_boot_all[, c("diff_Female_est", "AUC_Female_est_GA", "AUC_Female_est_GB", "D43_Female_est_GA", "D43_Female_est_GB", 
                                        "D202_Female_est_GA", "D202_Female_est_GB", "D43_Female_est_diff", "D202_Female_est_diff",
                                        "rate_Female_GA", "rate_Female_GB", "rate_diff_Female")]
    bkstar <- AUC_boot_Female
    b <- c(AUC_Female_i[1], AUCcA$est[1], AUCcB$est[1], D43fitA["Est_Female"], D43fitB["Est_Female"], D202fitA["Est_Female"], D202fitB["Est_Female"],
           D43diff[1], D202diff[1], rate_A[1], rate_B[1], rate_A[1] - rate_B[1])
    
   
    colsc <- c(j, 13+j, 26+j, 39+j, 52+j, 65+j, 78+j, 91+j, 104+j, 117 + j, 130+j, 143+j ) 
    
    aucLL <- apply(bkstar, 2, function(x)quantile(x, probs = 0.025))
    aucUL <- apply(bkstar, 2, function(x)quantile(x, probs = 0.975))
    
 
    AUC_Female[c, colsc] = as.vector(b)
    AUC_Female_CL[c, colsc] = aucLL
    AUC_Female_CU[c, colsc] = aucUL
    
    n_Female[c, j] <- length(ptidstoIncludeA)
    n_Female[c, j+13] <- length(ptidstoIncludeB)
    
    #male 
    AUC_boot_Male <- AUC_boot_all[, c("diff_Female_est", "AUC_Female_est_GA", "AUC_Female_est_GB", "D43_Female_est_GA", "D43_Female_est_GB", 
                                      "D202_Female_est_GA", "D202_Female_est_GB", "D43_Female_est_diff", "D202_Female_est_diff",
                                      "rate_Female_GA", "rate_Female_GB", "rate_diff_Female")]
    
    b <- c(AUC_Male_i[1], AUCcA$est[1], AUCcB$est[1], D43fitA["Est_Male"], D43fitB["Est_Male"], D202fitA["Est_Male"], D202fitB["Est_Male"],
           D43diff[2], D202diff[2], rate_A[2], rate_B[2], rate_A[2] - rate_B[2])
   
    colsc <- c(j, 13+j, 26+j, 39+j, 52+j, 65+j, 78+j, 91+j, 104+j, 117 + j, 130+j, 143+j ) 
    
    aucLL <- apply(bkstar, 2, function(x)quantile(x, probs = 0.025))
    aucUL <- apply(bkstar, 2, function(x)quantile(x, probs = 0.975))
    
    AUC_Male_CL[c, colsc] = aucLL
    AUC_Male_CU[c, colsc] = aucUL
    AUC_Male[c, colsc] = b
    
    
   
    j = j + 1
    l = l+1
  }
  
}


colnames(AUC_Female_CL) <- c(paste(markers, "diff", sep = "_"), paste(markers, "GA", sep = "_"), paste(markers, "GB", sep = "_"),
                             paste(markers, "D43_GA", sep = "_"), paste(markers, "D43_GB", sep = "_"), 
                             paste(markers, "D202_GA", sep = "_"), paste(markers, "D202_GB", sep = "_"),
                             paste(markers, "D43_diff", sep = "_"), paste(markers, "D202_diff", sep = "_"),
                             paste(markers, "rate_GA", sep = "_"), paste(markers, "rate_GB", sep = "_"),
                             paste(markers, "rate_diff", sep = "_"))
colnames(AUC_Female_CU) <- colnames(AUC_Female_CL)
colnames(AUC_Female) <- colnames(AUC_Female_CL)
colnames(n_Female) <- c(paste(markers, "A", sep = "_"), paste(markers, "B", sep = "_"))

colnames(AUC_Male_CL) <- colnames(AUC_Female_CL)
colnames(AUC_Male_CU) <- colnames(AUC_Female_CL)
colnames(AUC_Male) <- colnames(AUC_Female_CL)

write.csv(AUC_Female_CL, file.path(tabDir, "est_Female_CL.csv"))
write.csv(AUC_Female_CU, file.path(tabDir, "est_Female_CU.csv"))
write.csv(AUC_Female, file.path(tabDir, "est_Female.csv"))
write.csv(n_Female, file.path(tabDir, "n_Female.csv"))

write.csv(AUC_Male_CL, file.path(tabDir, "est_Male_CL.csv"))
write.csv(AUC_Male_CU, file.path(tabDir, "est_Male_CU.csv"))
write.csv(AUC_Male, file.path(tabDir, "est_Male.csv"))

