rm(list=ls(all=TRUE))
here::i_am("README.md")
repoDir <- here::here()
datDir <- file.path(repoDir,"data")
codeDir <- file.path(repoDir, "code")
outputDir <- file.path(repoDir, "output")
figDir <- file.path(repoDir, "figures")
tabDir <- file.path(repoDir, "tables")

library(tidyverse)
library(ggplot2)
library(plyr)
source(file.path(codeDir, "common.R"))
source(file.path(codeDir, "summary/summaryUtils.R"))
options(scipen = 0)
sex = "Female"


#formatting functions
format.est <- function(x){format(round(x,2), nsmall = 2, digits = 2)}
format.ratio <- function(x){format(round(x,2), nsmall = 2, digits = 2)}
format.est2 <- function(x){round(x,0)}
format.p <- function(x){ifelse(x<0.001, "<0.001", ifelse(x == 1, "1", format(x, nsmall = 3, digits = 3)))}
format.antibodylevel <- function(x){
  if(x < 100){round(x,1)}
  else{as.numeric(format(round(x, 0), scienfitic = FALSE))}
  return(x)
}

library(scales)
scientific_10 <- function(x) {
  parse(text=gsub("1e\\+*", " 10^", scales::scientific_format()(x)))
}

AUC_p_value<- read.csv(file.path(outputDir, paste0("AUC_p_value_", sex, ".csv")))
D43_p_value<- read.csv(file.path(outputDir, paste0("D43_p_value_", sex, ".csv")))
D202_p_value<- read.csv(file.path(outputDir, paste0("D202_p_value_", sex, ".csv")))
rate_p_value <- read.csv(file.path(outputDir, paste0("rate_p_value_", sex, ".csv")))

AUC_CL <- read.csv(file.path(outputDir, paste0("AUC_", sex, "_CL.csv")))
AUC_CU <- read.csv(file.path(outputDir, paste0("AUC_", sex, "_CU.csv")))
AUC <- read.csv(file.path(outputDir, paste0("AUC_", sex, ".csv")))
n<- read.csv(file.path(outputDir, paste0("n_", sex, ".csv")))

#FWER p-value adjustment 
#excluding stage 1 nonnaive placebo vs. stage 2 nonnaive placebo, stage 2 nonnaive group A vs. stage 2 nonnaive group B
AUC_fwer_p_value <- AUC_p_value[, 1:14]
D43_fwer_p_value <- D43_p_value[, 1:14]
rate_fwer_p_value <- rate_p_value[, 1:14]
for(i in 1:13){
  #contrasts_adjust <- seq(1, 9, 1)[-c(2, 4, 9)]
  contrasts_adjust <- seq(1, 11, 1)[-c(2, 4, 8, 11)]
  AUC_fwer_p_value[contrasts_adjust, i+1] <- p.adjust(AUC_fwer_p_value[contrasts_adjust, i+1], method = "holm")
  D43_fwer_p_value[contrasts_adjust, i+1] <- p.adjust(D43_fwer_p_value[contrasts_adjust, i+1], method = "holm")
  rate_fwer_p_value[contrasts_adjust, i+1] <- p.adjust(rate_fwer_p_value[contrasts_adjust, i+1], method = "holm")
  
}



markers <- c("bindSpike", "bindSpike_beta", "bindSpike_alpha", "bindSpike_gamma", "bindSpike_delta1", "bindSpike_delta2", "bindSpike_delta3", "bindSpike_omicron",
             "pseudoneutid50", "pseudoneutid50_B.1.351", "pseudoneutid50_BA.1", "pseudoneutid50_BA.2", "pseudoneutid50_BA.4.5")

output <- tibble("Marker" = character(), "Comparison" = character(), "nA" = numeric(), "nB" = numeric(),
                 "AUC_diff" = character(), "Pvalue" = character(), "FWER_Pvalue" = character(),
                 "AUC_diff_est" = numeric(), "AUC_diff_CL" = numeric(), "AUC_diff_CU" = numeric(),
                 "AUC_GA" = character(), "AUC_GB" = character(),  
                 "AUC_GA_est" = numeric(), "AUC_GA_CL" = numeric(), "AUC_GA_CU" = numeric(),
                 "AUC_GB_est" = numeric(), "AUC_GB_CL" = numeric(), "AUC_GB_CU" = numeric(),
                 
                 "D43_diff" = character(), "D43_Pvalue" = character(), "D43_FWER_Pvalue" = character(),
                 "D43_diff_est" = numeric(), "D43_diff_CL" = numeric(), "D43_diff_CU" = numeric(),
                 "D43_GA" = character(), "D43_GB" = character(),
                 "D43_GA_est" = numeric(), "D43_GA_CL" = numeric(), "D43_GA_CU" = numeric(),
                 "D43_GB_est" = numeric(), "D43_GB_CL" = numeric(), "D43_GB_CU" = numeric(),
                 
                 "D202_diff" = character(),  "D202_Pvalue" = character(),
                 "D202_diff_est" = numeric(), "D202_diff_CL" = numeric(), "D202_diff_CU" = numeric(),
                 "D202_GA" = character(), "D202_GB" = character(),
                 "D202_GA_est" = numeric(), "D202_GA_CL" = numeric(), "D202_GA_CU" = numeric(),
                 "D202_GB_est" = numeric(), "D202_GB_CL" = numeric(), "D202_GB_CU" = numeric(),
                 
                 "rate_diff" = character(),  "rate_Pvalue" = character(), "rate_FWER_Pvalue" = character(),
                 "rate_diff_est" = numeric(), "rate_diff_CL" = numeric(), "rate_diff_CU" = numeric(),
                 "rate_GA" = character(), "rate_GB" = character(),
                 "rate_GA_est" = numeric(), "rate_GA_CL" = numeric(), "rate_GA_CU" = numeric(),
                 "rate_GB_est" = numeric(), "rate_GB_CL" = numeric(), "rate_GB_CU" = numeric()
)

for(mark in markers){
  for(c in 1:dim(contrasts)[1]){
    groupA <- filter(groups, Group == as.numeric(contrasts[c, 1]))$Label
    groupB <- filter(groups, Group == as.numeric(contrasts[c, 2]))$Label
    
    #AUC analysis
    est <- as.numeric(AUC[c, paste(mark, "diff", sep = "_")])
    cl <- as.numeric(AUC_CL[c, paste(mark, "diff", sep = "_")])
    cu <- as.numeric(AUC_CU[c, paste(mark, "diff", sep = "_")])
    p <- format.p(as.numeric(AUC_p_value[c, mark]))
    fwer_p <- format.p(as.numeric(AUC_fwer_p_value[c, mark]))
    nA_c<- n[c, paste(mark, "A", sep="_")]
    nB_c<- n[c, paste(mark, "B", sep="_")]
    diff <- paste0(format.ratio(10^est), " (", format.ratio(10^cl), ", ", format.ratio(10^cu), ")")
    
    est_GA <- as.numeric(AUC[c, paste(mark, "GA", sep = "_")])
    cl_GA <- as.numeric(AUC_CL[c, paste(mark, "GA", sep = "_")])
    cu_GA <- as.numeric(AUC_CU[c, paste(mark, "GA", sep = "_")])
    AUC_GA <- paste0(format.est2(10^est_GA), " (", format.est2(10^cl_GA), ", ", format.est2(10^cu_GA), ")")
    est_GB <- as.numeric(AUC[c, paste(mark, "GB", sep = "_")])
    cl_GB <- as.numeric(AUC_CL[c, paste(mark, "GB", sep = "_")])
    cu_GB <- as.numeric(AUC_CU[c, paste(mark, "GB", sep = "_")])
    AUC_GB <- paste0(format.est2(10^est_GB), " (", format.est2(10^cl_GB), ", ", format.est2(10^cu_GB), ")")
    
    #D43 analysis
    
    D43_est <- 10^as.numeric(AUC[c, paste(mark, "D43_diff", sep = "_")])
    D43_cl <- 10^as.numeric(AUC_CL[c, paste(mark, "D43_diff", sep = "_")])
    D43_cu <- 10^as.numeric(AUC_CU[c, paste(mark, "D43_diff", sep = "_")])
    D43_p <- format.p(as.numeric(D43_p_value[c, mark]))
    D43_fwer_p <- format.p(as.numeric(D43_fwer_p_value[c, mark]))
    D43_est_GA <- 10^as.numeric(AUC[c, paste(mark, "D43_GA", sep = "_")])
    D43_cl_GA <- 10^as.numeric(AUC_CL[c, paste(mark, "D43_GA", sep = "_")])
    D43_cu_GA <- 10^as.numeric(AUC_CU[c, paste(mark, "D43_GA", sep = "_")])
    D43_est_GB <- 10^as.numeric(AUC[c, paste(mark, "D43_GB", sep = "_")])
    D43_cl_GB <- 10^as.numeric(AUC_CL[c, paste(mark, "D43_GB", sep = "_")])
    D43_cu_GB <- 10^as.numeric(AUC_CU[c, paste(mark, "D43_GB", sep = "_")])
    
    D43_diff <- paste0(format.ratio(D43_est), " (", format.ratio(D43_cl), ", ", format.ratio(D43_cu), ")")
    D43_GA <- paste0(format.est2(D43_est_GA), " (", format.est2(D43_cl_GA), ", ", format.est2(D43_cu_GA), ")")
    D43_GB <- paste0(format.est2(D43_est_GB), " (", format.est2(D43_cl_GB), ", ", format.est2(D43_cu_GB), ")")
    
    #D202 analysis
    
    D202_est <- 10^as.numeric(AUC[c, paste(mark, "D202_diff", sep = "_")])
    D202_cl <- 10^as.numeric(AUC_CL[c, paste(mark, "D202_diff", sep = "_")])
    D202_cu <- 10^as.numeric(AUC_CU[c, paste(mark, "D202_diff", sep = "_")])
    D202_p <- format.p(as.numeric(D202_p_value[c, mark]))
    
    
    D202_est_GA <- 10^as.numeric(AUC[c, paste(mark, "D202_GA", sep = "_")])
    D202_cl_GA <- 10^as.numeric(AUC_CL[c, paste(mark, "D202_GA", sep = "_")])
    D202_cu_GA <- 10^as.numeric(AUC_CU[c, paste(mark, "D202_GA", sep = "_")])
    
    
    D202_est_GB <- 10^as.numeric(AUC[c, paste(mark, "D202_GB", sep = "_")])
    D202_cl_GB <- 10^as.numeric(AUC_CL[c, paste(mark, "D202_GB", sep = "_")])
    D202_cu_GB <- 10^as.numeric(AUC_CU[c, paste(mark, "D202_GB", sep = "_")])
    
    D202_diff <- paste0(format.est(D202_est), " (", format.est(D202_cl), ", ", format.est(D202_cu), ")")
    D202_GA <- paste0(format.est2(D202_est_GA), " (", format.est2(D202_cl_GA), ", ", format.est2(D202_cu_GA), ")")
    D202_GB <- paste0(format.est2(D202_est_GB), " (", format.est2(D202_cl_GB), ", ", format.est2(D202_cu_GB), ")")
    
    #rate analysis
    
    rate_est <- 10^as.numeric(AUC[c, paste(mark, "rate_diff", sep = "_")])
    rate_cl <- 10^as.numeric(AUC_CL[c, paste(mark, "rate_diff", sep = "_")])
    rate_cu <- 10^as.numeric(AUC_CU[c, paste(mark, "rate_diff", sep = "_")])
    rate_p <- format.p(as.numeric(rate_p_value[c, mark]))
    rate_fwer_p <- format.p(as.numeric(rate_fwer_p_value[c, mark]))
    
    rate_est_GA <- 10^as.numeric(AUC[c, paste(mark, "rate_GA", sep = "_")])
    rate_cl_GA <- 10^as.numeric(AUC_CL[c, paste(mark, "rate_GA", sep = "_")])
    rate_cu_GA <- 10^as.numeric(AUC_CU[c, paste(mark, "rate_GA", sep = "_")])
    
    
    rate_est_GB <- 10^as.numeric(AUC[c, paste(mark, "rate_GB", sep = "_")])
    rate_cl_GB <- 10^as.numeric(AUC_CL[c, paste(mark, "rate_GB", sep = "_")])
    rate_cu_GB <- 10^as.numeric(AUC_CU[c, paste(mark, "rate_GB", sep = "_")])
    
    rate_diff <- paste0(format.ratio(rate_est), " (", format.ratio(rate_cl), ", ", format.ratio(rate_cu), ")")
    rate_GB <- paste0(format.est(rate_est_GB), " (", format.est(rate_cl_GB), ", ", format.est(rate_cu_GB), ")")
    rate_GA <- paste0(format.est(rate_est_GA), " (", format.est(rate_cl_GA), ", ", format.est(rate_cu_GA), ")")
    
    
    
    output<- add_row(.data = output, "Comparison" = paste0(groupA, " vs. ", groupB), "Marker" = labf(mark), 
                     "nA" = nA_c, "nB" = nB_c,
                     "AUC_diff" = diff, "Pvalue" = p, "FWER_Pvalue" = fwer_p,
                     "AUC_GA" = AUC_GA, 
                     "AUC_GB" = AUC_GB, 
                     
                     "AUC_diff_est" = est, "AUC_diff_CL" = cl, "AUC_diff_CU" = cu,
                     "AUC_GA_est" = est_GA, "AUC_GA_CL" = cl_GA, "AUC_GA_CU" = cu_GA,
                     "AUC_GB_est" = est_GB, "AUC_GB_CL" = cl_GB, "AUC_GB_CU" = cu_GB,
                     
                     "D43_diff" = D43_diff, "D43_Pvalue" = D43_p, "D43_FWER_Pvalue" = D43_fwer_p,
                     "D43_diff_est" = log10(D43_est), "D43_diff_CL" = log10(D43_cl), "D43_diff_CU" = log10(D43_cu), 
                     "D43_GA" = D43_GA, 
                     "D43_GB" = D43_GB, 
                     "D43_GA_est" = D43_est_GA, "D43_GA_CL" = D43_cl_GA, "D43_GA_CU" = D43_cu_GA,
                     "D43_GB_est" = D43_est_GB, "D43_GB_CL" = D43_cl_GB, "D43_GB_CU" = D43_cu_GB,
                     
                     
                     "D202_diff" = D202_diff, "D202_Pvalue" = D202_p,
                     "D202_diff_est" = log10(D202_est), "D202_diff_CL" = log10(D202_cl), "D202_diff_CU" = log10(D202_cu), 
                     "D202_GA" = D202_GA, 
                     "D202_GB" = D202_GB, 
                     "D202_GA_est" = D202_est_GA, "D202_GA_CL" = D202_cl_GA, "D202_GA_CU" = D202_cu_GA,
                     "D202_GB_est" = D202_est_GB, "D202_GB_CL" = D202_cl_GB, "D202_GB_CU" = D202_cu_GB,
                     
                     "rate_diff" = rate_diff, "rate_Pvalue" = rate_p, "rate_FWER_Pvalue" = rate_fwer_p,
                     "rate_diff_est" = log10(rate_est), "rate_diff_CL" = log10(rate_cl), "rate_diff_CU" = log10(rate_cu), 
                     "rate_GA" = rate_GA, 
                     "rate_GB" = rate_GB, 
                     "rate_GA_est" = rate_est_GA, "rate_GA_CL" = rate_cl_GA, "rate_GA_CU" = rate_cu_GA,
                     "rate_GB_est" = rate_est_GB, "rate_GB_CL" = rate_cl_GB, "rate_GB_CU" = rate_cu_GB
    )
    
  }
}

output$groupA <- laply(output$Comparison, function(x)strsplit(x, split = " vs. ")[[1]][1])
output$groupB <- laply(output$Comparison, function(x)strsplit(x, split = " vs. ")[[1]][2])

#Figures 2, 3
comparisons <- c("1-V-NN vs. 1-P-NN", "2-V-NN vs. 2-P-NN")
plotSummary_D43(comparisons, labels.x = c("Non-naïve MV", "Stage 1 Non-\nnaïve Placebo", "Non-naïve BV", "Stage 2 Non-\nnaïve Placebo"), figureLabel = "VvsP")
plotSummary_AUC(comparisons, labels.x = c("Non-naïve MV", "Stage 1 Non-\nnaïve Placebo", "Non-naïve BV", "Stage 2 Non-\nnaïve Placebo"), figureLabel = "VvsP")
plotSummary_GMR(comparisons, labels.x = c("Non-naïve MV", "Stage 1 Non-\nnaïve Placebo", "Non-naïve BV", "Stage 2 Non-\nnaïve Placebo"), figureLabel = "VvsP")

#Figures 4, 5
comparisons <- c("2-V-NN vs. 1-V-NN", "2-V-N vs. 1-V-N")
plotSummary_D43(comparisons, labels.x = c("Non-naïve MV","Non-naïve BV",  "Naïve MV", "Naïve BV"), figureLabel = "BVvsMV",
                IgG_ylim = c(3.5, 6.5), ID50_ylim = c(0.8, 5.1), reverseGroup = TRUE)
plotSummary_AUC(comparisons, labels.x = c("Non-naïve MV","Non-naïve BV",  "Naïve MV", "Naïve BV"), figureLabel = "BVvsMV",
                IgG_ylim = c(3.5, 6), ID50_ylim = c(0.8, 4.6), reverseGroup = TRUE)
plotSummary_GMR(comparisons, labels.x = c("Non-naïve MV","Non-naïve BV",  "Naïve MV", "Naïve BV"), figureLabel = "BVvsMV",
                IgG_ylim = c(0, 0.6), ID50_ylim = c(-0.05, 1.03), reverseGroup = TRUE)



#Table 1
permuTests_bootCI_j <- filter(output, Comparison %in% c("1-V-NN vs. 1-V-N", "2-V-NN vs. 2-V-N") 
                              & Marker %in% c("bAb-IgG Spike Omicron", "bAb-IgG Spike Index", "nAb-ID50 Omicron-BA.4/BA.5",
                                              "nAb-ID50 Reference"))

tab <- tibble("Measure" = character(), "Mark" = character(), "1-V-NN vs. 1-V-N" = character(), "P1" = character(), "FWER P1" = character(),
              "2-V-NN vs. 2-V-N" = character(), "P2" = character(), "FWER P2" = character())
for(markeri in c("nAb-ID50 Reference", "nAb-ID50 Omicron-BA.4/BA.5", "bAb-IgG Spike Index", "bAb-IgG Spike Omicron")){
  tmp1 <- filter(permuTests_bootCI_j, Marker == markeri & Comparison == "1-V-NN vs. 1-V-N")
  tmp2 <- filter(permuTests_bootCI_j, Marker == markeri & Comparison == "2-V-NN vs. 2-V-N")
  D43_diff_1 <- tmp1$D43_diff
  D43_diff_2 <- tmp2$D43_diff
  
  
  tab <- add_row(.data = tab,"Measure" = ifelse(markeri == "nAb-ID50 Reference", "Day 43", ""), "Mark" = markeri, 
                 "1-V-NN vs. 1-V-N" = D43_diff_1, "P1" = tmp1$D43_Pvalue, "FWER P1" = tmp1$D43_FWER_Pvalue,
                 "2-V-NN vs. 2-V-N" = D43_diff_2, "P2" = tmp2$D43_Pvalue, "FWER P2" = tmp2$D43_FWER_Pvalue)
}
for(markeri in c("nAb-ID50 Reference", "nAb-ID50 Omicron-BA.4/BA.5", "bAb-IgG Spike Index", "bAb-IgG Spike Omicron")){
  tmp1 <- filter(permuTests_bootCI_j, Marker == markeri & Comparison == "1-V-NN vs. 1-V-N")
  tmp2 <- filter(permuTests_bootCI_j, Marker == markeri & Comparison == "2-V-NN vs. 2-V-N")
  AUC_diff_1 <- tmp1$AUC_diff
  AUC_diff_2 <- tmp2$AUC_diff
  
  tab <- add_row(.data = tab,"Measure" = ifelse(markeri == "nAb-ID50 Reference", "AUC*", ""), "Mark" = markeri, 
                 "1-V-NN vs. 1-V-N" = AUC_diff_1, "P1" = tmp1$Pvalue, "FWER P1" = tmp1$FWER_Pvalue,
                 "2-V-NN vs. 2-V-N" = AUC_diff_2, "P2" = tmp2$Pvalue, "FWER P2" = tmp2$FWER_Pvalue)
}

for(markeri in c("nAb-ID50 Reference", "nAb-ID50 Omicron-BA.4/BA.5", "bAb-IgG Spike Index", "bAb-IgG Spike Omicron")){
  tmp1 <- filter(permuTests_bootCI_j, Marker == markeri & Comparison == "1-V-NN vs. 1-V-N")
  tmp2 <- filter(permuTests_bootCI_j, Marker == markeri & Comparison == "2-V-NN vs. 2-V-N")
  rate_diff_1 <- tmp1$rate_diff
  rate_diff_2 <- tmp2$rate_diff
  
  tab <- add_row(.data = tab,"Measure" = ifelse(markeri == "nAb-ID50 Reference", "rate", ""), "Mark" = markeri, 
                 "1-V-NN vs. 1-V-N" = rate_diff_1, "P1" = tmp1$rate_Pvalue, "FWER P1" = tmp1$rate_FWER_Pvalue,
                 "2-V-NN vs. 2-V-N" = rate_diff_2, "P2" = tmp2$rate_Pvalue, "FWER P2" = tmp2$rate_FWER_Pvalue)
}

write.csv(tab, file.path(tabDir,"naiveVsnonnaive_mainMarkers.csv"), row.names = FALSE)
