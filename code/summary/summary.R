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

#functions
format.est <- function(x){format(round(x,2), nsmall = 1, digits = 1)}
format.est2 <- function(x){round(x,0)}
format.p <- function(x){ifelse(x<0.001, "<0.001", ifelse(x == 1, "1", format(x, nsmall = 1, digits = 3)))}
format.antibodylevel <- function(x){
  if(x < 100){round(x,1)}
  else{as.numeric(format(round(x, 0), scienfitic = FALSE))}
  return(x)
}

library(scales)
scientific_10 <- function(x) {
  y <- NULL
  for(i in 1:length(x)){
    y[i] <- parse(text=gsub("1e\\+*", " 10^", scales::scientific_format()(x[i])))
  }
  return(y)
}

durability_p_value<- read.csv(file.path(tabDir, paste0("durability_p_value_", sex, ".csv")))
D43_p_value<- read.csv(file.path(tabDir, paste0("D43_p_value_", sex, ".csv")))
D202toD43GMR_p_value <- read.csv(file.path(tabDir, paste0("D202toD43GMR_p_value_", sex, ".csv")))

est_CL <- read.csv(file.path(tabDir, paste0("est_", sex, "_CL.csv")))
est_CU <- read.csv(file.path(tabDir, paste0("est_", sex, "_CU.csv")))
est <- read.csv(file.path(tabDir, paste0("est_", sex, ".csv")))
n<- read.csv(file.path(tabDir, paste0("n_", sex, ".csv")))

#FWER p-value adjustment 
#excluding stage 1 nonnaive placebo vs. stage 2 nonnaive placebo, stage 2 nonnaive group A vs. stage 2 nonnaive group B
durability_fwer_p_value <- durability_p_value[, 1:14]
D43_fwer_p_value <- D43_p_value[, 1:14]
D202toD43GMR_FWER_p_value <- D202toD43GMR_p_value[, 1:14]
for(i in 1:13){
  contrasts_adjust <- seq(1, 11, 1)[-c(2, 4, 8, 11)]
  durability_fwer_p_value[contrasts_adjust, i+1] <- p.adjust(durability_fwer_p_value[contrasts_adjust, i+1], method = "holm")
  D43_fwer_p_value[contrasts_adjust, i+1] <- p.adjust(D43_fwer_p_value[contrasts_adjust, i+1], method = "holm")
  D202toD43GMR_FWER_p_value[contrasts_adjust, i+1] <- p.adjust(D202toD43GMR_FWER_p_value[contrasts_adjust, i+1], method = "holm")
  
}



markers <- c("bindSpike", "bindSpike_beta", "bindSpike_alpha", "bindSpike_gamma", "bindSpike_delta1", "bindSpike_delta2", "bindSpike_delta3", "bindSpike_omicron",
             "pseudoneutid50", "pseudoneutid50_B.1.351", "pseudoneutid50_BA.1", "pseudoneutid50_BA.2", "pseudoneutid50_BA.4.5")

output <- tibble("Marker" = character(), "Comparison" = character(), "nA" = numeric(), "nB" = numeric(),
                 "durability_diff" = character(), "Pvalue" = character(), "FWER_Pvalue" = character(),
                 "durability_diff_est" = numeric(), "durability_diff_CL" = numeric(), "durability_diff_CU" = numeric(),
                 "durability_GA" = character(), "durability_GB" = character(),  
                 "durability_GA_est" = numeric(), "durability_GA_CL" = numeric(), "durability_GA_CU" = numeric(),
                 "durability_GB_est" = numeric(), "durability_GB_CL" = numeric(), "durability_GB_CU" = numeric(),
  
                 "D43_diff" = character(), "D43_Pvalue" = character(), "D43_FWER_Pvalue" = character(),
                 "D43_diff_est" = numeric(), "D43_diff_CL" = numeric(), "D43_diff_CU" = numeric(),
                 "D43_GA" = character(), "D43_GB" = character(),
                 "D43_GA_est" = numeric(), "D43_GA_CL" = numeric(), "D43_GA_CU" = numeric(),
                 "D43_GB_est" = numeric(), "D43_GB_CL" = numeric(), "D43_GB_CU" = numeric(),
               
                 "D202toD43GMR_diff" = character(),  "D202toD43GMR_pvalue" = character(), "D202toD43GMR_FWER_Pvalue" = character(),
                 "D202toD43GMR_diff_est" = numeric(), "D202toD43GMR_diff_CL" = numeric(), "D202toD43GMR_diff_CU" = numeric(),
                 "D202toD43GMR_GA" = character(), "D202toD43GMR_GB" = character(),
                 "D202toD43GMR_GA_est" = numeric(), "D202toD43GMR_GA_CL" = numeric(), "D202toD43GMR_GA_CU" = numeric(),
                 "D202toD43GMR_GB_est" = numeric(), "D202toD43GMR_GB_CL" = numeric(), "D202toD43GMR_GB_CU" = numeric()
                 )

for(mark in markers){
  for(c in 1:dim(contrasts)[1]){
    groupA <- filter(groups, Group == as.numeric(contrasts[c, 1]))$Label
    groupB <- filter(groups, Group == as.numeric(contrasts[c, 2]))$Label
  
    #durability results
    durability_est <- as.numeric(est[c, paste(mark, "diff", sep = "_")])
    durability_cl <- as.numeric(est_CL[c, paste(mark, "diff", sep = "_")])
    durability_cu <- as.numeric(est_CU[c, paste(mark, "diff", sep = "_")])
    durability_p <- format.p(as.numeric(durability_p_value[c, mark]))
    durability_fwer_p <- format.p(as.numeric(durability_fwer_p_value[c, mark]))
    nA_c<- n[c, paste(mark, "A", sep="_")]
    nB_c<- n[c, paste(mark, "B", sep="_")]
    durability_diff <- paste0(format.est(10^durability_est), " (", format.est(10^durability_cl), ", ", format.est(10^durability_cu), ")")
   
    durability_est_GA <- as.numeric(est[c, paste(mark, "GA", sep = "_")])
    durability_cl_GA <- as.numeric(est_CL[c, paste(mark, "GA", sep = "_")])
    durability_cu_GA <- as.numeric(est_CU[c, paste(mark, "GA", sep = "_")])
    durability_GA <- paste0(format.est2(10^durability_est_GA), " (", format.est2(10^durability_cl_GA), ", ", format.est2(10^durability_cu_GA), ")")
    durability_est_GB <- as.numeric(est[c, paste(mark, "GB", sep = "_")])
    durability_cl_GB <- as.numeric(est_CL[c, paste(mark, "GB", sep = "_")])
    durability_cu_GB <- as.numeric(est_CU[c, paste(mark, "GB", sep = "_")])
    durability_GB <- paste0(format.est2(10^durability_est_GB), " (", format.est2(10^durability_cl_GB), ", ", format.est2(10^durability_cu_GB), ")")
    
    #D43 results
    D43_est <- 10^as.numeric(est[c, paste(mark, "D43_diff", sep = "_")])
    D43_cl <- 10^as.numeric(est_CL[c, paste(mark, "D43_diff", sep = "_")])
    D43_cu <- 10^as.numeric(est_CU[c, paste(mark, "D43_diff", sep = "_")])
    D43_p <- format.p(as.numeric(D43_p_value[c, mark]))
    D43_fwer_p <- format.p(as.numeric(D43_fwer_p_value[c, mark]))
    D43_est_GA <- 10^as.numeric(est[c, paste(mark, "D43_GA", sep = "_")])
    D43_cl_GA <- 10^as.numeric(est_CL[c, paste(mark, "D43_GA", sep = "_")])
    D43_cu_GA <- 10^as.numeric(est_CU[c, paste(mark, "D43_GA", sep = "_")])
    D43_est_GB <- 10^as.numeric(est[c, paste(mark, "D43_GB", sep = "_")])
    D43_cl_GB <- 10^as.numeric(est_CL[c, paste(mark, "D43_GB", sep = "_")])
    D43_cu_GB <- 10^as.numeric(est_CU[c, paste(mark, "D43_GB", sep = "_")])
    
    D43_diff <- paste0(format.est(D43_est), " (", format.est(D43_cl), ", ", format.est(D43_cu), ")")
    D43_GA <- paste0(format.est2(D43_est_GA), " (", format.est2(D43_cl_GA), ", ", format.est2(D43_cu_GA), ")")
    D43_GB <- paste0(format.est2(D43_est_GB), " (", format.est2(D43_cl_GB), ", ", format.est2(D43_cu_GB), ")")
    
    
    #D202 to D43 GMR results
    D202toD43GMR_est <- 10^as.numeric(est[c, paste(mark, "rate_diff", sep = "_")])
    D202toD43GMR_cl <- 10^as.numeric(est_CL[c, paste(mark, "rate_diff", sep = "_")])
    D202toD43GMR_cu <- 10^as.numeric(est_CU[c, paste(mark, "rate_diff", sep = "_")])
    D202toD43GMR_p <- format.p(as.numeric(D202toD43GMR_p_value[c, mark]))
    D202toD43GMR_FWER_p <- format.p(as.numeric(D202toD43GMR_FWER_p_value[c, mark]))
    
    D202toD43GMR_est_GA <- 10^as.numeric(est[c, paste(mark, "rate_GA", sep = "_")])
    D202toD43GMR_cl_GA <- 10^as.numeric(est_CL[c, paste(mark, "rate_GA", sep = "_")])
    D202toD43GMR_cu_GA <- 10^as.numeric(est_CU[c, paste(mark, "rate_GA", sep = "_")])
    
    
    D202toD43GMR_est_GB <- 10^as.numeric(est[c, paste(mark, "rate_GB", sep = "_")])
    D202toD43GMR_cl_GB <- 10^as.numeric(est_CL[c, paste(mark, "rate_GB", sep = "_")])
    D202toD43GMR_cu_GB <- 10^as.numeric(est_CU[c, paste(mark, "rate_GB", sep = "_")])
    
    D202toD43GMR_diff <- paste0(format.est(D202toD43GMR_est), " (", format.est(D202toD43GMR_cl), ", ", format.est(D202toD43GMR_cu), ")")
    D202toD43GMR_GB <- paste0(format.est(D202toD43GMR_est_GB), " (", format.est(D202toD43GMR_cl_GB), ", ", format.est(D202toD43GMR_cu_GB), ")")
    D202toD43GMR_GA <- paste0(format.est(D202toD43GMR_est_GA), " (", format.est(D202toD43GMR_cl_GA), ", ", format.est(D202toD43GMR_cu_GA), ")")
    
    
    
    output<- add_row(.data = output, "Comparison" = paste0(groupA, " vs. ", groupB), "Marker" = labf(mark), 
                     "nA" = nA_c, "nB" = nB_c,
                     "durability_diff" = durability_diff, "Pvalue" = durability_p, "FWER_Pvalue" = durability_fwer_p,
                     "durability_GA" = durability_GA, 
                     "durability_GB" = durability_GB, 
                     
                     "durability_diff_est" = durability_est, "durability_diff_CL" = durability_cl, "durability_diff_CU" = durability_cu,
                     "durability_GA_est" = durability_est_GA, "durability_GA_CL" = durability_cl_GA, "durability_GA_CU" = durability_cu_GA,
                     "durability_GB_est" = durability_est_GB, "durability_GB_CL" = durability_cl_GB, "durability_GB_CU" = durability_cu_GB,
                   
                     "D43_diff" = D43_diff, "D43_Pvalue" = D43_p, "D43_FWER_Pvalue" = D43_fwer_p,
                     "D43_diff_est" = log10(D43_est), "D43_diff_CL" = log10(D43_cl), "D43_diff_CU" = log10(D43_cu), 
                     "D43_GA" = D43_GA, 
                     "D43_GB" = D43_GB, 
                     "D43_GA_est" = D43_est_GA, "D43_GA_CL" = D43_cl_GA, "D43_GA_CU" = D43_cu_GA,
                     "D43_GB_est" = D43_est_GB, "D43_GB_CL" = D43_cl_GB, "D43_GB_CU" = D43_cu_GB,
                     
                     "D202toD43GMR_diff" = D202toD43GMR_diff, "D202toD43GMR_pvalue" = D202toD43GMR_p, "D202toD43GMR_FWER_Pvalue" = D202toD43GMR_FWER_p,
                     "D202toD43GMR_diff_est" = log10(D202toD43GMR_est), "D202toD43GMR_diff_CL" = log10(D202toD43GMR_cl), "D202toD43GMR_diff_CU" = log10(D202toD43GMR_cu), 
                     "D202toD43GMR_GA" = D202toD43GMR_GA, 
                     "D202toD43GMR_GB" = D202toD43GMR_GB, 
                     "D202toD43GMR_GA_est" = D202toD43GMR_est_GA, "D202toD43GMR_GA_CL" = D202toD43GMR_cl_GA, "D202toD43GMR_GA_CU" = D202toD43GMR_cu_GA,
                     "D202toD43GMR_GB_est" = D202toD43GMR_est_GB, "D202toD43GMR_GB_CL" = D202toD43GMR_cl_GB, "D202toD43GMR_GB_CU" = D202toD43GMR_cu_GB
                    )
   
  }
}

output$groupA <- laply(output$Comparison, function(x)strsplit(x, split = " vs. ")[[1]][1])
output$groupB <- laply(output$Comparison, function(x)strsplit(x, split = " vs. ")[[1]][2])
#output$Marker[output$Marker == "bAb-IgG Spike Reference"] = "bAb-IgG Spike Index"

#plot the summaries for the main text 
comparisons <- c("1-V-NN vs. 1-P-NN", "2-V-NN vs. 2-P-NN")
plotSummary_D43(comparisons, labels.x = c("Non-naïve MV", "Stage 1 Non-\nnaïve Placebo", "Non-naïve BV", "Stage 2 Non-\nnaïve Placebo"), 
                figureLabel = "_stage1and2_Vaccine_vs_Placebo_")
plotSummary_durability(comparisons, labels.x = c("Non-naïve MV", "Stage 1 Non-\nnaïve Placebo", "Non-naïve BV", "Stage 2 Non-\nnaïve Placebo"), 
                       figureLabel = "_stage1and2_Vaccine_vs_Placebo_")
plotSummary_GMR(comparisons, labels.x = c("Non-naïve MV", "Stage 1 Non-\nnaïve Placebo", "Non-naïve BV", "Stage 2 Non-\nnaïve Placebo"), 
                figureLabel = "_stage1and2_Vaccine_vs_Placebo_")

comparisons <- c("2-V-NN vs. 1-V-NN","2-V-N vs. 1-V-N")
plotSummary_D43(comparisons, labels.x = c("Non-naïve MV","Non-naïve BV",  "Naïve MV", "Naïve BV"), 
                figureLabel = "stage2_vs_stage1",
                IgG_ylim = c(3.5, 6.5), ID50_ylim = c(0.8, 5.1), reverseGroup = TRUE)
plotSummary_durability(comparisons, labels.x = c("Non-naïve MV","Non-naïve BV",  "Naïve MV", "Naïve BV"), 
                       figureLabel = "stage2_vs_stage1",
                IgG_ylim = c(3.5, 6), ID50_ylim = c(0.8, 4.6), reverseGroup = TRUE)
plotSummary_GMR(comparisons, labels.x = c("Non-naïve MV","Non-naïve BV",  "Naïve MV", "Naïve BV"), 
                figureLabel = "stage2_vs_stage1",
                IgG_ylim = c(0, 0.6), ID50_ylim = c(-0.05, 1.03), reverseGroup = TRUE)


##################################################################################################################################
#table summaries for all markers 
##################################################################################################################################
#vaccine vs placebo
allMarkers<- c("nAb-ID50 Reference", "nAb-ID50 Beta", "nAb-ID50 Omicron-BA.1", "nAb-ID50 Omicron-BA.2","nAb-ID50 Omicron-BA.4/BA.5", "bAb-IgG Spike Index", "bAb-IgG Spike Alpha", "bAb-IgG Spike Beta",
                  "bAb-IgG Spike Delta-AY.4", "bAb-IgG Spike Delta-B.1.617.2", "bAb-IgG Spike Delta-B.1.617.2/AY.4", "bAb-IgG Spike Gamma", "bAb-IgG Spike Omicron")
permuTests_bootCI_j <- filter(output, Comparison== c("1-V-NN vs. 1-P-NN") & Marker %in% allMarkers)

tab_stage1 <- tibble("Stage" = character(),"Measure" = character(), "Mark" = character(), "nA" = numeric(), "Nonnaive Vaccine" = character(), "nB" = numeric(), "Nonnaive Placebo" = character(), 
                     "Ratio" = character(), "Unadjusted P" = character(), "FWER P" = character())
for(markeri in allMarkers){
  tmp <- filter(permuTests_bootCI_j, Marker == markeri)
  tab_stage1 <- add_row(.data = tab_stage1,"Stage" = "Stage 1","Measure" = "Day 43", "Mark" = markeri, "nA" = tmp$nA, "Nonnaive Vaccine" = tmp$D43_GA, "nB" = tmp$nB, "Nonnaive Placebo" = tmp$D43_GB,
                        "Ratio" = tmp$D43_diff, "Unadjusted P" = tmp$D43_Pvalue, "FWER P" = tmp$D43_FWER_Pvalue)
}
for(markeri in allMarkers){
  tmp <- filter(permuTests_bootCI_j, Marker == markeri)
  tab_stage1 <- add_row(.data = tab_stage1,"Stage" = "Stage 1","Measure" = "Durability", "Mark" = markeri, "nA" = tmp$nA, "Nonnaive Vaccine" = tmp$durability_GA, "nB" = tmp$nB, 
                        "Nonnaive Placebo" = tmp$durability_GB,
                        "Ratio" = tmp$durability_diff, "Unadjusted P" = tmp$Pvalue, "FWER P" = tmp$FWER_Pvalue)
}
for(markeri in allMarkers){
  tmp <- filter(permuTests_bootCI_j, Marker == markeri)
  tab_stage1 <- add_row(.data = tab_stage1,"Stage" = "Stage 1","Measure" = "D202 to D43 GMR", "Mark" = markeri, "nA" = tmp$nA, "Nonnaive Vaccine" = tmp$D202toD43GMR_GA, "nB" = tmp$nB, 
                        "Nonnaive Placebo" = tmp$D202toD43GMR_GB,
                        "Ratio" = tmp$D202toD43GMR_diff, "Unadjusted P" = tmp$D202toD43GMR_pvalue, "FWER P" = tmp$D202toD43GMR_FWER_Pvalue)
}

permuTests_bootCI_j <- filter(output, Comparison== c("2-V-NN vs. 2-P-NN") & Marker %in% allMarkers)

tab_stage2 <- tibble("Stage" = character(), "Measure" = character(), "Mark" = character(), "nA" = numeric(), "Nonnaive Vaccine" = character(), "nB" = numeric(), "Nonnaive Placebo" = character(), 
                     "Ratio" = character(), "Unadjusted P" = character(), "FWER P" = character())
for(markeri in allMarkers){
  tmp <- filter(permuTests_bootCI_j, Marker == markeri)
  tab_stage2 <- add_row(.data = tab_stage2,"Stage" = "Stage 2","Measure" = "Day 43", "Mark" = markeri, "nA" = tmp$nA, "Nonnaive Vaccine" = tmp$D43_GA, "nB" = tmp$nB, "Nonnaive Placebo" = tmp$D43_GB,
                        "Ratio" = tmp$D43_diff, "Unadjusted P" = tmp$D43_Pvalue, "FWER P" = tmp$D43_FWER_Pvalue)
}
for(markeri in allMarkers){
  tmp <- filter(permuTests_bootCI_j, Marker == markeri)
  tab_stage2 <- add_row(.data = tab_stage2,"Stage" = "Stage 2","Measure" = "Durability", "Mark" = markeri, "nA" = tmp$nA, "Nonnaive Vaccine" = tmp$durability_GA, "nB" = tmp$nB, 
                        "Nonnaive Placebo" = tmp$durability_GB,
                        "Ratio" = tmp$durability_diff, "Unadjusted P" = tmp$Pvalue, "FWER P" = tmp$FWER_Pvalue)
}
for(markeri in allMarkers){
  tmp <- filter(permuTests_bootCI_j, Marker == markeri)
  tab_stage2 <- add_row(.data = tab_stage2,"Stage" = "Stage 2","Measure" = "D202 to D43 GMR", "Mark" = markeri, "nA" = tmp$nA, "Nonnaive Vaccine" = tmp$D202toD43GMR_GA, "nB" = tmp$nB, 
                        "Nonnaive Placebo" = tmp$D202toD43GMR_GB,
                        "Ratio" = tmp$D202toD43GMR_diff, "Unadjusted P" = tmp$D202toD43GMR_pvalue, "FWER P" = tmp$D202toD43GMR_FWER_Pvalue)
}


write.csv(rbind(tab_stage1,tab_stage2) , file.path(tabDir,"stage1and2_nonnaive_vaccine_vs_placebo.csv"), row.names = FALSE)



##################################################################################################################################
permuTests_bootCI_j <- filter(output, Comparison== "2-V-N vs. 1-V-N" & Marker %in% allMarkers)

tab <- tibble("Measure" = character(), "Mark" = character(), "nA" = numeric(), "Stage 1 Naive" = character(), "nB" = numeric(), 
              "Stage 2 Naive" = character(), "Ratio" = character(), "Unadjusted P" = character(), "FWER P" = character())
for(markeri in allMarkers){
  tmp <- filter(permuTests_bootCI_j, Marker == markeri)
  tab <- add_row(.data = tab,"Measure" =  "Day 43", "Mark" = markeri, "nA" = tmp$nB, "Stage 1 Naive" = tmp$D43_GB, "nB" = tmp$nA, "Stage 2 Naive" = tmp$D43_GA,
                 "Ratio" = tmp$D43_diff, "Unadjusted P" = tmp$D43_Pvalue, "FWER P" = tmp$D43_FWER_Pvalue)
}
for(markeri in allMarkers){
  tmp <- filter(permuTests_bootCI_j, Marker == markeri)
  tab <- add_row(.data = tab,"Measure" = "Durability", "Mark" = markeri, "nA" = tmp$nB, "Stage 1 Naive" = tmp$durability_GB, "nB" = tmp$nA, 
                 "Stage 2 Naive" = tmp$durability_GA,
                 "Ratio" = tmp$durability_diff, "Unadjusted P" = tmp$Pvalue, "FWER P" = tmp$FWER_Pvalue)
}
for(markeri in allMarkers){
  tmp <- filter(permuTests_bootCI_j, Marker == markeri)
  tab <- add_row(.data = tab,"Measure" = "D202 to D43 GMR", "Mark" = markeri, "nA" = tmp$nB, "Stage 1 Naive" = tmp$D202toD43GMR_GB, "nB" = tmp$nA, 
                 "Stage 2 Naive" = tmp$D202toD43GMR_GA,
                 "Ratio" = tmp$D202toD43GMR_diff, "Unadjusted P" = tmp$D202toD43GMR_pvalue, "FWER P" = tmp$D202toD43GMR_FWER_Pvalue)
}
write.csv(tab, file.path(tabDir,"naive_stage2_vs_stage1.csv"), row.names = FALSE)


###############################################################################################################################################################################
permuTests_bootCI_j <- filter(output, Comparison %in% c("1-V-NN vs. 1-V-N", "2-V-NN vs. 2-V-N") 
                              & Marker %in% allMarkers)

tab <- tibble("Measure" = character(), "Mark" = character(), "1-V-NN vs. 1-V-N" = character(), "Unadjusted P1" = character(), "FWER P1" = character(),
              "2-V-NN vs. 2-V-N" = character(), "Unadjusted P2" = character(), "FWER P2" = character())
for(markeri in allMarkers){
  tmp1 <- filter(permuTests_bootCI_j, Marker == markeri & Comparison == "1-V-NN vs. 1-V-N")
  tmp2 <- filter(permuTests_bootCI_j, Marker == markeri & Comparison == "2-V-NN vs. 2-V-N")
  D43_diff_1 <- tmp1$D43_diff
  D43_diff_2 <- tmp2$D43_diff
  
  
  tab <- add_row(.data = tab,"Measure" = "Day 43", "Mark" = markeri, 
                 "1-V-NN vs. 1-V-N" = D43_diff_1, "Unadjusted P1" = tmp1$D43_Pvalue, "FWER P1" = tmp1$D43_FWER_Pvalue,
                 "2-V-NN vs. 2-V-N" = D43_diff_2, "Unadjusted P2" = tmp2$D43_Pvalue, "FWER P2" = tmp2$D43_FWER_Pvalue)
}
for(markeri in allMarkers){
  tmp1 <- filter(permuTests_bootCI_j, Marker == markeri & Comparison == "1-V-NN vs. 1-V-N")
  tmp2 <- filter(permuTests_bootCI_j, Marker == markeri & Comparison == "2-V-NN vs. 2-V-N")
  durability_diff_1 <- tmp1$durability_diff
  durability_diff_2 <- tmp2$durability_diff
  
  tab <- add_row(.data = tab,"Measure" = "Durability", "Mark" = markeri, 
                 "1-V-NN vs. 1-V-N" = durability_diff_1, "Unadjusted P1" = tmp1$Pvalue,  "FWER P1" = tmp1$FWER_Pvalue,
                 "2-V-NN vs. 2-V-N" = durability_diff_2, "Unadjusted P2" = tmp2$Pvalue, "FWER P2" = tmp2$FWER_Pvalue)
}
for(markeri in allMarkers){
  tmp1 <- filter(permuTests_bootCI_j, Marker == markeri & Comparison == "1-V-NN vs. 1-V-N")
  tmp2 <- filter(permuTests_bootCI_j, Marker == markeri & Comparison == "2-V-NN vs. 2-V-N")
  D202toD43GMR_diff_1 <- tmp1$D202toD43GMR_diff
  D202toD43GMR_diff_2 <- tmp2$D202toD43GMR_diff
  
  tab <- add_row(.data = tab,"Measure" = "D202 to D43 GMR", "Mark" = markeri, 
                 "1-V-NN vs. 1-V-N" = D202toD43GMR_diff_1, "Unadjusted P1" = tmp1$D202toD43GMR_pvalue,  "FWER P1" = tmp1$D202toD43GMR_FWER_Pvalue,
                 "2-V-NN vs. 2-V-N" = D202toD43GMR_diff_2, "Unadjusted P2" = tmp2$D202toD43GMR_pvalue, "FWER P2" = tmp2$D202toD43GMR_FWER_Pvalue,)
}
write.csv(tab, file.path(tabDir,"stage1and2_nonnaive_vs_naive.csv"), row.names = FALSE)


#compare stage 1 vs stage 2 nonnaive
permuTests_bootCI_j <- filter(output, Comparison %in% c("2-V-NN vs. 1-V-NN", "2-P-NN vs. 1-P-NN") 
                              & Marker %in% allMarkers)

tab <- tibble("Measure" = character(), "Mark" = character(), "2-V-NN vs. 1-V-NN" = character(), "Unadjusted P1" = character(), "FWER P1" = character(),
              "2-P-NN vs. 1-P-NN" = character(), "Unadjusted P2" = character(), "FWER P2" = character())
for(markeri in allMarkers){
  tmp1 <- filter(permuTests_bootCI_j, Marker == markeri & Comparison == "2-V-NN vs. 1-V-NN")
  tmp2 <- filter(permuTests_bootCI_j, Marker == markeri & Comparison == "2-P-NN vs. 1-P-NN")
  
  tab <- add_row(.data = tab,"Measure" = "Day 43", "Mark" = markeri, 
                 "2-V-NN vs. 1-V-NN" = tmp1$D43_diff, "Unadjusted P1" = tmp1$D43_Pvalue, "FWER P1" = tmp1$D43_FWER_Pvalue,
                 "2-P-NN vs. 1-P-NN" = tmp2$D43_diff, "Unadjusted P2" = tmp2$D43_Pvalue, "FWER P2" = tmp2$D43_FWER_Pvalue)
}
for(markeri in allMarkers){
  tmp1 <- filter(permuTests_bootCI_j, Marker == markeri & Comparison == "2-V-NN vs. 1-V-NN")
  tmp2 <- filter(permuTests_bootCI_j, Marker == markeri & Comparison == "2-P-NN vs. 1-P-NN")
  
  tab <- add_row(.data = tab,"Measure" = "Durability", "Mark" = markeri, 
                 "2-V-NN vs. 1-V-NN" = tmp1$durability_diff, "Unadjusted P1" = tmp1$Pvalue, "FWER P1" = tmp1$FWER_Pvalue,
                 "2-P-NN vs. 1-P-NN" = tmp2$durability_diff, "Unadjusted P2" = tmp2$Pvalue, "FWER P2" = tmp2$FWER_Pvalue,)
}
for(markeri in allMarkers){
  tmp1 <- filter(permuTests_bootCI_j, Marker == markeri & Comparison == "2-V-NN vs. 1-V-NN")
  tmp2 <- filter(permuTests_bootCI_j, Marker == markeri & Comparison == "2-P-NN vs. 1-P-NN")
  
  tab <- add_row(.data = tab,"Measure" = "D202 to D43 GMR", "Mark" = markeri, 
                 "2-V-NN vs. 1-V-NN" = tmp1$D202toD43GMR_diff, "Unadjusted P1" = tmp1$D202toD43GMR_pvalue, "FWER P1" = tmp1$D202toD43GMR_FWER_Pvalue,
                 "2-P-NN vs. 1-P-NN" = tmp2$D202toD43GMR_diff, "Unadjusted P2" = tmp2$D202toD43GMR_pvalue, "FWER P2" = tmp1$D202toD43GMR_FWER_Pvalue)
}
write.csv(tab, file.path(tabDir,"nonnaive_stage2_vs_stage1.csv"), row.names = FALSE)




permuTests_bootCI_j <- filter(output, Comparison %in% c("2-V-NN1 vs. 2-V-NN2", "2-P-NN1 vs. 2-P-NN2") 
                              & Marker %in% allMarkers)

tab <- tibble("Measure" = character(), "Mark" = character(), "2-V-NN1 vs. 2-V-NN2" = character(), "Unadjusted P1" = character(),
              "2-P-NN1 vs. 2-P-NN2" = character(), "Unadjusted P2" = character())
for(markeri in allMarkers){
  tmp1 <- filter(permuTests_bootCI_j, Marker == markeri & Comparison == "2-V-NN1 vs. 2-V-NN2")
  tmp2 <- filter(permuTests_bootCI_j, Marker == markeri & Comparison == "2-P-NN1 vs. 2-P-NN2")
  tab <- add_row(.data = tab,"Measure" = "Day 43", "Mark" = markeri, 
                 "2-V-NN1 vs. 2-V-NN2" = tmp1$D43_diff, "Unadjusted P1" = tmp1$D43_Pvalue,
                 "2-P-NN1 vs. 2-P-NN2" = tmp2$D43_diff, "Unadjusted P2" = tmp2$D43_Pvalue)
  
  
}
for(markeri in allMarkers){
  tmp1 <- filter(permuTests_bootCI_j, Marker == markeri & Comparison == "2-V-NN1 vs. 2-V-NN2")
  tmp2 <- filter(permuTests_bootCI_j, Marker == markeri & Comparison == "2-P-NN1 vs. 2-P-NN2")
  tab <- add_row(.data = tab,"Measure" = "Durability", "Mark" = markeri, 
                 "2-V-NN1 vs. 2-V-NN2" = tmp1$durability_diff, "Unadjusted P1" = tmp1$Pvalue,
                 "2-P-NN1 vs. 2-P-NN2" = tmp2$durability_diff, "Unadjusted P2" = tmp2$Pvalue)
}
for(markeri in allMarkers){
  tmp1 <- filter(permuTests_bootCI_j, Marker == markeri & Comparison == "2-V-NN1 vs. 2-V-NN2")
  tmp2 <- filter(permuTests_bootCI_j, Marker == markeri & Comparison == "2-P-NN1 vs. 2-P-NN2")
  tab <- add_row(.data = tab,"Measure" = "D202 to D43 GMR", "Mark" = markeri, 
                 "2-V-NN1 vs. 2-V-NN2" = tmp1$D202toD43GMR_diff, "Unadjusted P1" = tmp1$D202toD43GMR_pvalue,
                 "2-P-NN1 vs. 2-P-NN2" = tmp2$D202toD43GMR_diff, "Unadjusted P2" = tmp2$D202toD43GMR_pvalue)
}
write.csv(tab, file.path(tabDir,"stage2_nonnaive_GA_vs_GB.csv"), row.names = FALSE)

