#Data processing including
 #a. Exclude participants with missing baseline Naive/Non-naive status
 #b. Impute missing antibody markers
 #c. Detect asymptomatic infections and exclude measurements post primary endpoints and asymptomatic infections
 

rm(list=ls(all=TRUE))
here::i_am("README.md")
repoDir <- here::here()
datDir <- file.path(repoDir,"data")
codeDir <- file.path(repoDir, "code")
outputDir <- file.path(repoDir, "output")

library(tidyverse)
library(plyr)
source(file.path(repoDir, "code/common.R"))


############################################################################################
ppdat <- rbind(read.csv(file.path(datDir, "COVID_Sanofi_stage1_20250312.csv")),
               read.csv(file.path(datDir, "COVID_Sanofi_stage2_20250312.csv")))

longitudinalSubcohortPtid <- read.csv(file.path(datDir, "longitudinalSubcohortPtid.csv"))

#a. exclude participants with missing Bserostatus
ppdat <- filter(ppdat, !is.na(Bserostatus) & Perprotocol == 1)
ppdat$longsubcohort <- 1*ppdat$Subjectid %in% longitudinalSubcohortPtid$x


#b. Impute immune markers to have all-or-none pattern at each time point, saperately for bAb-IgG and nAb-ID50 markers
library(mice) # version 3.17.0
timepoints <- c(0, 22, 43, 78, 134, 202, 292, 387)
imputedppdat <- ppdat


missingTab <- tibble("time" = numeric(),"assay" = character(), "stage" = numeric(), "Bserostatus" = numeric(), "treatment" = numeric(), "total" = numeric(), "occasionalMiss" = numeric(),
                     "ReftiterOnly" = numeric(), "oneMarkerObserved" = numeric())

for(time in timepoints){
  if(time == 0){
    vars1 <- paste("B", c("pseudoneutid50", "pseudoneutid50_BA.1", "pseudoneutid50_BA.2", "pseudoneutid50_BA.4.5",
                                  "pseudoneutid50_B.1.351"), sep = "")
    
    vars2 <- paste("B", c("bindSpike", "bindSpike_beta", "bindSpike_alpha", "bindSpike_gamma", "bindSpike_delta1", "bindSpike_delta2",
                                  "bindSpike_delta3", "bindSpike_omicron"), sep = "")
    
  }else{
    vars1 <- paste("Day", time, c("pseudoneutid50", "pseudoneutid50_BA.1", "pseudoneutid50_BA.2", "pseudoneutid50_BA.4.5",
                                  "pseudoneutid50_B.1.351"), sep = "")
    
    vars2 <- paste("Day", time, c("bindSpike", "bindSpike_beta", "bindSpike_alpha", "bindSpike_gamma", "bindSpike_delta1", "bindSpike_delta2",
                                  "bindSpike_delta3", "bindSpike_omicron"), sep = "")
    
  }
  
  
  for(step in c(1, 2)){
    if(step == 1){
      vars <- vars1
    }else{
      vars <- vars2
    }
    imputeD <- select(ppdat, all_of(c(vars, "Subjectid","Trialstage", "Bserostatus", "Trt", "longsubcohort")))
    newImputeD <- data.frame()
    for(s in c(1, 2)){
      for(b in c(0, 1)){
        for(trt in c(0, 1)){
          datasub <- select(filter(imputeD, Trialstage == s & Bserostatus == b & Trt == trt), all_of(c("Subjectid",vars)))
          longsubcohort <- filter(imputeD, Trialstage == s & Bserostatus == b & Trt == trt)$longsubcohort
          nmiss <- apply(datasub, 1, function(x)sum(is.na(x)))
          nmiss_longsubcohort <- apply(datasub[longsubcohort==1, ], 1, function(x)sum(is.na(x)))
          
          #remove rows that are completely missing
          datasub <- filter(datasub, nmiss < length(vars))
          ntotal <- dim(datasub)[1]
          
          longsubcohort2 <- longsubcohort[nmiss < length(vars)]
          ntotal_longsubcohort <- dim(datasub[longsubcohort2==1, ])[1]
          # if there is no variability, fill in NA with constant values
          for (a in vars) {
            if (all(datasub[, a]==min(datasub[, a], na.rm=TRUE), na.rm=TRUE)) datasub[, a]=min(datasub[, a], na.rm=TRUE)
          } 
          
          if(step == 1  ){
            #identify rows with titer against reference only (99% of rows with one titer)
            kp <- (!is.na(datasub[, 2]) & is.na(datasub[, 3]) & is.na(datasub[, 4]) & is.na(datasub[, 5]) & is.na(datasub[, 6]))
            missingdf <- data.frame(!is.na(datasub[, 2]), !is.na(datasub[, 3]), !is.na(datasub[, 4]), !is.na(datasub[, 5]), !is.na(datasub[, 6]))
            n1 <- apply(missingdf, 1, sum)
            ncomplete <-  sum((!is.na(datasub[, 2]) & !is.na(datasub[, 3]) & !is.na(datasub[, 4]) & !is.na(datasub[, 5]) & !is.na(datasub[, 6])))
            
            kp_longsubcohort <-  (!is.na(datasub[longsubcohort2==1, 2]) & is.na(datasub[longsubcohort2==1, 3]) & 
                                    is.na(datasub[longsubcohort2==1, 4]) & is.na(datasub[longsubcohort2==1, 5]) & 
                                    is.na(datasub[longsubcohort2==1, 6]))
            n1_longsubcohort <- apply(missingdf[longsubcohort2==1,], 1, sum)
            ncomplete_longsubcohort <- sum((!is.na(datasub[longsubcohort2==1, 2]) & !is.na(datasub[longsubcohort2==1, 3]) & 
                                              !is.na(datasub[longsubcohort2==1, 4]) & !is.na(datasub[longsubcohort2==1, 5]) & 
                                              !is.na(datasub[longsubcohort2==1, 6])))
            if(sum(kp < 10)){
              #use MICE to impute all rows
              micey <- mice(datasub[,-1], printFlag = FALSE, diagnostics = FALSE , remove_collinear = FALSE, seed = 2, m = 1)
              datasub[, vars] <- complete(micey, 1)[, vars] 
              
            }else{
              #use MICE to impute rows with more than one titer
              micey <- mice(datasub[!kp,-1], printFlag = FALSE, diagnostics = FALSE , remove_collinear = FALSE, seed = 2, m = 1)
              datasub[!kp, vars] <- complete(micey, 1)[, vars] 
              
              #use Deming regression to impute rows with only one titer
              for(var in c("pseudoneutid50_BA.1", "pseudoneutid50_BA.2", "pseudoneutid50_BA.4.5","pseudoneutid50_B.1.351") ){
                if(time == 0){
                  marker <- paste0("B","pseudoneutid50")
                  nvar <- paste0("B", var)
                }else{
                  marker <- paste0("Day",time,"pseudoneutid50")
                  nvar <- paste0("Day",time, var)
                }
                
                x <- datasub[, marker]
                y <- datasub[, nvar]
                #the cond removes outliers and focus on where the new x will be
                fit = kyotil::Deming(x[x>3 & x-y<2], y[x>3 & x-y<2], boot=F)
                set.seed(1)
                newvalue = predict(fit, datasub[, marker]) + rnorm(mean=0,sd=fit$coef[3],n=nrow(datasub[, ]))
                #previously observed will not be imputed again
                newvalue[!is.na(datasub[, nvar])] = datasub[!is.na(datasub[, nvar]), nvar]
                if (any(is.na(newvalue))) {
                  # if imputation failed, just fill it with median
                  newvalue=y
                  newvalue[is.na(newvalue)]  = median(y,na.rm=T)
                }
                newvalue[newvalue<log10(40)] = log10(20)
                datasub[, nvar] <- newvalue
              }
            }
            missingTab <- add_row(.data = missingTab, "time" = time,"assay" = ifelse(step==1, "nAb", "bAb"), "stage" = s, 
                                  "Bserostatus" = b, "treatment" = trt, "total" = ntotal_longsubcohort, 
                                  "occasionalMiss" = ntotal_longsubcohort - sum(n1_longsubcohort==1) - ncomplete_longsubcohort,
                                 "ReftiterOnly" = sum(kp_longsubcohort), "oneMarkerObserved" = sum(n1_longsubcohort==1))
            
            
          }else{
            ncomplete <-  sum((!is.na(datasub[, 2]) & !is.na(datasub[, 3]) & !is.na(datasub[, 4]) & !is.na(datasub[, 5]) & !is.na(datasub[, 6])
                               & !is.na(datasub[, 7]) & !is.na(datasub[, 8]) & !is.na(datasub[, 9])))
            missingdf <- data.frame(!is.na(datasub[, 2]), !is.na(datasub[, 3]), !is.na(datasub[, 4]), !is.na(datasub[, 5]), !is.na(datasub[, 6]), !is.na(datasub[, 7]), !is.na(datasub[, 8]) , !is.na(datasub[, 9]))
            n1 <- apply(missingdf, 1, sum)
            
            ncomplete_longsubcohort <- sum((!is.na(datasub[longsubcohort2==1, 2]) & !is.na(datasub[longsubcohort2==1, 3]) & 
                                              !is.na(datasub[longsubcohort2==1, 4]) & !is.na(datasub[longsubcohort2==1, 5]) & 
                                              !is.na(datasub[longsubcohort2==1, 6]) & !is.na(datasub[longsubcohort2==1, 7]) &
                                              !is.na(datasub[longsubcohort2==1, 8]) & !is.na(datasub[longsubcohort2==1, 9])))
            n1_longsubcohort <- apply(missingdf[longsubcohort2==1,], 1, sum)
            missingTab <- add_row(.data = missingTab, "time" = time,"assay" = ifelse(step==1, "nAb", "bAb"), "stage" = s, 
                                  "Bserostatus" = b, "treatment" = trt, "total" = ntotal_longsubcohort, 
                                  "occasionalMiss" = ntotal_longsubcohort  - ncomplete_longsubcohort - sum(n1_longsubcohort == 1),
                                  "ReftiterOnly" = NA, "oneMarkerObserved"= sum(n1_longsubcohort == 1))
            micey <- mice(datasub[,-1], printFlag = FALSE, diagnostics = FALSE , remove_collinear = FALSE, seed = 2, m = 1)
            datasub[, vars] <- complete(micey, 1)[, vars]
            
          }
          newImputeD <- rbind(newImputeD, datasub)
        }
      }
    }
    #replace columns with imputed data
    rows <- laply(newImputeD$Subjectid, function(x)c(1:length(imputedppdat$Subjectid))[imputedppdat$Subjectid == x])
    imputedppdat[rows, vars] <- newImputeD[, vars]
    
  }
 
  
}

write.csv(missingTab, file.path(outputDir, "missingTab.csv"))



#c. create a longitudinal dataset for asymptomatic infection detection
bindingAntibodyMarkers <- c("BbindSpike","Day22bindSpike", "Day43bindSpike","Day78bindSpike", "Day134bindSpike", "Day202bindSpike", "Day292bindSpike", "Day387bindSpike",
                            "BbindSpike_beta","Day22bindSpike_beta","Day43bindSpike_beta","Day78bindSpike_beta", "Day134bindSpike_beta", "Day202bindSpike_beta", "Day292bindSpike_beta", "Day387bindSpike_beta",
                            "BbindSpike_alpha","Day22bindSpike_alpha","Day43bindSpike_alpha","Day78bindSpike_alpha", "Day134bindSpike_alpha", "Day202bindSpike_alpha", "Day292bindSpike_alpha", "Day387bindSpike_alpha",
                            "BbindSpike_gamma","Day22bindSpike_gamma","Day43bindSpike_gamma","Day78bindSpike_gamma", "Day134bindSpike_gamma", "Day202bindSpike_gamma", "Day292bindSpike_gamma", "Day387bindSpike_gamma",
                            "BbindSpike_delta1","Day22bindSpike_delta1","Day43bindSpike_delta1","Day78bindSpike_delta1", "Day134bindSpike_delta1", "Day202bindSpike_delta1", "Day292bindSpike_delta1", "Day387bindSpike_delta1",
                            "BbindSpike_delta2","Day22bindSpike_delta2","Day43bindSpike_delta2","Day78bindSpike_delta2", "Day134bindSpike_delta2", "Day202bindSpike_delta2", "Day292bindSpike_delta2", "Day387bindSpike_delta2",
                            "BbindSpike_delta3","Day22bindSpike_delta3","Day43bindSpike_delta3","Day78bindSpike_delta3", "Day134bindSpike_delta3", "Day202bindSpike_delta3", "Day292bindSpike_delta3", "Day387bindSpike_delta3",
                            "BbindSpike_omicron","Day22bindSpike_omicron","Day43bindSpike_omicron","Day78bindSpike_omicron", "Day134bindSpike_omicron", "Day202bindSpike_omicron", "Day292bindSpike_omicron", "Day387bindSpike_omicron")

neutAntibodyMarkers <- c("Bpseudoneutid50", "Day22pseudoneutid50", "Day43pseudoneutid50", "Day78pseudoneutid50","Day134pseudoneutid50", "Day202pseudoneutid50", "Day292pseudoneutid50", "Day387pseudoneutid50",
                         "Bpseudoneutid50_BA.1", "Day22pseudoneutid50_BA.1","Day43pseudoneutid50_BA.1", "Day78pseudoneutid50_BA.1","Day134pseudoneutid50_BA.1", "Day202pseudoneutid50_BA.1", "Day292pseudoneutid50_BA.1", "Day387pseudoneutid50_BA.1",
                         "Bpseudoneutid50_BA.2", "Day22pseudoneutid50_BA.2","Day43pseudoneutid50_BA.2", "Day78pseudoneutid50_BA.2","Day134pseudoneutid50_BA.2", "Day202pseudoneutid50_BA.2", "Day292pseudoneutid50_BA.2", "Day387pseudoneutid50_BA.2",
                         "Bpseudoneutid50_BA.4.5", "Day22pseudoneutid50_BA.4.5","Day43pseudoneutid50_BA.4.5", "Day78pseudoneutid50_BA.4.5","Day134pseudoneutid50_BA.4.5", "Day202pseudoneutid50_BA.4.5", "Day292pseudoneutid50_BA.4.5", "Day387pseudoneutid50_BA.4.5",
                         "Bpseudoneutid50_B.1.351", "Day22pseudoneutid50_B.1.351","Day43pseudoneutid50_B.1.351", "Day78pseudoneutid50_B.1.351","Day134pseudoneutid50_B.1.351", "Day202pseudoneutid50_B.1.351", "Day292pseudoneutid50_B.1.351", "Day387pseudoneutid50_B.1.351")
                      
antibodyMarkers <- c(bindingAntibodyMarkers, neutAntibodyMarkers)
vars <- c("Ptid","Trt", "Trialstage", "Sex", "Age", "BMI", "Age", "Country",  "Bserostatus", 
          "EventTimeFirstInfectionD1", "EventIndFirstInfectionD1",
          "EventTimeSecondInfectionD1", "EventIndSecondInfectionD1",
          "EventTimeFirstInfectionD43", "EventIndFirstInfectionD43", 
          "FirstEnrollmentDate", "Perprotocol", 
          "D01_S_pos_only_in_non_naive_group","NumberdaysD1toD22", "NumberdaysD1toD43", 
          "NumberdaysD1toD78", "NumberdaysD1toD134", "NumberdaysD1toD202", "NumberdaysD1toD292",
          "NumberdaysD1toD387", "D01_S_pos", "D01_N_pos", "D01_NAAT_pos", "D22_N_pos", "D22_NAAT_pos")

allvars <- antibodyMarkers
colnames(imputedppdat)[1] <- "Ptid"
df <- select(imputedppdat, all_of(c(vars, antibodyMarkers)))
df1 <- filter(df, !if_all(allvars, is.na))

longdf <- tibble("Ptid" = character(), "visitn" = character(), "time" = numeric(), "bindSpike" = numeric(), "bindSpike_beta" = numeric(), "bindSpike_alpha" = numeric(), 
                 "bindSpike_gamma" = numeric(), "bindSpike_delta1" = numeric(), "bindSpike_delta2" = numeric(), "bindSpike_delta3" = numeric(), "bindSpike_omicron" = numeric(),
                 "pseudoneutid50" = numeric(), "pseudoneutid50_BA.1" = numeric(), "pseudoneutid50_BA.2" = numeric(), "pseudoneutid50_BA.4.5" = numeric(), "pseudoneutid50_B.1.351" = numeric())

for(id in unique(df1$Ptid)){
  for(t in c(0, 22, 43, 78, 134, 202, 292, 387)){
    markernames <- antibodyMarkers[grepl(t, antibodyMarkers)]
    dfi <- filter(df1, Ptid == id)
    visit <- paste0("Day", t)
    if(t == 0){
      time <- 0
    }else{
      time <- dfi[, paste0("NumberdaysD1toD", t)]
    }
    
    if(t == 0){
      visit = "B"
      longdf <- add_row(longdf, Ptid  = id, visitn = "Baseline", time = time, 
                        bindSpike = dfi[, paste0(visit, "bindSpike")], "bindSpike_beta" = dfi[, paste0(visit, "bindSpike_beta")], "bindSpike_alpha" = dfi[, paste0(visit, "bindSpike_alpha")], 
                        "bindSpike_gamma" = dfi[, paste0(visit, "bindSpike_gamma")], "bindSpike_delta1" = dfi[, paste0(visit, "bindSpike_delta1")], 
                        "bindSpike_delta2" = dfi[, paste0(visit, "bindSpike_delta2")], "bindSpike_delta3" = dfi[, paste0(visit, "bindSpike_delta3")], 
                        "bindSpike_omicron" = dfi[, paste0(visit, "bindSpike_omicron")],
                        "pseudoneutid50" = dfi[, paste0(visit, "pseudoneutid50")], "pseudoneutid50_BA.1" = dfi[, paste0(visit, "pseudoneutid50_BA.1")], 
                        "pseudoneutid50_BA.2" = dfi[, paste0(visit, "pseudoneutid50_BA.2")], "pseudoneutid50_BA.4.5" = dfi[, paste0(visit, "pseudoneutid50_BA.4.5")], 
                        "pseudoneutid50_B.1.351" = dfi[, paste0(visit, "pseudoneutid50_B.1.351")])
    }else{
      longdf <- add_row(longdf, Ptid  = id, visitn = visit, time = time, 
                        bindSpike = dfi[, paste0(visit, "bindSpike")], "bindSpike_beta" = dfi[, paste0(visit, "bindSpike_beta")], "bindSpike_alpha" = dfi[, paste0(visit, "bindSpike_alpha")], 
                        "bindSpike_gamma" = dfi[, paste0(visit, "bindSpike_gamma")], "bindSpike_delta1" = dfi[, paste0(visit, "bindSpike_delta1")], 
                        "bindSpike_delta2" = dfi[, paste0(visit, "bindSpike_delta2")], "bindSpike_delta3" = dfi[, paste0(visit, "bindSpike_delta3")], 
                        "bindSpike_omicron" = dfi[, paste0(visit, "bindSpike_omicron")],
                        "pseudoneutid50" = dfi[, paste0(visit, "pseudoneutid50")], "pseudoneutid50_BA.1" = dfi[, paste0(visit, "pseudoneutid50_BA.1")], 
                        "pseudoneutid50_BA.2" = dfi[, paste0(visit, "pseudoneutid50_BA.2")], "pseudoneutid50_BA.4.5" = dfi[, paste0(visit, "pseudoneutid50_BA.4.5")], 
                        "pseudoneutid50_B.1.351" = dfi[, paste0(visit, "pseudoneutid50_B.1.351")])
    }
    
  }
}

write.csv(df1, file.path(datDir, "antibodyDurabilityWide_imputed.csv"), row.names = FALSE)
write.csv(longdf, file.path(datDir, "antibodyDurabilityLong_imputed.csv"), row.names = FALSE)


#remove observations after symptomatic and asymptomatic infections
fulllongdf <- left_join(longdf, df1, by = "Ptid")
fulllongdf$visitn[fulllongdf$visitn=="Baseline"]<-"Day0"
fulllongdf$time2 <- as.numeric(gsub("Day","",fulllongdf$visitn))

fulllongdf$BserostatusD01_S <- ifelse(fulllongdf$Bserostatus==0, "naive", 
                                      ifelse(fulllongdf$D01_S_pos_only_in_non_naive_group == 1, 
                                             "nonnaiveD01SposOnly", "nonnaiveElse"))


#x: a matrix of individual longitudinal data for bAb-IgG Delta B.1.617.2/AY.4, bAb-IgG Omicron, and nAb-ID50 Omicron BA.4/5
#time: a vector of time corresponding to the longitudinal data
#time.begin: the beginning time point to start asymptomatic infection detection 
#markers: the variable name of the markers
#fold: a vector of scalars; the first three are for fold-increases in bAb-IgG Delta B.1.617.2/AY.4, bAb-IgG Omicron, and 
#nAb-ID50 Omicron BA.4/5 needed to define an asymptomatic infection when starting from above LOD/LLOQ
#the fourth scalar is for fold-increase in serum conversion 
#output: indicators of whether the marker has had an significant increase according to the criterion at each time point; 
#0 for not asymptomatic infection detected and hence excluded
#1 for asymptomatic infection not detected and hence included
markerIncrease<-function(x, time, time.begin, markers, fold){
  increaseInd <- matrix(ncol = dim(x)[2], nrow = length(time))
  includeInd <- matrix(ncol = dim(x)[2], nrow = length(time))
  seq <- 1:length(time)
  increaseInd[seq[time <= time.begin], ] <- 0
  i.max <- max(seq[time <= time.begin])
  for(i in 1:i.max){
    includeInd[i, ] <- ifelse(is.na(x[i, ]), NA, 1)
  }
  
  if(i.max < length(time)){
    if(length(time) > 1){
      for(j in 1:dim(x)[2]){
        for(i in (i.max+1):length(time)){
          # if at time i, the marker is missing, then indication of inclusion is 0
          # if at time i, the marker is not missing but the marker is missing before time i, include time i
          # if at time i, the marker is not missing and there is observed marker before time i, exclude time i if there is x fold increase, or cross from before LLOQ to above fold[4]*LLOQ comparing time i 
          # (if there were below LLOQ, they were set as approximately LLOQ/2)
          # and the minimum before time i. It could be that marker at time i is much greater than it at time i-2 but not much than it at time i-1. This still indicates an infection already happened.
          # if at time i, the marker is not missing and there is observed marker before time i, exclude time i if before time i, there is x fold increase, or cross from before LLOQ to above LLOQ; 
          # i.e., if asymptomatic infections happened, the inclusion indicator is zero 
          
          if(is.na(x[i,j])){ 
            increaseInd[i, j] <- NA
            includeInd[i, j] <- NA
          }else if(sum(!is.na(x[1:(i-1), j]), na.rm = TRUE) >0){ # there is observation in the past
            
            startObs <- min(x[1:(i-1), j], na.rm = TRUE)
            LLOQmarker <- LLOQf(markers[j])
            if((startObs < LLOQmarker & x[i,j] >= fold[4]*LLOQmarker) | (startObs >= LLOQmarker & x[i, j]> (startObs + log10(fold[j])))){
              increaseInd[i, j] <- 1 
              includeInd[i, j] <- 0
            }else{
              increaseInd[i, j] <- 0
              includeInd[i, j] <- ifelse(sum(!is.na(increaseInd[1:(i-1),j]))>0 && sum(increaseInd[1:(i-1),j], na.rm = TRUE) >0, 0, 1) 
            }
            
          }else{ 
            increaseInd[i, j] <- 0
            includeInd[i, j] <- 1
          }
          
        }
      }
    }
  }
  
  
  includeInd <- data.frame(includeInd)
  colnames(includeInd) <- paste(colnames(x), "ind", sep = "_")
  return(includeInd)
}

#tmpobs: a matrix of time and individual longitudinal data for bAb-IgG Delta B.1.617.2/AY.4, bAb-IgG Omicron, and nAb-ID50 Omicron BA.4/5
#dat_proc_tmp: the dataset to which indicators of asymptimatic detection be added to
#time.begin: the beginning time of screening for asymptomatic infections
#foldx: a vector of scalars; the first three are for fold-increases in bAb-IgG Delta B.1.617.2/AY.4, bAb-IgG Omicron, and 
#nAb-ID50 Omicron BA.4/5 needed to define an asymptomatic infection when starting from above LOD/LLOQ
asymptomaticInf <- function(tmpobs, dat_proc_tmp, time.begin, foldx){
  if(dim(tmpobs)[1] >=1){
    #see if any sample should be included at each time point
    includeIndAllmarkers <- markerIncrease(select(tmpobs,all_of(c("bindSpike_delta3", "bindSpike_omicron", "pseudoneutid50_BA.4.5"))), tmpobs$time2, time.begin,
                                           markers = c("bindSpike_delta3", "bindSpike_omicron", "pseudoneutid50_BA.4.5"), fold = foldx)
    
    includeInd <- includeIndAllmarkers
    includeSum <- apply(includeInd, 1, function(x) sum(x == 1, na.rm = TRUE))
    nmark <- apply(includeInd, 1, function(x) sum(!is.na(x)))
    tmpobsfinal <- tmpobs[includeSum !=0 & includeSum/nmark == 1, ] #if all available markers approve including, then include, i.e, if any available marker indicates asymptomatic infection, exclude
    
    
    if(dim(tmpobsfinal)[1] >=1){
      dat_proc_tmp$missingDueToAsymInfection <- NA
      dat_proc_tmp$missingDueToAsymInfection[dat_proc_tmp$time %in% tmpobsfinal$time] <- 0
      dat_proc_tmp$missingDueToAsymInfection[!dat_proc_tmp$time %in% tmpobsfinal$time & 
                                                dat_proc_tmp$time %in% tmpobs$time] <- 1
    }else{
      dat_proc_tmp$missingDueToAsymInfection <- NA
      dat_proc_tmp$missingDueToAsymInfection[dat_proc_tmp$time %in% tmpobs$time] <- 1
      
    }
    
  }else{
    dat_proc_tmp$missingDueToAsymInfection <- NA
  }
  return(dat_proc_tmp)
}


#The four version of the algorithm considered
foldx <- c(4, 4, 4, 1)


dffinal <- data.frame()
dat_proc <- data.frame()
for(ptid in unique(fulllongdf$Ptid)){
  tmp <- filter(fulllongdf, Ptid == ptid)
  dat_proc_tmp <- tmp
  #label observations became cases after visit D22 
  if(tmp$EventIndFirstInfectionD1[1] == 1 & tmp$EventTimeFirstInfectionD1[1] > filter(tmp, visitn == "Day22")$time){
    tmp <- filter(tmp, time <= EventTimeFirstInfectionD1)
  }
  
  dat_proc_tmp$missingDueToSymInfection <- NA
  dat_proc_tmp$missingDueToSymInfection[dat_proc_tmp$time %in% tmp$time] <- 0
  dat_proc_tmp$missingDueToSymInfection[!dat_proc_tmp$time %in% tmp$time] <- 1
  
  #remove obs post asymptomatic infections after D43 for vaccine group
  #remove obs post asymptomatic infections after D22 for naive placebo group
  #remove obs post asymptomatic infections after D43 for non-naive placebo group
  
  if(tmp$Trt[1] == 1){
    tmpobs <- filter(tmp, time2>=43)
    dat_proc_tmp <- asymptomaticInf(tmpobs, dat_proc_tmp, time.begin = 43,foldx = foldx)
    dat_proc <- rbind(dat_proc, dat_proc_tmp)
  }
  
  if(tmp$Trt[1] == 0 & tmp$Bserostatus[1] == 1){
    tmpobs <- filter(tmp, time2>=43)
    dat_proc_tmp <- asymptomaticInf(tmpobs, dat_proc_tmp, time.begin = 43, foldx = foldx)
    dat_proc <- rbind(dat_proc, dat_proc_tmp)
  }
  
  if(tmp$Trt[1] == 0 & tmp$Bserostatus[1] == 0){
    tmpobs <- filter(tmp, time2>=22)
    dat_proc_tmp <- asymptomaticInf(tmpobs, dat_proc_tmp, time.begin= 22, foldx = foldx)
    dat_proc <- rbind(dat_proc, dat_proc_tmp)
  }
  
}

dat_proc <- distinct(dat_proc)
dat_proc$missingDueToSymInfection_nAb <- 1*(dat_proc$missingDueToSymInfection & !is.na(dat_proc$pseudoneutid50_BA.4.5))
dat_proc$missingDueToSymInfection_bAb <- 1*(dat_proc$missingDueToSymInfection & !is.na(dat_proc$bindSpike_omicron))
dat_proc$missingDueToAsymInfection_nAb <- 1*(dat_proc$missingDueToAsymInfection & !is.na(dat_proc$pseudoneutid50_BA.4.5))
dat_proc$missingDueToAsymInfection_bAb <- 1*(dat_proc$missingDueToAsymInfection & !is.na(dat_proc$bindSpike_omicron))


#remove observations post detection of infections
markers1 <- c("bindSpike", "bindSpike_beta", "bindSpike_alpha", "bindSpike_gamma", "bindSpike_delta1", 
              "bindSpike_delta2", "bindSpike_delta3", "bindSpike_omicron")
markers2 <- c("pseudoneutid50", "pseudoneutid50_B.1.351", "pseudoneutid50_BA.1", "pseudoneutid50_BA.2", "pseudoneutid50_BA.4.5")

longitudinalfinal_nAb <- select(filter(dat_proc, !(missingDueToSymInfection_nAb == 1 | missingDueToAsymInfection_nAb == 1 | is.na(pseudoneutid50_BA.4.5)) & time2 >=43), 
                                all_of(c(vars, "time2", markers2, paste("B", markers2, sep = ""), paste("Day22", markers2, sep = "") )))
longitudinalfinal_bAb <- select(filter(dat_proc, !(missingDueToSymInfection_bAb == 1 | missingDueToAsymInfection_bAb == 1 | is.na(bindSpike_omicron))& time2 >=43), 
                                all_of(c(vars, "time2", markers1, paste("B", markers1, sep = ""), paste("Day22", markers1, sep = "") )))


#remove observations before D43 and observations with less than 3 obs
longitudinalfinal_nAb <- filter(longitudinalfinal_nAb, time2 >=43)
longitudinalfinal_bAb <- filter(longitudinalfinal_bAb, time2 >=43)
ptids_nAb <- table(longitudinalfinal_nAb$Ptid)
ptids_bAb <- table(longitudinalfinal_bAb$Ptid)
longitudinalfinal_nAb <- filter(longitudinalfinal_nAb, Ptid %in% names(ptids_nAb)[ptids_nAb >=3])
longitudinalfinal_bAb <- filter(longitudinalfinal_bAb, Ptid %in% names(ptids_bAb)[ptids_bAb >=3])


#standardize continuous covariates
nAbdata <- longitudinalfinal_nAb 
bAbdata <- longitudinalfinal_bAb 
nAbdata0 <- filter(longitudinalfinal_nAb, Bserostatus == 1)
bAbdata0 <- filter(longitudinalfinal_bAb, Bserostatus == 1) 
nAbdata1 <- nAbdata
bAbdata1 <- bAbdata
scalef <- function(var, varname, data1, data2 = NULL){
  if(!is.null(data2)){
    data <- rbind(data1[, c("Ptid", varname)], data2[, c("Ptid", varname)])
    uniqueptids <- unique(data$Ptid)
    varunique <- laply(uniqueptids, function(x) as.numeric(unique(filter(data, Ptid == x)[, varname])))
    mean <- mean(varunique, na.rm = TRUE)
    sd <- sd(varunique, na.rm = TRUE)
    return(list((var - mean)/sd, mean, sd))
  }else{
    data <- data1[, c("Ptid", varname)]
    uniqueptids <- unique(data$Ptid)
    varunique <- laply(uniqueptids, function(x) as.numeric(unique(filter(data, Ptid == x)[, varname])))
    mean <- mean(varunique, na.rm = TRUE)
    sd <- sd(varunique, na.rm = TRUE)
    return(list((var - mean)/sd, mean, sd))
  }
  
}
agesummary <- scalef(nAbdata1$Age,"Age", nAbdata1,bAbdata1)[2:3]
nAbdata$Age0 <- nAbdata$Age
nAbdata$BMI0 <- nAbdata$BMI
bAbdata$Age0 <- bAbdata$Age
bAbdata$BMI0 <- bAbdata$BMI

nAbdata$Age <- as.vector(scalef(nAbdata1$Age,"Age", nAbdata1,bAbdata1)[[1]])
nAbdata$BMI <- as.vector(scalef(nAbdata1$BMI,"BMI", nAbdata1,bAbdata1)[[1]])

bAbdata$Age <- as.vector(scalef(bAbdata1$Age,"Age", nAbdata1,bAbdata1)[[1]])
bAbdata$BMI <- as.vector(scalef(bAbdata1$BMI,"BMI", nAbdata1,bAbdata1)[[1]])


#keep the unstandardized day 1 marker values with variables appended by _0
nAbdata$Bpseudoneutid50_0 <- nAbdata$Bpseudoneutid50
nAbdata$Bpseudoneutid50_B.1.351_0 <-  nAbdata$Bpseudoneutid50_B.1.351
nAbdata$Bpseudoneutid50_BA.1_0 <- nAbdata$Bpseudoneutid50_BA.1
nAbdata$Bpseudoneutid50_BA.2_0 <- nAbdata$Bpseudoneutid50_BA.2
nAbdata$Bpseudoneutid50_BA.4.5_0 <- nAbdata$Bpseudoneutid50_BA.4.5

nAbdata$Bpseudoneutid50 <- as.vector(scalef(nAbdata$Bpseudoneutid50, "Bpseudoneutid50", nAbdata0)[[1]])
nAbdata$Bpseudoneutid50_B.1.351 <- as.vector(scalef(nAbdata$Bpseudoneutid50_B.1.351, "Bpseudoneutid50_B.1.351", nAbdata0)[[1]])
nAbdata$Bpseudoneutid50_BA.1 <- as.vector(scalef(nAbdata$Bpseudoneutid50_BA.1, "Bpseudoneutid50_BA.1", nAbdata0)[[1]])
nAbdata$Bpseudoneutid50_BA.2 <- as.vector(scalef(nAbdata$Bpseudoneutid50_BA.2, "Bpseudoneutid50_BA.2", nAbdata0)[[1]])
nAbdata$Bpseudoneutid50_BA.4.5 <- as.vector(scalef(nAbdata$Bpseudoneutid50_BA.4.5, "Bpseudoneutid50_BA.4.5", nAbdata0)[[1]])


bAbdata$BbindSpike_0 <- bAbdata$BbindSpike
bAbdata$BbindSpike_alpha_0 <- bAbdata$BbindSpike_alpha
bAbdata$BbindSpike_beta_0 <- bAbdata$BbindSpike_beta
bAbdata$BbindSpike_gamma_0 <- bAbdata$BbindSpike_gamma
bAbdata$BbindSpike_delta1_0 <- bAbdata$BbindSpike_delta1
bAbdata$BbindSpike_delta2_0 <- bAbdata$BbindSpike_delta2
bAbdata$BbindSpike_delta3_0 <- bAbdata$BbindSpike_delta3
bAbdata$BbindSpike_omicron_0 <- bAbdata$BbindSpike_omicron


bAbdata$BbindSpike <- as.vector(scalef(bAbdata$BbindSpike, "BbindSpike", bAbdata0)[[1]])
bAbdata$BbindSpike_alpha <- as.vector(scalef(bAbdata$BbindSpike_alpha, "BbindSpike_alpha", bAbdata0)[[1]])
bAbdata$BbindSpike_beta <- as.vector(scalef(bAbdata$BbindSpike_beta, "BbindSpike_beta", bAbdata0)[[1]])
bAbdata$BbindSpike_gamma <- as.vector(scalef(bAbdata$BbindSpike_gamma, "BbindSpike_gamma", bAbdata0)[[1]])
bAbdata$BbindSpike_delta1 <- as.vector(scalef(bAbdata$BbindSpike_delta1, "BbindSpike_delta1", bAbdata0)[[1]])
bAbdata$BbindSpike_delta2 <- as.vector(scalef(bAbdata$BbindSpike_delta2, "BbindSpike_delta2", bAbdata0)[[1]])
bAbdata$BbindSpike_delta3 <- as.vector(scalef(bAbdata$BbindSpike_delta3, "BbindSpike_delta3", bAbdata0)[[1]])
bAbdata$BbindSpike_omicron <- as.vector(scalef(bAbdata$BbindSpike_omicron, "BbindSpike_omicron", bAbdata0)[[1]])
write.csv(dat_proc, file.path(datDir, paste0("vat08_combined_data_processed_withmissingIndicator.csv")), row.names = FALSE)
write.csv(nAbdata, file.path(datDir, paste0("vat08_combined_data_processed_longitudinal_nAb.csv")), row.names = FALSE)
write.csv(bAbdata, file.path(datDir, paste0("vat08_combined_data_processed_longitudinal_bAb.csv")), row.names = FALSE)

