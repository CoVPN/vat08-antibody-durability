#function to fit the best fraction polynomial model 
#return AUC/(tmax - tmin) and the variance of it
BFpolynomialFit <- function(lmmdata, baselineAdj = FALSE){

  nptid <- length(unique(lmmdata$Ptid))
  op <- options(warn=2)
  lmmdata$time2 <- (lmmdata$time2/43)
  aictable <- tibble("p" = numeric(), "aicv" = numeric(), "random effect" = character(), "warning" = numeric())
  fraction <- c(-2, -1.5, -1, -0.5, 0, 0.5, 1)
  for(i in 1:length(fraction)){
    piter = fraction[i]
    if(piter == 0){
      lmmdata$tftime <- log(lmmdata$time2)
    }else{
      lmmdata$tftime <- lmmdata$time2^piter
    }
    lmmdata$Bmark1 <- ifelse(lmmdata$Bmark <=0, lmmdata$Bmark, 0)
    if(baselineAdj){
      if(min(table(lmmdata$Sex)) <10){ # roughly 3 participants with minority sex
        fm1 <- try(lmer(mark ~ tftime  + (1 |Ptid) + Bmark  + Age + BMI , lmmdata, REML=FALSE))
        fm2 <- try(lmer(mark ~ tftime  + (tftime |Ptid)+ Bmark  + Age + BMI  , lmmdata, REML=FALSE))
      }else{
        fm1 <- try(lmer(mark ~ tftime  + (1 |Ptid) + Bmark  + Age + Sex + BMI, lmmdata, REML=FALSE))
        fm2 <- try(lmer(mark ~ tftime  + (tftime |Ptid)+ Bmark + Age + Sex + BMI , lmmdata, REML=FALSE))
      }
      
    }else{
      if(min(table(lmmdata$Sex)) <10){
        fm1 <- try(lmer(mark ~ tftime  + (1 |Ptid) + Age + BMI , lmmdata, REML=FALSE))
        fm2 <- try(lmer(mark ~ tftime  + (tftime |Ptid) + Age + BMI , lmmdata, REML=FALSE))
      }else{
        fm1 <- try(lmer(mark ~ tftime  + (1 |Ptid) + Age + Sex + BMI, lmmdata, REML=FALSE))
        fm2 <- try(lmer(mark ~ tftime  + (tftime |Ptid) + Age + Sex + BMI , lmmdata, REML=FALSE))
      }
      
    }  
      
    
    if((is(fm1, "try-error"))){
      aictable <- add_row(.data = aictable, "p" = piter, "random effect" = "random intercept", "aicv" = NA,
                          "warning" = 1)
    }else{
      fmsummary <- summary(fm1)
      aictable <- add_row(.data = aictable, "p" = piter, "random effect" = "random intercept", "aicv" = fmsummary$AICtab[1],
                          "warning" = 0)
    }
    if((is(fm2, "try-error"))){
      aictable <- add_row(.data = aictable, "p" = piter, "random effect" = "random intercept and slope", "aicv" = NA,
                          "warning" = 1)
    }else{
      fmsummary <- summary(fm2)
      aictable <- add_row(.data = aictable, "p" = piter, "random effect" = "random intercept and slope", "aicv" = fmsummary$AICtab[1],
                          "warning" = 0)
    }
  }
  options(warn = 0)
  a <- filter(aictable, aicv ==(min(aictable$aicv, na.rm = TRUE)))
  
  if(length(a$p) >0) {
    if(a$p[1] == 0){
      lmmdata$tftime <- log(lmmdata$time2)
    }else{
      lmmdata$tftime <- lmmdata$time2^(a$p[1])
    }
    if(a$`random effect`[1] == "random intercept"){
      if(baselineAdj){
        if(min(table(lmmdata$Sex)) <10){
          fm <- try(lmer(mark ~ tftime  + (1|Ptid)+ Bmark  + Age + BMI, lmmdata, REML=TRUE))
        }else{
          fm <- try(lmer(mark ~ tftime  + (1|Ptid)+ Bmark  + Age + Sex + BMI , lmmdata, REML=TRUE))
        }
        
      }else{
        if(min(table(lmmdata$Sex)) <10){
          fm <- try(lmer(mark ~ tftime  + (1 |Ptid) + Age + BMI, lmmdata, REML=TRUE))
        }else{
          fm <- try(lmer(mark ~ tftime  + (1 |Ptid) + Age + Sex + BMI, lmmdata, REML=TRUE))
        }
        
      } 
    }else{
      if(baselineAdj){
        if(min(table(lmmdata$Sex)) <10){
          fm <- try(lmer(mark ~ tftime  + (tftime |Ptid)+ Bmark  + Age + BMI, lmmdata, REML=TRUE))
        }else{
          fm <- try(lmer(mark ~ tftime  + (tftime |Ptid)+ Bmark   + Age + Sex + BMI, lmmdata, REML=TRUE))
        }
        
      }else{
        if(min(table(lmmdata$Sex)) <10){
          fm <- try(lmer(mark ~ tftime  + (tftime |Ptid) + Age + BMI, lmmdata, REML=TRUE))
        }else{
          fm <- try(lmer(mark ~ tftime  + (tftime |Ptid) + Age + Sex + BMI, lmmdata, REML=TRUE))
        }
        
      } 
    }
   
    fmsummary <- summary(fm)
    
    return(list(fmsummary = fmsummary, p = a$p[1], nptid = nptid,
                fm = fm))
    
  }else{
    return(list(fmsummary = NA, p = NA, nptid = nptid, fm = NA))
  }
}


AUCsummaryConditional<- function(fmsummary, p){
  if(!is.na(p)){
    coefest <- as.numeric(fmsummary$coefficients[, "Estimate"])
    names(coefest) <- rownames(fmsummary$coefficient)
    vcov <- fmsummary$vcov
    names(coefest) <- rownames(fmsummary$coefficient)
    rownames(vcov) <- rownames(fmsummary$coefficient)
    colnames(vcov) <- rownames(fmsummary$coefficient)
    
    if(!is.na(coefest["Sex"])){
      beta0 = coefest["(Intercept)"]
      beta1 = coefest["tftime"]
      beta2 = coefest["Sex"]
      covbeta = vcov[c("(Intercept)", "tftime", "Sex"), c("(Intercept)", "tftime", "Sex")]
      # varbeta0 = vcov["(Intercept)", "(Intercept)"]
      # varbeta1 = vcov["tftime", "tftime"]
      # covbeta0beta1 = vcov["(Intercept)", "tftime"]
      t <- seq(43, 387, 1)
      plott <- t/43
      t1 = 43
      t2 = 202
      if(p == 0){
        predFemale <- beta0 + beta1*log(plott)
        predMale <- beta0 + beta2 + beta1*log(plott)
        a1 = (t2 - t1) 
        a2 = (t2*log(t2) - t1*log(t1) - t2 + t1 - log(t1)*(t2 - t1))
        
      }else if (p == (-1)){
        predFemale <- beta0 + beta1*plott^p
        predMale <- beta0 + beta2 + beta1*plott^p
        a1 = (t2 - t1) 
        a2 = t1*log(t2/t1)
      }else {
        predFemale <- beta0 + beta1*plott^p
        predMale <- beta0 + beta2 + beta1*plott^p
        a1 = (t2 - t1) 
        a2 = 1/t1^p*(t2^(p+1)/(p+1) - t1^(p+1)/(p+1))
      }
      
      AUCestFemale = beta0*a1+ beta1*a2
      AUCestnormalizedFemale = AUCestFemale/(t2-t1)
      AUCnormalizedVarFemale = (covbeta[1, 1]*a1^2 + covbeta[2, 2]*a2^2 + 2*covbeta[1, 2]*a1*a2)/(t2-t1)^2
      
      AUCestMale = (beta0+beta2)*a1+ beta1*a2
      AUCestnormalizedMale = AUCestMale/(t2-t1)
      AUCnormalizedVarMale = ((covbeta[1, 1] + covbeta[3, 3] + 2*covbeta[1, 3])*a1^2 + covbeta[2, 2]*a2^2 + 2*(covbeta[1, 2] + covbeta[2, 3])*a1*a2)/(t2-t1)^2
      est <- c(AUCestnormalizedFemale, AUCnormalizedVarFemale, AUCestnormalizedMale, AUCnormalizedVarMale)
      names(est) <- c("AUC_Female", "Var_Female", "AUC_Male", "Var_Male")
      pred = cbind(t, predFemale, predMale)
      colnames(pred) <- c("time","Female", "Male")
      return(list(est = est, pred = pred))
    }else{
      beta0 = coefest["(Intercept)"]
      beta1 = coefest["tftime"]
      varbeta0 = vcov["(Intercept)", "(Intercept)"]
      varbeta1 = vcov["tftime", "tftime"]
      covbeta0beta1 = vcov["(Intercept)", "tftime"]
      t <- seq(43, 387, 1)
      plott <- t/43
      t1 = 43
      t2 = 202
      if(p == 0){
        predF <- beta0 + beta1*log(plott)
        a1 = (t2 - t1) 
        a2 = (t2*log(t2) - t1*log(t1) - t2 + t1 - log(t1)*(t2 - t1))
        
      }else if (p == (-1)){
        predF <- beta0 + beta1*plott^p
        a1 = (t2 - t1) 
        a2 = t1*log(t2/t1)
      }else {
        predF <- beta0 + beta1*plott^p
        a1 = (t2 - t1) 
        a2 = 1/t1^p*(t2^(p+1)/(p+1) - t1^(p+1)/(p+1))
      }
      
      AUCest = beta0*a1+ beta1*a2
      AUCestnormalized = AUCest/(t2-t1)
      AUCnormalizedVar = (varbeta0*a1^2 + varbeta1*a2^2 + 2*covbeta0beta1*a1*a2)/(t2-t1)^2
      est = c(AUCestnormalized, AUCnormalizedVar, AUCestnormalized, AUCnormalizedVar)
      names(est) <- c("AUC_Female", "Var_Female", "AUC_Male", "Var_Male")
      npredF <- cbind(t, predF, predF)
      colnames(npredF) <- c("time","Female", "Male")
      return(list(est = est, pred = npredF))
    }
    
  }else{
    est = c(NA, NA, NA, NA)
    names(est) <- c("AUC_Female", "Var_Female", "AUC_Male", "Var_Male")
    return(list(est = est, pred = NA))
  }
  
}


AntibodyLevelsummaryConditional<- function(fmsummary, p, t){
  if(!is.na(p)){
    coefest <- as.numeric(fmsummary$coefficients[, "Estimate"])
    names(coefest) <- rownames(fmsummary$coefficient)
    vcov <- fmsummary$vcov
    names(coefest) <- rownames(fmsummary$coefficient)
    rownames(vcov) <- rownames(fmsummary$coefficient)
    colnames(vcov) <- rownames(fmsummary$coefficient)
    
    if(!is.na(coefest["Sex"])){
      beta0 = coefest["(Intercept)"]
      beta1 = coefest["tftime"]
      beta2 = coefest["Sex"]
      covbeta = vcov[c("(Intercept)", "tftime", "Sex"), c("(Intercept)", "tftime", "Sex")]
      
      if(p == 0){
        a2 = log(t/43)
      }else{
        a2 = (t/43)^p
      }
      
      AntibodyLevelEstFemale = beta0+ beta1*a2
      AntibodyLevelVarFemale = (covbeta[1, 1] + covbeta[2, 2]*a2^2 + 2*covbeta[1, 2]*a2)
      
      AntibodyLevelEstMale = beta0+beta2+ beta1*a2
      AntibodyLevelVarMale = covbeta[1, 1] + covbeta[3, 3] + 2*covbeta[1, 3]+ covbeta[2, 2]*a2^2 + 2*(covbeta[1, 2] + covbeta[2, 3])*a2
      est <- c(AntibodyLevelEstFemale, AntibodyLevelVarFemale, AntibodyLevelEstMale, AntibodyLevelVarMale)
      names(est) <- c("Est_Female", "Var_Female", "Est_Male", "Var_Male")
      
      return(est)
    }else{
      beta0 = coefest["(Intercept)"]
      beta1 = coefest["tftime"]
      varbeta0 = vcov["(Intercept)", "(Intercept)"]
      varbeta1 = vcov["tftime", "tftime"]
      covbeta0beta1 = vcov["(Intercept)", "tftime"]
      if(p == 0){
        a2 = log(t/43)
      }else{
        a2 = (t/43)^p
      }
      
      AntibodyLevelEst = beta0+ beta1*a2
      AntibodyLevelVar = varbeta0 + varbeta1*a2^2 + 2*covbeta0beta1*a2
      est = c(AntibodyLevelEst, AntibodyLevelVar, AntibodyLevelEst, AntibodyLevelVar)
      names(est) <- c("Est_Female", "Var_Female", "Est_Male", "Var_Male")
      
      return(est)
    }
  }else{
    est = c(NA, NA, NA, NA)
    names(est) <- c("Est_Female", "Var_Female", "Est_Male", "Var_Male")
    
    return(est)
  }
  
  
}
