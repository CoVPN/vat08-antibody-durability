rm(list=ls(all=TRUE))
here::i_am("Readme.txt")
repoDir <- here::here()
datDir <- file.path(repoDir,"data")
codeDir <- file.path(repoDir, "code")
outputDir <- file.path(repoDir, "output")
figDir <- file.path(repoDir, "figures")

library(tidyverse)
library(ggplot2)
library(plyr)
library(lme4)
source(file.path(codeDir, "common.R"))
source(file.path(codeDir, "LMMmodeling/utils.R"))

nAbdata <- read.csv(file.path(datDir, "vat08_combined_data_processed_longitudinal_nAb.csv"))
bAbdata <- read.csv(file.path(datDir, "vat08_combined_data_processed_longitudinal_bAb.csv"))

plotFitted = TRUE

markers <- c("bindSpike", "bindSpike_beta", "bindSpike_alpha", "bindSpike_gamma", "bindSpike_delta1", "bindSpike_delta2", "bindSpike_delta3", "bindSpike_omicron",
             "pseudoneutid50", "pseudoneutid50_B.1.351", "pseudoneutid50_BA.1", "pseudoneutid50_BA.2", "pseudoneutid50_BA.4.5")



#Fit the linear mixed models
df1pool <- data.frame()
df2pool <- data.frame()
for(g in 1:dim(groups)[1]){
  rawdatapool <- data.frame()
  
  setup <- groups[g, ]
  stage <- setup$Stage
  Bs <- ifelse(setup$Baseline=="naive", 0, 1)
  antiSonly <- ifelse(setup$Baseline == "nonnaiveD01SposOnly", 1, 0)
  tr <- setup$Treatment
  
  
  if(setup$Baseline == "nonnaive"){
    nAbdatag = filter(nAbdata, Trialstage == stage & Trt == tr & Bserostatus == 1)
    bAbdatag = filter(bAbdata, Trialstage == stage & Trt == tr & Bserostatus == 1)
  }else{
    nAbdatag = filter(nAbdata, Trialstage == stage & Trt == tr & Bserostatus == Bs & D01_S_pos_only_in_non_naive_group == antiSonly)
    bAbdatag = filter(bAbdata, Trialstage == stage & Trt == tr & Bserostatus == Bs & D01_S_pos_only_in_non_naive_group == antiSonly)
  }
  
  if(g %in% c(1, 4)){ #naive groups
    baselineAdj = FALSE
  }else{
    baselineAdj = TRUE
  }
  
  
  for(marker in markers){
    if(marker %in% c("pseudoneutid50", "pseudoneutid50_B.1.351", "pseudoneutid50_BA.1", "pseudoneutid50_BA.2", "pseudoneutid50_BA.4.5")){
      lmmdata <- nAbdatag 
    }else{
      lmmdata <- bAbdatag 
    }
    lmmdata$mark <- lmmdata[, marker]
    lmmdata$Bmark <- lmmdata[, paste0("B", marker)]
  
    
    #find mean follow-up visit at each time
    
    yall <- lmmdata$mark
    df1 <- data.frame()
    for(timex in c(43, 78, 134, 202, 292, 387)){
      y <- yall[lmmdata$time2 == timex]
      ybar <- mean(y, na.rm = TRUE)
      ll <- quantile(y, probs = 0.025, na.rm = TRUE)
      ul <- quantile(y, probs = 0.975, na.rm = TRUE)
      nobs <- sum(!is.na(y))
      ans <- c(timex, ybar, ll, ul, nobs, min(y, na.rm = TRUE), max(y, na.rm = TRUE))
      df1 <- rbind(df1, ans)
    }
    names(df1) <- c("time","mean", "low", "high", "nobs", "min", "max")
    if(tr == 1){
      df1$plottime <- df1$time 
    }else{
      df1$plottime <- df1$time +7
    }
    
    df1$markerlabel <-labf(marker)
    df1$marker <- marker
    df1$LLOQ <- LLOQf(marker)
    df1$treatment <- ifelse(tr == 1, "Vaccine", "Placebo")
    df1$stage <- stage
    df1$Bgroup <- ifelse(setup$Baseline=="naive", "Naive", 
                         ifelse(setup$Baseline=="nonnaiveD01SposOnly", "Non-naive Group A",
                                ifelse(setup$Baseline == "nonnaive", "Non-naive", "Non-naive Group B")))
    stage <- setup$Stage
    Bs <- ifelse(setup$Baseline=="naive", 0, 1)
    antiSonly <- ifelse(setup$Baseline == "nonnaiveD01SposOnly", 1, 0)
    df1$Bgroupfile <- setup$Baseline
    df1pool <- rbind(df1pool, df1)
    #fractional polynomial
    if(plotFitted){
      fit <- BFpolynomialFit(lmmdata, baselineAdj)
      fitsummary <- AUCsummaryMarginal(fit$fmsummary, lmmdata, fit$p)
      df2 <- data.frame(fitsummary$pred)
    }
    df2$markerlabel <-labf(marker)
    df2$marker <- marker
    df2$LLOQ <- LLOQf(marker)
    df2$treatment <- ifelse(tr == 1, "Vaccine", "Placebo")
    df2$stage <- stage
    df2$Bgroup <- ifelse(setup$Baseline=="naive", "Naive", 
                         ifelse(setup$Baseline=="nonnaiveD01SposOnly", "Non-naive Group A",
                                ifelse(setup$Baseline == "nonnaive", "Non-naive", "Non-naive Group B")))
    df2$Bgroupfile <- setup$Baseline
    df2pool <- rbind(df2pool, df2)
    
  }}


df1pool$markerlabel <- factor(df1pool$markerlabel, levels = c("bAb-IgG Spike Reference", "bAb-IgG Spike Omicron-B.1.1.529", "bAb-IgG Spike Alpha", 
                                                              "bAb-IgG Spike Beta", "bAb-IgG Spike Delta-AY.4", "bAb-IgG Spike Delta-B.1.617.2", 
                                                              "bAb-IgG Spike Delta-B.1.617.2/AY.4", "bAb-IgG Spike Gamma",
                                                              "nAb-ID50 Reference", "nAb-ID50 Beta", "nAb-ID50 Omicron-BA.1", 
                                                              "nAb-ID50 Omicron-BA.2", "nAb-ID50 Omicron-BA.4/BA.5"))


df2pool$markerlabel <- factor(df2pool$markerlabel, levels = c("bAb-IgG Spike Reference", "bAb-IgG Spike Omicron-B.1.1.529", "bAb-IgG Spike Alpha", 
                                                              "bAb-IgG Spike Beta", "bAb-IgG Spike Delta-AY.4", "bAb-IgG Spike Delta-B.1.617.2", 
                                                              "bAb-IgG Spike Delta-B.1.617.2/AY.4", "bAb-IgG Spike Gamma",
                                                              "nAb-ID50 Reference", "nAb-ID50 Beta", "nAb-ID50 Omicron-BA.1", 
                                                              "nAb-ID50 Omicron-BA.2", "nAb-ID50 Omicron-BA.4/BA.5"))

mainmarkers <- c("bAb-IgG Spike Reference","nAb-ID50 Reference", "bAb-IgG Spike Omicron-B.1.1.529","nAb-ID50 Omicron-BA.4/BA.5" )

library(scales)
scientific_10 <- function(x) {
  y <- NULL
  for(i in 1:length(x)){
    y[i] <- parse(text=gsub("1e\\+*", " 10^", scales::scientific_format()(x[i])))
  }
  return(y)
}

#plot the trajectories for each group and marker
for(markeri in markers){
  if(markeri %in% c("pseudoneutid50", "pseudoneutid50_B.1.351", "pseudoneutid50_BA.1", "pseudoneutid50_BA.2", "pseudoneutid50_BA.4.5")){
    ylim <- c(1.25, 5.5)
  }else{
    ylim <- c(2.3, 6.1)
  }
  LLOQmarker <- LLOQf (markeri)
  for(s in c(1, 2)){
    for(bs in c("Naive", "Non-naive","Non-naive Group A","Non-naive Group B" )){
      df1pooli <- filter(df1pool, stage == s & Bgroup == bs & marker == markeri)
      df2pooli <- filter(df2pool, stage == s & Bgroup == bs & marker == markeri)
      if(dim(df1pooli)[1] > 0){
        if(plotFitted){
          vlabel <- ifelse(s == 1, "Non-naïve MV", "Non-naïve BV")
          plabel <- ifelse(s == 1, "Stage 1 Non-naïve Placebo", "Stage 2 Non-naïve Placebo")
          p1 <- ggplot() +
            geom_point(data = df1pooli, aes(x = plottime, y = mean, color = treatment), size = 5) +
            geom_errorbar(data = df1pooli, aes(x = plottime, ymin =low, ymax = high, color = treatment), linewidth = 1.7, width = 10) +
            geom_line(aes(x = time, y = Female, color = treatment), data = df2pooli, linewidth = 1.5) +
            scale_y_continuous(limits = ylim, 
                               breaks = seq(floor(ylim[1]), ceiling(ylim[2]), 1), 
                               labels = scientific_10(10^seq(floor(ylim[1]), ceiling(ylim[2]), 1))) +
            #facet_wrap(vars(markerlabel)) + 
            scale_x_continuous(breaks = c(43, 78, 134, 202, 292, 387 ), limits = c(0, 420), 
                               expand = expansion(0, 0), 
                               labels = c("43", "78", "134", "202", "292", "387" )) +
            geom_hline(yintercept = LLOQmarker, linetype = 'dotted')+
            annotate(geom = "text", x = 18, y = LLOQmarker+0.1, label = ifelse(markeri%in% c("pseudoneutid50", "pseudoneutid50_B.1.351", "pseudoneutid50_BA.1", 
                                                                                             "pseudoneutid50_BA.2", "pseudoneutid50_BA.4.5"), "LOD", "LLOQ"), size = 5) + 
            theme_bw()+
            ylab(plot_labf(markeri)) +
            xlab("Days since Enrollment")+
            scale_color_manual(breaks = c("Vaccine", "Placebo"), 
                               values = c( "darkorange", "#619CFF"), labels = c(vlabel, plabel)) +
            #ggtitle(paste0("Stage ", s, " ",gsub("i","ï",bs)))+
            theme(legend.title = element_blank(),
                  plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "lines"),
                  legend.direction = "vertical",
                  legend.key.width=unit(1,"cm"),
                  legend.text = element_text (size = 22),
                  legend.position = c(0.6, 0.93),
                  legend.background = element_rect(fill = NA),
                  title = element_text(size = 23),
                  strip.text.y =  element_text (size = 9),
                  strip.text.x =  element_text (size = 9),
                  axis.title.x =  element_text (size = 24, vjust = -0.8),
                  axis.title.y =  element_text (size = 24),
                  axis.text.x =  element_text (size = 23,  vjust = 0.5),
                  axis.text.y =  element_text (size = 23))
        }else{
          p1 <- ggplot() +
            geom_point(data = df1pooli, aes(x = plottime, y = mean, color = treatment), size = 5) +
            geom_errorbar(data = df1pooli, aes(x = plottime, ymin =low, ymax = high, color = treatment), linewidth = 1.7, width = 10) +
            scale_y_continuous(limits = ylim, 
                               breaks = seq(floor(ylim[1]), ceiling(ylim[2]), 1), 
                               labels = scientific_10(10^seq(floor(ylim[1]), ceiling(ylim[2]), 1))) +
            #facet_wrap(vars(markerlabel)) + 
            scale_x_continuous(breaks = c(43, 78, 134, 202, 292, 387 ), limits = c(0, 420), 
                               expand = expansion(0, 0), 
                               labels = c("43", "78", "134", "202", "292", "387" )) +
            geom_hline(yintercept = LLOQmarker, linetype = 'dotted')+
            annotate(geom = "text", x = 18, y = LLOQmarker+0.1, label = ifelse(markeri%in% c("pseudoneutid50", "pseudoneutid50_B.1.351", "pseudoneutid50_BA.1", 
                                                                                             "pseudoneutid50_BA.2", "pseudoneutid50_BA.4.5"), "LOD", "LLOQ"), size = 5) + 
            theme_bw()+
            ylab(plot_labf(markeri)) +
            xlab("Days since Enrollment")+
            scale_color_manual(breaks = c("Vaccine", "Placebo"), 
                               values = c( "darkorange", "#619CFF")) +
            #ggtitle(paste0("Stage ", s, " ",gsub("i","ï",bs)))+
            theme(legend.title = element_blank(),
                  plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "lines"),
                  legend.direction = "vertical",
                  legend.key.width=unit(1,"cm"),
                  legend.text = element_text (size = 22),
                  legend.position = c(0.6, 0.9),
                  legend.background = element_rect(fill = NA),
                  title = element_text(size = 23),
                  strip.text.y =  element_text (size = 9),
                  strip.text.x =  element_text (size = 9),
                  axis.title.x =  element_text (size = 23, vjust = -0.8),
                  axis.title.y =  element_text (size = 23),
                  axis.text.x =  element_text (size = 22,  vjust = 0.5),
                  axis.text.y =  element_text (size = 22))
        }
      
        ggsave(filename = paste0("stage",s,"_", unique(df1pooli$Bgroupfile), unique(df1pooli$marker),"_fitted.pdf"), 
               plot = p1, path = figDir, width = 8, height = 8, units = "in") 
        
        
      }
      
    }
  }   
  
}

#Combine stage 1 naive and stage 2 naive

for(markeri in markers){
  if(markeri %in% c("pseudoneutid50", "pseudoneutid50_B.1.351", "pseudoneutid50_BA.1", "pseudoneutid50_BA.2", "pseudoneutid50_BA.4.5")){
    ylim <- c(1.25, 5.5)
  }else{
    ylim <- c(2.5, 6.1)
  }
  LLOQmarker <- LLOQf (markeri)
  df1pooli <- filter(df1pool, Bgroup == "Naive" & marker == markeri & treatment == "Vaccine")
  df2pooli <- filter(df2pool, Bgroup == "Naive" & marker == markeri & treatment == "Vaccine")
  df1pooli$stage <- as.character(df1pooli$stage)
  df1pooli$plottime[df1pooli$stage == "2"] <- df1pooli$plottime[df1pooli$stage == "2"] +7
 
  df2pooli$stage <- as.character(df2pooli$stage)
  df2pooli$plottime[df2pooli$stage == "2"] <- df2pooli$plottime[df2pooli$stage == "2"] +7
  
  if(dim(df1pooli)[1] > 0){
    p1 <- ggplot() +
      geom_point(data = df1pooli, aes(x = plottime, y = mean, color = stage), size = 5) +
      geom_errorbar(data = df1pooli, aes(x = plottime, ymin =low, ymax = high, color = stage), linewidth = 1.7, width = 10) +
      geom_line(aes(x = time, y = Female, color = stage), data = df2pooli, linewidth = 1.5) +
      scale_y_continuous(limits = ylim, 
                         breaks = seq(floor(ylim[1]), ceiling(ylim[2]), 1), 
                         labels = scientific_10(10^seq(floor(ylim[1]), ceiling(ylim[2]), 1))) +
      #facet_wrap(vars(markerlabel)) + 
      scale_x_continuous(breaks = c(43, 78, 134, 202, 292, 387 ), limits = c(0, 420), 
                         expand = expansion(0, 0), 
                         labels = c("43", "78", "134", "202", "292", "387" )) +
      geom_hline(yintercept = LLOQmarker, linetype = 'dotted')+
      annotate(geom = "text", x = 18, y = LLOQmarker+0.1, label = ifelse(markeri%in% c("pseudoneutid50", "pseudoneutid50_B.1.351", "pseudoneutid50_BA.1", 
                                                                                       "pseudoneutid50_BA.2", "pseudoneutid50_BA.4.5"), "LOD", "LLOQ"), size = 5) + 
      
      theme_bw()+
      ylab(plot_labf(markeri)) +
      xlab("Days since Enrollment")+
      scale_color_manual(breaks = c("1", "2"), 
                         values = c( "darkorchid", "limegreen"), labels = c("Naïve MV", "Naïve BV")) +
    
      theme(legend.title = element_blank(),
            plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "lines"),
            legend.direction = "vertical",
            legend.key.width=unit(1,"cm"),
            legend.text = element_text (size = 22),
            legend.position = c(0.68, 0.93),
            legend.background = element_rect(fill = NA),
            title = element_text(size = 23),
            strip.text.y =  element_text (size = 9),
            strip.text.x =  element_text (size = 9),
            axis.title.x =  element_text (size = 23, vjust = -0.8),
            axis.title.y =  element_text (size = 23),
            axis.text.x =  element_text (size = 22,  vjust = 0.5),
            axis.text.y =  element_text (size = 22))
  ggsave(filename = paste0("naive", unique(df1pooli$marker),"_fitted.pdf"), 
           plot = p1, path = figDir, width = 8, height = 8, units = "in") 
    
    
  }
  
  
}


#Combine stage 1 nonnaive and stage 2 nonnaive

for(markeri in markers){
  if(markeri %in% c("pseudoneutid50", "pseudoneutid50_B.1.351", "pseudoneutid50_BA.1", "pseudoneutid50_BA.2", "pseudoneutid50_BA.4.5")){
    ylim <- c(1.25, 5.5)
  }else{
    ylim <- c(2.7, 6.1)
  }
  LLOQmarker <- LLOQf (markeri)
  df1pooli <- filter(df1pool, Bgroup == "Non-naive" & marker == markeri & treatment == "Vaccine")
  df2pooli <- filter(df2pool, Bgroup == "Non-naive" & marker == markeri & treatment == "Vaccine")
  df1pooli$stage <- as.character(df1pooli$stage)
  df1pooli$plottime[df1pooli$stage == "2"] <- df1pooli$plottime[df1pooli$stage == "2"] +7

  df2pooli$stage <- as.character(df2pooli$stage)
  df2pooli$plottime[df2pooli$stage == "2"] <- df2pooli$plottime[df2pooli$stage == "2"] +7

  if(dim(df1pooli)[1] > 0){
    p1 <- ggplot() +
      geom_point(data = df1pooli, aes(x = plottime, y = mean, color = stage), size = 4) +
      geom_errorbar(data = df1pooli, aes(x = plottime, ymin =low, ymax = high, color = stage), linewidth = 1.7, width = 8) +
      geom_line(aes(x = time, y = Female, color = stage), data = df2pooli, linewidth = 1.5) +
      scale_y_continuous(limits = ylim,
                         breaks = seq(floor(ylim[1]), ceiling(ylim[2]), 1),
                         labels = scientific_10(10^seq(floor(ylim[1]), ceiling(ylim[2]), 1))) +
      #facet_wrap(vars(markerlabel)) +
      scale_x_continuous(breaks = c(43, 78, 134, 202, 292, 387 ), limits = c(0, 420),
                         expand = expansion(0, 0),
                         labels = c("43", "78", "134", "202", "292", "387" )) +
      geom_hline(yintercept = LLOQmarker, linetype = 'dotted')+
      annotate(geom = "text", x = 18, y = LLOQmarker+0.1, label = ifelse(markeri%in% c("pseudoneutid50", "pseudoneutid50_B.1.351", "pseudoneutid50_BA.1",
                                                                                       "pseudoneutid50_BA.2", "pseudoneutid50_BA.4.5"), "LOD", "LLOQ"), size = 5) +
      theme_bw()+
      ylab(plot_labf(markeri)) +
      xlab("Day since Enrollment")+
      scale_color_manual(breaks = c("1", "2"),
                         values = c( "navyblue", "coral"), labels = c("Non-naïve MV", "Non-naïve BV")) +
      theme(legend.title = element_blank(),
            plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "lines"),
            legend.direction = "vertical",
            legend.key.width=unit(1,"cm"),
            legend.text = element_text (size = 22),
            legend.position = c(0.68, 0.93),
            legend.background = element_rect(fill = NA),
            title = element_text(size = 23),
            strip.text.y =  element_text (size = 9),
            strip.text.x =  element_text (size = 9),
            axis.title.x =  element_text (size = 23, vjust = -0.8),
            axis.title.y =  element_text (size = 23),
            axis.text.x =  element_text (size = 22,  vjust = 0.5),
            axis.text.y =  element_text (size = 22))

    ggsave(filename = paste0("nonnaive", unique(df1pooli$marker),"_fitted.pdf"),
           plot = p1, path = figDir, width = 8, height = 8, units = "in")
  }
}


#Combine stage 2 nonnaive GA vaccine and stage 2 nonnaive GB vaccine
for(markeri in markers){
  if(markeri %in% c("pseudoneutid50", "pseudoneutid50_B.1.351", "pseudoneutid50_BA.1", "pseudoneutid50_BA.2", "pseudoneutid50_BA.4.5")){
    ylim <- c(1.25, 5.5)
  }else{
    ylim <- c(2.7, 6.1)
  }
  LLOQmarker <- LLOQf (markeri)
  df1pooli <- filter(df1pool, Bgroup %in% c("Non-naive Group A", "Non-naive Group B") &stage == 2 & marker == markeri & treatment == "Vaccine")
  df2pooli <- filter(df2pool, Bgroup %in% c("Non-naive Group A", "Non-naive Group B")&stage == 2 &marker == markeri & treatment == "Vaccine")
  df1pooli$plottime[df1pooli$Bgroup == "Non-naive Group B"] <- df1pooli$plottime[df1pooli$Bgroup == "Non-naive Group B"] +7
  df2pooli$plottime[df2pooli$Bgroup == "Non-naive Group B"] <- df2pooli$plottime[df2pooli$Bgroup == "Non-naive Group B"] +7
  
  
  if(dim(df1pooli)[1] > 0){
    p1 <- ggplot() +
      geom_point(data = df1pooli, aes(x = plottime, y = mean, color = Bgroup), size = 4) +
      geom_errorbar(data = df1pooli, aes(x = plottime, ymin =low, ymax = high, color = Bgroup), linewidth = 1.7, width = 8) +
      geom_line(aes(x = time, y = Female, color = Bgroup), data = df2pooli, linewidth = 1.5) +
      scale_y_continuous(limits = ylim,
                         breaks = seq(floor(ylim[1]), ceiling(ylim[2]), 1),
                         labels = scientific_10(10^seq(floor(ylim[1]), ceiling(ylim[2]), 1))) +
      #facet_wrap(vars(markerlabel)) +
      scale_x_continuous(breaks = c(43, 78, 134, 202, 292, 387 ), limits = c(0, 420),
                         expand = expansion(0, 0),
                         labels = c("43", "78", "134", "202", "292", "387" )) +
      geom_hline(yintercept = LLOQmarker, linetype = 'dotted')+
      annotate(geom = "text", x = 18, y = LLOQmarker+0.1, label = ifelse(markeri%in% c("pseudoneutid50", "pseudoneutid50_B.1.351", "pseudoneutid50_BA.1",
                                                                                       "pseudoneutid50_BA.2", "pseudoneutid50_BA.4.5"), "LOD", "LLOQ"), size = 5) +
      theme_bw()+
      ylab(plot_labf(markeri)) +
      xlab("Day since Enrollment")+
      scale_color_manual(breaks = c("Non-naive Group A", "Non-naive Group B"),
                         values = c( "navyblue", "coral"), labels = c("Non-naïve BV Group A", "Non-naïve BV Group B")) +
      #ggtitle(paste0("Stage ", s, " ",gsub("i","ï",bs)))+
      theme(legend.title = element_blank(),
            plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "lines"),
            legend.direction = "vertical",
            legend.key.width=unit(1,"cm"),
            legend.text = element_text (size = 22),
            legend.position = c(0.68, 0.93),
            legend.background = element_rect(fill = NA),
            title = element_text(size = 23),
            strip.text.y =  element_text (size = 9),
            strip.text.x =  element_text (size = 9),
            axis.title.x =  element_text (size = 23, vjust = -0.8),
            axis.title.y =  element_text (size = 23),
            axis.text.x =  element_text (size = 22,  vjust = 0.5),
            axis.text.y =  element_text (size = 22))
    
   
    ggsave(filename = paste0("nonnaive_stage2_vac", unique(df1pooli$marker),"_fitted.pdf"),
           plot = p1, path = figDir, width = 8, height = 8, units = "in")
  }
}

#Combine stage 2 nonnaive GA placebo and stage 2 nonnaive GB placebo
for(markeri in markers){
  if(markeri %in% c("pseudoneutid50", "pseudoneutid50_B.1.351", "pseudoneutid50_BA.1", "pseudoneutid50_BA.2", "pseudoneutid50_BA.4.5")){
    ylim <- c(1.25, 5.5)
  }else{
    ylim <- c(2.3, 6.3)
  }
  LLOQmarker <- LLOQf (markeri)
  df1pooli <- filter(df1pool, Bgroup %in% c("Non-naive Group A", "Non-naive Group B")&stage == 2 & marker == markeri & treatment != "Vaccine")
  df2pooli <- filter(df2pool, Bgroup %in% c("Non-naive Group A", "Non-naive Group B")&stage == 2 &marker == markeri & treatment != "Vaccine")
  df1pooli$plottime[df1pooli$Bgroup == "Non-naive Group B"] <- df1pooli$plottime[df1pooli$Bgroup == "Non-naive Group B"] +7
  df2pooli$plottime[df2pooli$Bgroup == "Non-naive Group B"] <- df2pooli$plottime[df2pooli$Bgroup == "Non-naive Group B"] +7
  
  if(dim(df1pooli)[1] > 0){
    p1 <- ggplot() +
      geom_point(data = df1pooli, aes(x = plottime, y = mean, color = Bgroup), size = 4) +
      geom_errorbar(data = df1pooli, aes(x = plottime, ymin =low, ymax = high, color = Bgroup), linewidth = 1.7, width = 8) +
      geom_line(aes(x = time, y = Female, color = Bgroup), data = df2pooli, linewidth = 1.5) +
      scale_y_continuous(limits = ylim,
                         breaks = seq(floor(ylim[1]), ceiling(ylim[2]), 1),
                         labels = scientific_10(10^seq(floor(ylim[1]), ceiling(ylim[2]), 1))) +

      scale_x_continuous(breaks = c(43, 78, 134, 202, 292, 387 ), limits = c(0, 420),
                         expand = expansion(0, 0),
                         labels = c("43", "78", "134", "202", "292", "387" )) +
      geom_hline(yintercept = LLOQmarker, linetype = 'dotted')+
      annotate(geom = "text", x = 18, y = LLOQmarker+0.1, label = ifelse(markeri%in% c("pseudoneutid50", "pseudoneutid50_B.1.351", "pseudoneutid50_BA.1",
                                                                                       "pseudoneutid50_BA.2", "pseudoneutid50_BA.4.5"), "LOD", "LLOQ"), size = 5) +
      theme_bw()+
      ylab(plot_labf(markeri)) +
      xlab("Day since Enrollment")+
      scale_color_manual(breaks = c("Non-naive Group A", "Non-naive Group B"),
                         values = c( "navyblue", "coral"), labels = c("Stage 2 Non-naïve Group A Placebo", "Stage 2 Non-naïve Group B Placebo")) +
      #ggtitle(paste0("Stage ", s, " ",gsub("i","ï",bs)))+
      theme(legend.title = element_blank(),
            plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "lines"),
            legend.direction = "vertical",
            legend.key.width=unit(1,"cm"),
            legend.text = element_text (size = 22),
            legend.position = c(0.5, 0.93),
            legend.background = element_rect(fill = NA),
            title = element_text(size = 23),
            strip.text.y =  element_text (size = 9),
            strip.text.x =  element_text (size = 9),
            axis.title.x =  element_text (size = 23, vjust = -0.8),
            axis.title.y =  element_text (size = 23),
            axis.text.x =  element_text (size = 22,  vjust = 0.5),
            axis.text.y =  element_text (size = 22))
    
   
    ggsave(filename = paste0("nonnaive_stage2_plac", unique(df1pooli$marker),"_fitted.pdf"),
           plot = p1, path = figDir, width = 8, height = 8, units = "in")
  }
}



