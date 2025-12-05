rm(list=ls(all=TRUE))
here::i_am("README.md")
repoDir <- here::here()
datDir <- file.path(repoDir,"data")
codeDir <- file.path(repoDir, "code")
outputDir <- file.path(repoDir, "output")
figDir <- file.path(repoDir, "figures")

library(tidyverse)
library(plyr)
source(file.path(codeDir, "common.R"))
source(file.path(codeDir, "desc/descFig2Utils.R"))

longitudinalSubcohortPtid <- read.csv(file.path(datDir, "longitudinalSubcohortPtid.csv"))
dat_proc0 <- read_csv(file.path(datDir, paste0("vat08_combined_data_processed_withmissingIndicator.csv")))
dat_proc0$Trt2 <- ifelse(dat_proc0$Trt == 0, "Placebo","Vaccine")
dat_proc0$Trialstage2 <- ifelse(dat_proc0$Trialstage == 1, "Stage 1","Stage 2")
dat_proc0$Trt2 <- factor(dat_proc0$Trt2, levels = c("Vaccine", "Placebo"))
dat_proc0$Trialstage2 <- factor(dat_proc0$Trialstage2, levels = c("Stage 1", "Stage 2"))
nAbdataLong <- read.csv(file.path(datDir, "vat08_combined_data_processed_longitudinal_nAb.csv"))
bAbdataLong <- read.csv(file.path(datDir, "vat08_combined_data_processed_longitudinal_bAb.csv"))


markers <- c("bindSpike", "bindSpike_beta", "bindSpike_alpha", "bindSpike_gamma", "bindSpike_delta1", "bindSpike_delta2", "bindSpike_delta3", "bindSpike_omicron",
             "pseudoneutid50", "pseudoneutid50_B.1.351", "pseudoneutid50_BA.1", "pseudoneutid50_BA.2", "pseudoneutid50_BA.4.5")
nAbID50titers <- c("pseudoneutid50", "pseudoneutid50_B.1.351", "pseudoneutid50_BA.1", "pseudoneutid50_BA.2", "pseudoneutid50_BA.4.5")
#restrict plotting to participants in the longitudinal sub-cohort
dat_proc0 <- filter(dat_proc0, Ptid %in% longitudinalSubcohortPtid$x)
for(markeri in markers){
  LLOQmarker <- LLOQf (markeri)
  dat_proc0$marker <- as.vector(dat_proc0[, markeri])[[1]]
  
  dat_proc1 <- filter(dat_proc0, Bserostatus == 0)
  dat_proc1$missing <- ifelse(dat_proc1$time2<=22, "Pre-Day 43 levels",
                              ifelse(dat_proc1$missingDueToSymInfection_bAb==1, 
                                     "After primary endpoint diagnosis date", 
                                     ifelse(dat_proc1$missingDueToAsymInfection_bAb==1, 
                                            "After asymptomatic infection criteria met", "Asymptomatic infection criteria not met")))
  colorCategory <- c("Pre-Day 43 levels", "After primary endpoint diagnosis date",
                     "After asymptomatic infection criteria met", "Asymptomatic infection criteria not met")
  dat_proc1$missing <- factor(dat_proc1$missing, levels = colorCategory)
  
  if(markeri %in% nAbID50titers){
    ylim <- c(1, 6)
  }else{
    ylim <- c(2, 7)
  }
  pnaive <-  ggplot(data = dat_proc1) +
    geom_line(aes(x = time, y = marker, group = Ptid, color = missing, alpha = missing, linewidth = missing)) +
    geom_point(aes(x = time, y = marker, group = Ptid, color = missing, size = missing, shape = missing), 
               stroke = 0.7, alpha = 1)+
    ylab(plot_labf(markeri)) +
    xlab("Days since Enrollment")+
    geom_hline(yintercept = LLOQmarker, linetype = 'dotted')+
    annotate(geom = "text", x = 380, y = LLOQmarker-0.15, 
             label = ifelse(markeri%in% c("pseudoneutid50", "pseudoneutid50_B.1.351", "pseudoneutid50_BA.1", 
                                          "pseudoneutid50_BA.2", "pseudoneutid50_BA.4.5"), "LOD", "LLOQ"), size = 3) +
    facet_grid(cols = vars(Trt2), rows = vars(Trialstage2)) +
    scale_y_continuous(limits = ylim, 
                       breaks = seq(floor(ylim[1]), ceiling(ylim[2]), 1), 
                       labels = scientific_10(10^seq(floor(ylim[1]), ceiling(ylim[2]), 1))) +
    scale_x_continuous(breaks = c(0, 22, 43, 78, 134, 202, 292, 387 ), limits = c(-15, 410), 
                       expand = expansion(0, 0), 
                       labels = c("0","22","43", "78", "134", "202", "292", "387" ),
                       minor_breaks = NULL) +
    scale_color_manual(breaks = colorCategory, 
                       values = c( "gray40", "#1f77b4", "#ff7f0e", "darkolivegreen")) +
    scale_size_manual(breaks = colorCategory, values = c(1.3, 1.8, 1.3, 1.3)) +
    scale_alpha_manual(breaks = colorCategory, values = c(1, 1, 0.7, 0.8)) +
    scale_shape_manual(breaks = colorCategory, values = c(16, 2, 4, 1)) +
    scale_linewidth_manual(breaks = colorCategory, values = c(0.1, 0.5, 0.3, 0.1)) +
    guides(color = guide_legend(nrow = 1))+
    theme_bw()+
    ggtitle("Naïve")+
    theme(legend.title = element_blank(),
          plot.margin = unit(c(1, 1, 1, 1), "lines"),
          legend.direction = "vertical",
          legend.key.size=unit(0.5, "cm"),
          legend.text = element_text (size = 12),
          legend.position = "top", legend.justification = c(1, 1),
          legend.spacing = unit(0.13, "cm"),
          legend.background = element_rect(fill = NA),
          title = element_text(vjust = -5, size = 20),
          strip.text.y =  element_text (size = 16),
          strip.text.x =  element_text (size = 16),
          axis.title.x =  element_text (size = 18, vjust = -0.8),
          axis.title.y =  element_text (size = 18),
          axis.text.x =  element_text (size = 14,  vjust = 0.5),
          axis.text.y =  element_text (size = 14))
  

  dat_proc2 <- filter(dat_proc0, Bserostatus == 1)
  dat_proc2$missing <- ifelse(dat_proc2$time2<=22, "Pre-Day 43 levels",
                              ifelse(dat_proc2$missingDueToSymInfection_bAb==1, 
                                     "After primary endpoint diagnosis date", 
                                     ifelse(dat_proc2$missingDueToAsymInfection_bAb==1, 
                                            "After asymptomatic infection criteria met", "Asymptomatic infection criteria not met")))
  
  
  colorCategory <- c("Pre-Day 43 levels", "After primary endpoint diagnosis date",
                     "After asymptomatic infection criteria met", "Asymptomatic infection criteria not met")
 
  dat_proc2$missing <- factor(dat_proc2$missing, levels = colorCategory)
  
 
  pnonnaive <-ggplot(data = dat_proc2) +
    geom_line(aes(x = time, y = marker, group = Ptid, color = missing, alpha = missing, linewidth = missing)) +
    geom_point(aes(x = time, y = marker, group = Ptid, color = missing, size = missing, shape = missing), 
               stroke = 0.7, alpha = 1)+
    ylab(plot_labf(markeri)) +
    xlab("Days since Enrollment")+
    geom_hline(yintercept = LLOQmarker, linetype = 'dotted')+
    annotate(geom = "text", x = 380, y = LLOQmarker-0.15, 
             label = ifelse(markeri%in% c("pseudoneutid50", "pseudoneutid50_B.1.351", "pseudoneutid50_BA.1", 
                                          "pseudoneutid50_BA.2", "pseudoneutid50_BA.4.5"), "LOD", "LLOQ"), size = 3) +
    facet_grid(cols = vars(Trt2), rows = vars(Trialstage2)) +
    scale_y_continuous(limits = ylim, 
                       breaks = seq(floor(ylim[1]), ceiling(ylim[2]), 1), 
                       labels = scientific_10(10^seq(floor(ylim[1]), ceiling(ylim[2]), 1))) +
    scale_x_continuous(breaks = c(0, 22, 43, 78, 134, 202, 292, 387 ), limits = c(-15, 410), 
                       expand = expansion(0, 0), 
                       labels = c("0","22","43", "78", "134", "202", "292", "387" ),
                       minor_breaks = NULL) +
    scale_color_manual(breaks = colorCategory, 
                       values = c( "gray40", "#1f77b4", "#ff7f0e", "darkolivegreen")) +
    scale_size_manual(breaks = colorCategory, values = c(1.3, 1.8, 1.3, 1.3)) +
    scale_alpha_manual(breaks = colorCategory, values = c(1, 1, 0.7, 0.8)) +
    scale_shape_manual(breaks = colorCategory, values = c(16, 2, 4, 1)) +
    scale_linewidth_manual(breaks = colorCategory, values = c(0.1, 0.5, 0.3, 0.1)) +
    guides(color = guide_legend(nrow = 1))+
    theme_bw()+
    ggtitle("Non-naïve")+
    theme(legend.title = element_blank(),
          plot.margin = unit(c(1, 1, 1, 1), "lines"),
          legend.direction = "vertical",
          legend.key.size=unit(0.5, "cm"),
          legend.text = element_text (size = 12),
          legend.spacing = unit(0.13, "cm"),
          legend.background = element_rect(fill = NA),
          legend.position = "top", legend.justification = c(1, 1),
          title = element_text(size = 20),
          strip.text.y =  element_text (size = 16),
          strip.text.x =  element_text (size = 16),
          axis.title.x =  element_text (size = 18, vjust = -0.8),
          axis.title.y =  element_text (size = 18),
          axis.text.x =  element_text (size = 14,  vjust = 0.5),
          axis.text.y =  element_text (size = 14))
  ggsave(filename = paste0("naive", markeri,"_trajectoryPlot_fulldata.pdf"), 
         plot = pnaive, path = figDir, width = 13, height = 8, units = "in") 
  ggsave(filename = paste0("nonnaive", markeri,"_trajectoryPlot_fulldata.pdf"), 
         plot = pnonnaive, path = figDir, width = 13, height = 8, units = "in")
}


#plot the trajectory for the longitudinal analysis cohorts before censoring
for(markeri in markers){
  LLOQmarker <- LLOQf (markeri)
  dat_proc0$marker <- as.vector(dat_proc0[, markeri])[[1]]
  
  dat_proc1 <- filter(dat_proc0, Bserostatus == 0 )
  dat_proc1$missing <- ifelse(dat_proc1$time2<=22, "Pre-Day 43 levels",
                              ifelse(dat_proc1$missingDueToSymInfection_bAb==1, 
                                     "At or after primary endpoint diagnosis date", 
                                     ifelse(dat_proc1$missingDueToAsymInfection_bAb==1, 
                                            "After asymptomatic infection criteria met", "Asymptomatic infection criteria not met")))
  
  colorCategory <- c("Pre-Day 43 levels", "At or after primary endpoint diagnosis date",
                     "After asymptomatic infection criteria met", "Asymptomatic infection criteria not met")
  dat_proc1$missing <- factor(dat_proc1$missing, levels = colorCategory)
  
  if(markeri %in% nAbID50titers){
    dat_proc1 <- filter(dat_proc1, Ptid %in% nAbdataLong$Ptid & missing == "Asymptomatic infection criteria not met")
  }else{
    dat_proc1 <- filter(dat_proc1, Ptid %in% bAbdataLong$Ptid & missing == "Asymptomatic infection criteria not met")
  }
  
  
 if(markeri %in% nAbID50titers){
    ylim <- c(1.25, 5.5)
  }else{
    ylim <- c(2.3, 6.1)
  }
 
  pnaive <-  ggplot(data = dat_proc1) +
    geom_line(aes(x = time, y = marker, group = Ptid, color = missing, alpha = missing, linewidth = missing)) +
    geom_point(aes(x = time, y = marker, group = Ptid, color = missing, size = missing, shape = missing), 
               stroke = 0.7, alpha = 1)+
    ylab(plot_labf(markeri)) +
    xlab("Days since Enrollment")+
    geom_hline(yintercept = LLOQmarker, linetype = 'dotted')+
    annotate(geom = "text", x = 380, y = LLOQmarker-0.1, 
             label = ifelse(markeri%in% c("pseudoneutid50", "pseudoneutid50_B.1.351", "pseudoneutid50_BA.1", 
                                          "pseudoneutid50_BA.2", "pseudoneutid50_BA.4.5"), "LOD", "LLOQ"), size = 4) +
    facet_grid(cols = vars(Trt2), rows = vars(Trialstage2)) +
    scale_y_continuous(limits = ylim, 
                       breaks = seq(floor(ylim[1]), ceiling(ylim[2]), 1), 
                       labels = scientific_10(10^seq(floor(ylim[1]), ceiling(ylim[2]), 1)),
                       minor_breaks = NULL) +
    scale_x_continuous(breaks = c(0, 22, 43, 78, 134, 202, 292, 387 ), limits = c(30, 410), 
                       expand = expansion(0, 0), 
                       labels = c("0","22","43", "78", "134", "202", "292", "387" ),
                       minor_breaks = NULL) +
    scale_color_manual(breaks = colorCategory, 
                       values = c( "gray40", "#1f77b4", "#ff7f0e", "black")) +
    scale_size_manual(breaks = colorCategory, values = c(1.3, 1.8, 1.3, 1.8)) +
    scale_alpha_manual(breaks = colorCategory, values = c(1, 1, 0.7, 1)) +
    scale_shape_manual(breaks = colorCategory, values = c(16, 2, 4, 1)) +
    scale_linewidth_manual(breaks = colorCategory, values = c(0.1, 0.5, 0.3, 0.2)) +
    guides(color = "none", shape = "none", size = "none", "linewidth" = "none", "alpha" = "none")+
  
    theme_bw()+
    ggtitle("Naïve")+
    theme(legend.title = element_blank(),
          plot.margin = unit(c(1, 1, 1, 1), "lines"),
          legend.direction = "vertical",
          legend.key.size=unit(0.5, "cm"),
          legend.text = element_text (size = 12),
          legend.position = "top", legend.justification = c(1, 1),
          legend.spacing = unit(0.13, "cm"),
          legend.background = element_rect(fill = NA),
          title = element_text(vjust = -5, size = 20),
          strip.text.y =  element_text (size = 18),
          strip.text.x =  element_text (size = 18),
          axis.title.x =  element_text (size = 20, vjust = -0.8),
          axis.title.y =  element_text (size = 20),
          axis.text.x =  element_text (size = 18,  vjust = 0.5),
          axis.text.y =  element_text (size = 18))
  
  
  dat_proc2 <- filter(dat_proc0, Bserostatus == 1 )
  dat_proc2$missing <- ifelse(dat_proc2$time2<=22, "Pre-Day 43 levels",
                              ifelse(dat_proc2$missingDueToSymInfection_bAb==1, 
                                     "At or after primary endpoint diagnosis date", 
                                     ifelse(dat_proc2$missingDueToAsymInfection_bAb==1, 
                                            "After asymptomatic infection criteria met", "Asymptomatic infection criteria not met")))
  
  
  colorCategory <- c("Pre-Day 43 levels", "At or after primary endpoint diagnosis date",
                     "After asymptomatic infection criteria met", "Asymptomatic infection criteria not met")
  
  dat_proc2$missing <- factor(dat_proc2$missing, levels = colorCategory)
  

  if(markeri %in% nAbID50titers){
    dat_proc2 <- filter(dat_proc2, Ptid %in% nAbdataLong$Ptid & missing == "Asymptomatic infection criteria not met")
  }else{
    dat_proc2 <- filter(dat_proc2, Ptid %in% bAbdataLong$Ptid & missing == "Asymptomatic infection criteria not met")
  }
  
 
  pnonnaive <-ggplot(data = dat_proc2) +
    geom_line(aes(x = time, y = marker, group = Ptid, color = missing, alpha = missing, linewidth = missing)) +
    geom_point(aes(x = time, y = marker, group = Ptid, color = missing, size = missing, shape = missing), 
               stroke = 0.7, alpha = 1)+
    ylab(plot_labf(markeri)) +
    xlab("Days since Enrollment")+
    geom_hline(yintercept = LLOQmarker, linetype = 'dotted')+
    annotate(geom = "text", x = 380, y = LLOQmarker-0.1, 
             label = ifelse(markeri%in% c("pseudoneutid50", "pseudoneutid50_B.1.351", "pseudoneutid50_BA.1", 
                                          "pseudoneutid50_BA.2", "pseudoneutid50_BA.4.5"), "LOD", "LLOQ"), size = 3) +
    facet_grid(cols = vars(Trt2), rows = vars(Trialstage2)) +
    scale_y_continuous(limits = ylim, 
                       breaks = seq(floor(ylim[1]), ceiling(ylim[2]), 1), 
                       labels = scientific_10(10^seq(floor(ylim[1]), ceiling(ylim[2]), 1)),
                       minor_breaks = NULL) +
    scale_x_continuous(breaks = c(0, 22, 43, 78, 134, 202, 292, 387 ), limits = c(30, 410), 
                       expand = expansion(0, 0), 
                       labels = c("0","22","43", "78", "134", "202", "292", "387" ),
                       minor_breaks = NULL) +
    scale_color_manual(breaks = colorCategory, 
                       values = c( "gray40", "#1f77b4", "#ff7f0e", "black")) +
    scale_size_manual(breaks = colorCategory, values = c(1.3, 1.8, 1.3, 1.8)) +
    scale_alpha_manual(breaks = colorCategory, values = c(1, 1, 0.7, 1)) +
    scale_shape_manual(breaks = colorCategory, values = c(16, 2, 4, 1)) +
    scale_linewidth_manual(breaks = colorCategory, values = c(0.1, 0.5, 0.3, 0.2)) +
    guides(color = "none", shape = "none", size = "none", "linewidth" = "none", "alpha" = "none")+
   
    theme_bw()+
    ggtitle("Non-naïve")+
    theme(legend.title = element_blank(),
          plot.margin = unit(c(1, 1, 1, 1), "lines"),
          legend.direction = "vertical",
          legend.key.size=unit(0.5, "cm"),
          legend.text = element_text (size = 12),
          legend.spacing = unit(0.13, "cm"),
          legend.background = element_rect(fill = NA),
          legend.position = "top", legend.justification = c(1, 1),
          title = element_text(vjust = -5, size = 20),
          strip.text.y =  element_text (size = 18),
          strip.text.x =  element_text (size = 18),
          axis.title.x =  element_text (size = 20, vjust = -0.8),
          axis.title.y =  element_text (size = 20),
          axis.text.x =  element_text (size = 18,  vjust = 0.5),
          axis.text.y =  element_text (size = 18))
  ggsave(filename = paste0("naive_", markeri,"_trajectoryPlot_longAnaCohort.pdf"), 
         plot = pnaive, path = figDir, width = 9, height = 8, units = "in") 
  ggsave(filename = paste0("nonnaive_", markeri,"_trajectoryPlot_longAnaCohort.pdf"), 
         plot = pnonnaive, path = figDir, width = 9, height = 8, units = "in")
}

#pairwise scatter plot and correlation among time points, within each relevant group and for a given antigen, 
for(stage in c(1, 2)){
  for(mark in c("pseudoneutid50", "pseudoneutid50_BA.4.5", "bindSpike", "bindSpike_omicron")){
    if(mark %in% c("pseudoneutid50", "pseudoneutid50_BA.4.5")){
     
      data <- pivot_wider(select(filter(nAbdataLong, Trialstage == stage), all_of(c("Ptid", "Bserostatus", "Trt", "Trialstage", "time2", mark))), 
                          names_from = time2, values_from = mark)
      
      range.mark <- range (filter(filter(nAbdataLong, Trialstage == stage), !(Bserostatus == 0 & Trt == 0))[, mark])
    }else{
      data <- pivot_wider(select(filter(bAbdataLong,Trialstage == stage), all_of(c("Ptid", "Bserostatus", "Trt", "Trialstage", "time2", mark))), 
                          names_from = time2, values_from = mark)
      range.mark <- range (filter(filter(bAbdataLong,Trialstage == stage), !(Bserostatus == 0 & Trt == 0))[, mark])
    }
    
    data$bs <- ifelse(data$Bserostatus == 0, "Naïve", "Non-naïve")
    data$Trt2 <- ifelse(data$Trt == 0, "Placebo","Vaccine")
    data$Trialstage2 <- ifelse(data$Trialstage == 1, "Stage 1","Stage 2")
    
    data$Trt2 <- factor(data$Trt2, levels = c("Vaccine", "Placebo"))
    data$Trialstage2 <- factor(data$Trialstage2, levels = c("Stage 1", "Stage 2"))
    
    
    plotdata <- select(data, all_of(c("Trialstage","Trt2", "bs", "43", "78","134", "202", "292", "387")))
    
    plotdata$group <- case_when(plotdata$Trt2 == "Vaccine" & plotdata$bs == "Non-naïve" ~ "Non-naïve Vaccine",
                                plotdata$Trt2 == "Placebo" & plotdata$bs == "Non-naïve" ~ "Non-naïve Placebo",
                                plotdata$Trt2 == "Vaccine" & plotdata$bs == "Naïve" ~ "Naïve Vaccine",
                                plotdata$Trt2 == "Placebo" & plotdata$bs == "Naïve" ~ "Naïve Placebo")
    plotdata$group <- factor(plotdata$group, levels = c("Non-naïve Vaccine", "Non-naïve Placebo", "Naïve Vaccine", "Naïve Placebo"))
    plotdata <- filter(plotdata, group != "Naïve Placebo") 
    labels = c("Day 43", "Day 78", "Day 134", "Day 202", "Day 292", "Day 387")
    
    
    library(GGally)
    
    if ("plyr" %in% .packages()) {
      detach("package:plyr", unload=TRUE) #this cause group_by to fail in dplyr
    }
    library(dplyr)
    
   
    
    if(mark %in% c("pseudoneutid50", "pseudoneutid50_BA.4.5")){
      
      scatterplot_lower <- function(data, mapping, ...) {
        ggplot(data = data, mapping = mapping, ...) +
          geom_point(...) +
          scale_y_continuous(limits = c(1.5, 5.5), breaks = c(2, 3, 4, 5), labels = c(2, 3, 4, 5)) + # Set desired y-axis limits
          scale_x_continuous(limits = c(1.5, 5.5), breaks = c(2, 3, 4, 5), labels = c(2, 3, 4, 5))   # Set desired x-axis limits
          
            }
      
      density_diag <- function(data, mapping, ...) {
        ggplot(data = data, mapping = mapping, ...) +
          geom_density(...) +
          scale_x_continuous(limits = c(1.5, 5.5), breaks = c(2, 3, 4, 5), labels = c(2, 3, 4, 5))
      }
    }else{
      scatterplot_lower <- function(data, mapping, ...) {
        ggplot(data = data, mapping = mapping, ...) +
          geom_point(...) +
          scale_y_continuous(limits = c(2.3, 6.2), breaks = c(2, 3, 4, 5,6), labels = c(2, 3, 4, 5, 6)) + # Set desired y-axis limits
          scale_x_continuous(limits = c(2.3, 6.2), breaks = c(2, 3, 4, 5,6), labels = c(2, 3, 4, 5, 6))   # Set desired x-axis limits
        
        }
      
      density_diag <- function(data, mapping, ...) {
        ggplot(data = data, mapping = mapping, ...) +
          geom_density(...) +
          scale_x_continuous(limits = c(2.3, 6.2), breaks = c(2, 3, 4, 5,6), labels = c(2, 3, 4, 5, 6))
        }
    }
    
    p <- ggpairs(plotdata, columns = 4:9, 
                 mapping = ggplot2::aes(shape = group, color = group, alpha = 0.6), 
                 columnLabels = labels,
                 upper = list(continuous = mycorrelations),
                 diag = list(continuous = density_diag ),
                 lower = list(continuous = scatterplot_lower)
    )
    p <- p + scale_color_manual(breaks = c("Non-naïve Vaccine", "Non-naïve Placebo", "Naïve Vaccine", "Naïve Placebo"),
                                values = c("gray20", "#1f77b4", "#ff7f0e", "darkolivegreen"))+
      scale_fill_manual(breaks = c("Non-naïve Vaccine", "Non-naïve Placebo", "Naïve Vaccine", "Naïve Placebo"),
                        values = c("gray20", "#1f77b4", "#ff7f0e", "darkolivegreen"))+
      ggtitle(labf(mark))+
      theme_bw()+
      theme(
        text = element_text(size = 9),   # reduce font size here
        axis.text = element_text(size = 9),
        strip.text = element_text(size = 9),
        legend.position = "top",                     # Move legend to the top
        legend.box = "horizontal"
      )
    ggsave(filename = paste0("ScatterPlot_among_visits_", mark, "_stage", stage,".pdf"), 
           plot = p, path = figDir, width = 10.5, height = 8, units = "in") 
  }
  
}



