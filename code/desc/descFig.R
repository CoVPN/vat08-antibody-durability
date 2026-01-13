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
source(file.path(codeDir, "desc/descFigUtils.R"))

nAbdataLong <- read.csv(file.path(datDir, "vat08_combined_data_processed_longitudinal_nAb.csv"))
bAbdataLong <- read.csv(file.path(datDir, "vat08_combined_data_processed_longitudinal_bAb.csv"))


markers <- c("bindSpike", "bindSpike_omicron", "pseudoneutid50",  "pseudoneutid50_BA.4.5")
nAbID50titers <- c("pseudoneutid50", "pseudoneutid50_BA.4.5")


#plot the trajectory for the longitudinal analysis cohorts before censoring
for(markeri in markers){
  LLOQmarker <- LLOQf (markeri)
  if(markeri %in% nAbID50titers){
    dat_proc0 <- nAbdataLong
  }else{
    dat_proc0 <- bAbdataLong
  }
  dat_proc0$marker <- dat_proc0[, markeri]
  
  dat_proc0$Trt2 <- ifelse(dat_proc0$Trt == 0, "Placebo","Vaccine")
  dat_proc0$Trialstage2 <- ifelse(dat_proc0$Trialstage == 1, "Stage 1","Stage 2")
  dat_proc0$Trt2 <- factor(dat_proc0$Trt2, levels = c("Vaccine", "Placebo"))
  dat_proc0$Trialstage2 <- factor(dat_proc0$Trialstage2, levels = c("Stage 1", "Stage 2"))
  
  dat_proc1 <- filter(dat_proc0, Bserostatus == 0 )
 
 if(markeri %in% nAbID50titers){
    ylim <- c(1.25, 5.5)
  }else{
    ylim <- c(2.3, 6.1)
  }
 
  pnaive <-  ggplot(data = dat_proc1) +
    geom_line(aes(x = time, y = marker, group = Ptid), color = "black", linewidth = 0.2) +
    geom_point(aes(x = time, y = marker, group = Ptid), size = 1.8, shape =1, stroke = 0.7, alpha = 1)+
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
  
  
  dat_proc1 <- filter(dat_proc0, Bserostatus == 1 )
  
  pnonnaive <-ggplot(data = dat_proc1) +
    geom_line(aes(x = time, y = marker, group = Ptid), color = "black", linewidth = 0.2) +
    geom_point(aes(x = time, y = marker, group = Ptid), size = 1.8, shape =1, stroke = 0.7, alpha = 1)+
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


