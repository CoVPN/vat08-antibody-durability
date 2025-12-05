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
source(file.path(codeDir, "desc/descFig1Utils.R"))
source(file.path(codeDir, "desc/descFig2Utils.R"))


df1 <- read.csv(file.path(datDir, "antibodyDurabilityWide_imputed.csv"))
longdf <- read.csv(file.path(datDir, "antibodyDurabilityLong_imputed.csv"))
fulllongdf <- left_join(longdf, df1, by = "Ptid")
fulllongdf$BserostatusD01_S <- ifelse(fulllongdf$Bserostatus==0, "naive", ifelse(fulllongdf$D01_S_pos_only_in_non_naive_group == 1, "nonnaiveD01SposOnly", "nonnaiveElse"))

#compare nonnaive group A/B D1 antibody levels
plotmarkers <- c("bindSpike", "bindSpike_beta", "bindSpike_alpha", "bindSpike_gamma", "bindSpike_delta1", "bindSpike_delta2", "bindSpike_delta3", "bindSpike_omicron",
                 "pseudoneutid50", "pseudoneutid50_B.1.351", "pseudoneutid50_BA.1", "pseudoneutid50_BA.2", "pseudoneutid50_BA.4.5")

plotmarkers1 <- c("bindSpike", "bindSpike_beta", "bindSpike_alpha", "bindSpike_gamma", "bindSpike_delta1", "bindSpike_delta2", "bindSpike_delta3", "bindSpike_omicron"
)

labf <- function(marker){
  case_when(marker == "bindSpike" ~ "bAb-IgG Spike Reference",
            marker == "bindSpike_beta" ~ "bAb-IgG Spike Beta",
            marker == "bindSpike_alpha" ~ "bAb-IgG Spike Alpha",
            marker == "bindSpike_gamma" ~ "bAb-IgG Spike Gamma",
            marker == "bindSpike_delta1" ~ "bAb-IgG Spike Delta-AY.4.2",
            marker == "bindSpike_delta2" ~ "bAb-IgG Spike Delta-B.1.617.2",
            marker == "bindSpike_delta3" ~ "bAb-IgG Spike Delta-B.1.617.2/AY.4",
            marker == "bindSpike_omicron" ~ "bAb-IgG Spike Omicron",
            marker == "pseudoneutid50" ~ "nAb-ID50 Reference",
            marker == "pseudoneutid50_B.1.351" ~ "nAb-ID50 Beta",
            marker == "pseudoneutid50_BA.1" ~ "nAb-ID50 Omicron-BA.1",
            marker == "pseudoneutid50_BA.2" ~ "nAb-ID50 Omicron-BA.2",
            marker == "pseudoneutid50_BA.4.5" ~ "nAb-ID50 Omicron-BA.4/BA.5") 
}
library(scales)
scientific_10 <- function(x) {
  parse(text=gsub("1e\\+*", " 10^", scales::scientific_format()(x)))
}

for(s in c(1, 2)){
  fulllongdfD1 <- filter(fulllongdf, visitn == "Baseline" & Bserostatus==1 & Trialstage==s)
  set.seed(1)
  fulllongdfD1 <- dplyr::select(fulllongdfD1, all_of(c("Ptid","BserostatusD01_S","Trt", plotmarkers)))
  fulllongdfD1_long <- pivot_longer(fulllongdfD1, cols = plotmarkers, names_to = "marker", values_to = "antibodyLevel")
  fulllongdfD1_long$group <- ifelse(fulllongdfD1_long$BserostatusD01_S == "nonnaiveD01SposOnly", "Group A", "Group B")
  fulllongdfD1_long$nmarker <- labf(fulllongdfD1_long$marker)
  
  fulllongdfD1_long$nmarker <- factor(fulllongdfD1_long$nmarker, levels = labf(plotmarkers))

  nAbdata <- read.csv(file.path(datDir, "vat08_combined_data_processed_longitudinal_nAb.csv"))
  bAbdata <- read.csv(file.path(datDir, "vat08_combined_data_processed_longitudinal_bAb.csv"))
  
  df_nAb <- filter(fulllongdfD1, BserostatusD01_S == "nonnaiveD01SposOnly" & Ptid%in%nAbdata$Ptid )
  df_bAb <- filter(fulllongdfD1, BserostatusD01_S == "nonnaiveD01SposOnly" & Ptid%in%bAbdata$Ptid )
  df_Ab <- filter(fulllongdfD1, BserostatusD01_S == "nonnaiveD01SposOnly" & Ptid%in%c(bAbdata$Ptid, nAbdata$Ptid))
  
  
  p1 <- ggplot(filter(fulllongdfD1_long, marker %in% plotmarkers1 & Ptid%in%bAbdata$Ptid) ) +
    geom_violin(aes(x = group, y = antibodyLevel, color = group)) + 
    geom_boxplot(aes(x = group, y = antibodyLevel, color = group), width = 0.2)+
    geom_jitter(aes(x = group, y = antibodyLevel, color = group), width = 0.2, height = 0, size = 0.5)+
    scale_y_continuous(limits = c(2, 6.2), 
                       breaks = seq(1, 6, 1), 
                       labels = scientific_10(10^seq(1, 6, 1)),
                       minor_breaks = NULL) +
    scale_color_manual(breaks = c("Group A", "Group B"), 
                       values = c("darkorange", "darkblue")) + 
    theme_bw() +
    guides(color = "none")+
    facet_wrap(vars(nmarker)) + 
    ylab("D1 bAb-IgG Spike Concentrations (AU/ml)")+
    xlab("")+
    ggtitle("Stage 2 Non-naïve")+
    theme(
      title = element_text (size = 16),
      strip.text.y =  element_text (size = 12),
      strip.text.x =  element_text (size = 12),
      axis.title.x =  element_text (size = 16),
      axis.title.y =  element_text (size = 16),
      axis.text.x =  element_text (size = 16),
      axis.text.y =  element_text (size = 16))
  
  
  
  p2 <- ggplot(filter(fulllongdfD1_long, !marker %in% plotmarkers1& Ptid%in%nAbdata$Ptid)) +
    geom_violin(aes(x = group, y = antibodyLevel, color = group)) + 
    geom_boxplot(aes(x = group, y = antibodyLevel, color = group), width = 0.2)+
    geom_jitter(aes(x = group, y = antibodyLevel, color = group), width = 0.2, height = 0, size = 0.5)+
    scale_color_manual(breaks = c("Group A", "Group B"), 
                       values = c("darkorange", "darkblue")) + 
    theme_bw()+
    guides(color = "none")+
    facet_wrap(vars(nmarker))+
    ylab("D1 nAb-ID50 titers")+
    scale_y_continuous(limits = c(1, 5), 
                       breaks = seq(1, 6, 1), 
                       labels = scientific_10(10^seq(1, 6, 1)),
                       minor_breaks = NULL) +
    xlab("")+
    ggtitle("Stage 2 Non-naïve")+
    theme(
      title = element_text (size = 16),
      strip.text.y =  element_text (size = 13),
      strip.text.x =  element_text (size = 13),
      axis.title.x =  element_text (size = 16),
      axis.title.y =  element_text (size = 16),
      axis.text.x =  element_text (size = 16),
      axis.text.y =  element_text (size = 16))

  ggsave(filename = paste0("boxplot_stage", s, "D1nonnaiveAandB_bAb_longitudinalcohort.pdf"), plot = p1, path = figDir, width = 9.5, height = 7, units = "in") 
  ggsave(filename = paste0("boxplot_stage", s, "D1nonnaiveAandB_nAb_longitudinalcohort.pdf"), plot = p2, path = figDir, width = 9, height = 7, units = "in") 
}



##########################################################################################
#Pearson correlation plot between pairs of markers before imputation
#Individual trajectory plot
ppdat <- rbind(read.csv(file.path(datDir, "COVID_Sanofi_stage1_20250312.csv")),
               read.csv(file.path(datDir, "COVID_Sanofi_stage2_20250312.csv")))

#exclude participants with missing Bserostatus
data <- filter(ppdat, !is.na(Bserostatus) & Perprotocol == 1)
data$bs <- ifelse(data$Bserostatus == 0, "Naïve", "Non-naïve")
data$Trt2 <- ifelse(data$Trt == 0, "Placebo","Vaccine")
data$Trialstage2 <- ifelse(data$Trialstage == 1, "Stage 1","Stage 2")
data$group <- case_when(data$Trt2 == "Vaccine" & data$bs == "Non-naïve" ~ "Non-naïve Vaccine",
                        data$Trt2 == "Placebo" & data$bs == "Non-naïve" ~ "Non-naïve Placebo",
                        data$Trt2 == "Vaccine" & data$bs == "Naïve" ~ "Naïve Vaccine",
                        data$Trt2 == "Placebo" & data$bs == "Naïve" ~ "Naïve Placebo")


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
data <- dplyr::select(data, all_of(c("Subjectid","bs", "Trt2", "Trialstage2", "group", antibodyMarkers)))

#make data long
long_df <- data %>%
  pivot_longer(
    cols = all_of(antibodyMarkers),
    names_to = "marker_name",
    values_to = "value"
  )%>% filter(., !is.na(value)) %>% 
  mutate( day = sub("^(B|Day[0-9]+).*", "\\1", marker_name),
          marker_type = case_when(
            grepl("bindSpike", marker_name) ~ "bind",
            grepl("pseudoneutid50", marker_name) ~ "pseudoneutid50",
            TRUE ~ NA_character_
          ),
          variant = case_when(
            grepl("bindSpike_", marker_name) ~ sub(".*bindSpike_", "", marker_name),
            grepl("pseudoneutid50_", marker_name) ~ sub(".*pseudoneutid50_", "", marker_name),
            TRUE ~ "original"
          ))
#make data wider by variant
wide_by_variant <- long_df %>%
  select(-marker_name) %>%
  pivot_wider(
    names_from = variant,
    values_from = value
  )
wide_by_variant_IgG <- filter(wide_by_variant, marker_type == "bind") 
wide_by_variant_IgG$group <- factor(wide_by_variant_IgG$group, levels = c("Non-naïve Vaccine", "Non-naïve Placebo", "Naïve Vaccine", "Naïve Placebo"))
wide_by_variant_titer <- filter(wide_by_variant, marker_type == "pseudoneutid50") 
wide_by_variant_titer$group <- factor(wide_by_variant_titer$group, levels = c("Non-naïve Vaccine", "Non-naïve Placebo", "Naïve Vaccine", "Naïve Placebo"))


#compute the Pearson correlations by time point and plot it
xmap <- function(x){
  y <- case_when(x=="B" ~ 1,
                 x=="Day22" ~ 2,
                 x=="Day43" ~ 3,
                 x=="Day78" ~ 4,
                 x=="Day134" ~ 5,
                 x=="Day202" ~ 6,
                 x=="Day292" ~ 7,
                 x=="Day387" ~ 8)
}

concs <- c("original", "alpha", "beta", "gamma", "delta1", "delta2", "delta3", "omicron")
names(concs) <- c("IgG Spike Index", "IgG Spike Alpha", "IgG Spike Beta", "IgG Spike Gamma", "IgG Spike Delta-AY.4.2", "IgG Spike Delta-B.1.617.2",
                  "IgG Spike Delta-B.1.617.2/AY.4", "IgG Spike Omicron")

for(s in c("Stage 1", "Stage 2")){
  corrMatrix_conc <- tibble("Time" = character(), "plotx" = numeric(),"group" = character(), "markx" = character(), "marky" = character(), "corr" = numeric())
  for(timex in c("B", "Day22", "Day43", "Day78", "Day134", "Day202", "Day292", "Day387")){
    wide_by_variant_IgGi <- filter(wide_by_variant_IgG, day == timex & Trialstage2 == s)
    for(groupx in c("Non-naïve Vaccine", "Non-naïve Placebo", "Naïve Vaccine", "Naïve Placebo")){
      df <- filter(wide_by_variant_IgGi, group == groupx)
      cori <- cor(df[, c("original", "alpha", "beta", "gamma", "delta1", "delta2", "delta3", "omicron")], 
                  method = "pearson", use = "complete.obs")
      for(kx in 1:(length(concs)-1)){
        for(ky in (kx+1):length(concs)){
          corrMatrix_conc <- add_row(.data = corrMatrix_conc, 
                                     "Time" = ifelse(timex=="B", "D1", gsub("Day","D", timex)),
                                     "plotx" = xmap(timex),
                                     "group" = groupx, "markx" = names(concs)[concs==concs[kx]],
                                     "marky" = names(concs)[concs==concs[ky]], "corr" = cori[concs[kx],concs[ky]])
        }
      }
    }
    
  }
  print(min(corrMatrix_conc$corr))
  corrMatrix_conc$markx <- factor(corrMatrix_conc$markx, levels = names(concs))
  corrMatrix_conc$marky <- factor(corrMatrix_conc$marky, levels = names(concs))
  
  corrMatrix_conc$markxy <- paste(gsub("IgG Spike ","", corrMatrix_conc$markx), "~", 
                                  gsub("IgG Spike ","", corrMatrix_conc$marky))
  p1 <- ggplot(data = corrMatrix_conc)+
    geom_point(aes(x = plotx, y = corr, color = group), size = 2, shape = 1, stroke = 1.5)+
    facet_wrap(~ markxy, ncol = 4)+
    scale_x_continuous(breaks = seq(1, 8,1),
                       labels = c("D1", "D22", "D43","D78","D134","D202","D292","D387"), 
                       minor_breaks = NULL)+
    scale_y_continuous(breaks = seq(0, 1, 0.25), limits = c(-0.05,1.05))+
    ylab("Pearson Correlation of log10 IgG Spike Concentrations")+
    xlab("")+
    scale_color_manual(breaks = c("Non-naïve Vaccine", "Non-naïve Placebo", "Naïve Vaccine", "Naïve Placebo"),
                       values = c("gray20", "#1f77b4", "#ff7f0e", "darkolivegreen"))+
    theme_bw()+
    theme(legend.title = element_blank(),
          plot.margin = unit(c(1, 1, 1, 1), "lines"),
          legend.direction = "horizontal",
          legend.key.size=unit(0.5, "cm"),
          legend.text = element_text (size = 12),
          legend.spacing = unit(0.13, "cm"),
          legend.background = element_rect(fill = NA),
          legend.position = "top", legend.justification = c(1, 1),
          title = element_text(vjust = -5, size = 20),
          strip.text.y =  element_text (size = 10),
          strip.text.x =  element_text (size = 10),
          axis.title.x =  element_text (size = 14, vjust = -0.8),
          axis.title.y =  element_text (size = 14),
          axis.text.x =  element_text (size = 8,  vjust = 0.8, angle = 45),
          axis.text.y =  element_text (size = 10))
  ggsave(filename = paste0("pearson_correlation_bAb_IgG_by_visit_in_stage", s, ".pdf"), 
         plot = p1, path = figDir, width = 11.2, height = 15, units = "in")
  
}




titers <- c("original", "B.1.351","BA.1", "BA.2", "BA.4.5")
names(titers) <- c("nAb-ID50 Reference", "nAb-ID50 Beta", "nAb-ID50 BA.1",
                   "nAb-ID50 BA.2", "nAb-ID50 BA.4/BA.5")

for(s in c("Stage 1", "Stage 2")){
  corrMatrix_titer <- tibble("Time" = character(),"plotx" = numeric(),"group" = character(), "markx" = character(), "marky" = character(), "corr" = numeric())
  for(timex in c("B", "Day22", "Day43", "Day78", "Day134", "Day202", "Day292", "Day387")){
    wide_by_variant_titeri <- filter(wide_by_variant_titer, day == timex & Trialstage2 == s)
    for(groupx in c("Non-naïve Vaccine", "Non-naïve Placebo", "Naïve Vaccine", "Naïve Placebo")){
      df <- filter(wide_by_variant_titeri, group == groupx)
      cori <- cor(df[, c("original", "B.1.351","BA.1", "BA.2", "BA.4.5")], 
                  method = "pearson", use = "complete.obs")
      for(kx in 1:(length(titers)-1)){
        for(ky in (kx+1):length(titers)){
          corrMatrix_titer <- add_row(.data = corrMatrix_titer, 
                                      "Time" = ifelse(timex=="B", "D1", gsub("Day","D", timex)),
                                      "plotx" = xmap(timex),
                                      "group" = groupx, "markx" = names(titers)[titers==titers[kx]],
                                      "marky" = names(titers)[titers==titers[ky]], "corr" = cori[titers[kx],titers[ky]])
          
        }
      }
    }
    
  }
  print(min(corrMatrix_titer$corr, na.rm = TRUE))
  corrMatrix_titer$markx <- factor(corrMatrix_titer$markx, levels = names(titers))
  corrMatrix_titer$marky <- factor(corrMatrix_titer$marky, levels = names(titers))
  corrMatrix_titer$markxy <- paste(gsub("nAb-ID50 ","", corrMatrix_titer$markx), "~", 
                                   gsub("nAb-ID50 ","", corrMatrix_titer$marky))
  p1 <- ggplot(data = corrMatrix_titer)+
    geom_point(aes(x = plotx, y = corr, color = group), size = 2, shape = 1, stroke = 1.5)+
    facet_wrap(~ markxy, ncol = 4)+
    scale_x_continuous(breaks = seq(1, 8,1),
                       labels = c("D1", "D22", "D43","D78","D134","D202","D292","D387"), 
                       minor_breaks = NULL)+
    scale_y_continuous(breaks = seq(0, 1, 0.25), limits = c(-0.05,1.05))+
    ylab("Pearson Correlation of log10 nAb-ID50 Titers")+
    xlab("")+
    scale_color_manual(breaks = c("Non-naïve Vaccine", "Non-naïve Placebo", "Naïve Vaccine", "Naïve Placebo"),
                       values = c("gray20", "#1f77b4", "#ff7f0e", "darkolivegreen"))+
    theme_bw()+
    theme(legend.title = element_blank(),
          plot.margin = unit(c(1, 1, 1, 1), "lines"),
          legend.direction = "horizontal",
          legend.key.size=unit(0.5, "cm"),
          legend.text = element_text (size = 12),
          legend.spacing = unit(0.13, "cm"),
          legend.background = element_rect(fill = NA),
          legend.position = "top", legend.justification = c(1, 1),
          title = element_text(vjust = -5, size = 20),
          strip.text.y =  element_text (size = 10),
          strip.text.x =  element_text (size = 10),
          axis.title.x =  element_text (size = 14, vjust = -0.8),
          axis.title.y =  element_text (size = 14),
          axis.text.x =  element_text (size = 8,  vjust = 0.8, angle = 45),
          axis.text.y =  element_text (size = 10))
  ggsave(filename = paste0("pearson_correlation_nAb_ID50_by_visit_in_stage", s, ".pdf"), 
         plot = p1, path = figDir, width = 10, height = 10, units = "in")
}

