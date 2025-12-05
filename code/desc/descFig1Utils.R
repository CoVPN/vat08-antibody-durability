
plotf <- function(data, type = "naive", ylim = c(0, 0.4)){
  data$D01D22OmicronFoldrise <- data$Day22bindSpike_omicron - data$BbindSpike_omicron
  data$D01D22DeltaFoldrise1 <- data$Day22bindSpike_delta1 - data$BbindSpike_delta1
  data$D01D22DeltaFoldrise2 <- data$Day22bindSpike_delta2 - data$BbindSpike_delta2
  data$D01D22DeltaFoldrise3 <- data$Day22bindSpike_delta3 - data$BbindSpike_delta3
  data$D01D22OmicronTiterFoldrise <- data$Day22pseudoneutid50_BA.1 - data$Bpseudoneutid50_BA.1
  data$D01D22RefTiterFoldrise <- data$Day22pseudoneutid50 - data$Bpseudoneutid50
  data$D01D22RefFoldrise <- data$Day22bindSpike - data$BbindSpike
  tab <- tibble("Assay" = character(), "n" = numeric(), "npos" = numeric(),"Fold" = numeric(), "Percentage" = numeric())
  for(f in c(1.25, 1.5, 2, 2.5, 3, 3.5, 4,5)){
    tab <- add_row(tab, "Assay" = "nAb Reference", "n" = sum(!is.na(data$D01D22RefTiterFoldrise)), "Fold" = f,
                   "npos" = sum(data$D01D22RefTiterFoldrise > log10(f), na.rm = TRUE))
    tab <- add_row(tab, "Assay" = "IgG Reference", "n" = sum(!is.na(data$D01D22RefFoldrise)), "Fold" = f,
                   "npos" = sum(data$D01D22RefFoldrise > log10(f), na.rm = TRUE))
    
    tab <- add_row(tab, "Assay" = "IgG Omicron", "n" = sum(!is.na(data$D01D22OmicronFoldrise)), "Fold" = f,
                   "npos" = sum(data$D01D22OmicronFoldrise > log10(f), na.rm = TRUE))
    tab <- add_row(tab, "Assay" = "nAb Omicron", "n" = sum(!is.na(data$D01D22OmicronTiterFoldrise)), "Fold" = f,
                   "npos" = sum(data$D01D22OmicronTiterFoldrise > log10(f), na.rm = TRUE))
    
    tab <- add_row(tab, "Assay" = "IgG Delta-B.1.617.2/AY.4", "n" = sum(!is.na(data$D01D22DeltaFoldrise3)), "Fold" = f,
                   "npos" = sum(data$D01D22DeltaFoldrise3 > log10(f), na.rm = TRUE))
    tab <- add_row(tab, "Assay" = "IgG Delta/IgG Omicron", 
                   "n" = sum(!is.na(data$D01D22DeltaFoldrise3)|!is.na(data$D01D22OmicronFoldrise)), 
                   "Fold" = f,
                   "npos" = sum((data$D01D22DeltaFoldrise3 > log10(f))|data$D01D22OmicronFoldrise > log10(f), na.rm = TRUE))
    tab <- add_row(tab, "Assay" = "IgG Delta/IgG Omicron/nAb Omicron", 
                   "n" = sum(!is.na(data$D01D22DeltaFoldrise3)|!is.na(data$D01D22OmicronFoldrise)|!is.na(data$D01D22OmicronTiterFoldrise)), 
                   "Fold" = f,
                   "npos" = sum((data$D01D22DeltaFoldrise3 > log10(f))|data$D01D22OmicronFoldrise > log10(f)|data$D01D22OmicronTiterFoldrise > log10(f), na.rm = TRUE))
  }
  tab$Percentage <- tab$`npos`/tab$`n`
  print(unique(tab$n))
  if(type == "naive"){
    p <- ggplot(data = tab) + 
      geom_line(aes(x = Fold, y = Percentage, color = Assay), linewidth = 1) +
      geom_point(aes(x = Fold, y = Percentage, shape = Assay), size = 2) + 
      scale_x_continuous(breaks = c(1.25, 1.5, 2, 2.5, 3, 3.5, 4, 5), minor_breaks = NULL) + 
      scale_y_continuous(limits = ylim)+
      xlab("Fold Rise") + 
      ylab("Percentage of Positives\nIncorrectly Diagnosed by Assay Fold Rise") +
      theme_bw()+
      theme(legend.title = element_blank(),
            legend.position = "right",
            legend.direction = "vertical",
            legend.text = element_text (size = 12),
            title = element_text(size = 18),
            strip.text.y =  element_text (size = 14),
            strip.text.x =  element_text (size = 14),
            axis.title.x =  element_text (size = 14),
            axis.title.y =  element_text (size = 14),
            axis.text.x =  element_text (size = 12),
            axis.text.y =  element_text (size = 12))
  }else{
    p <- ggplot(data = tab) + 
      geom_line(aes(x = Fold, y = Percentage, color = Assay), linewidth = 1) +
      geom_point(aes(x = Fold, y = Percentage, shape = Assay), size = 2) + 
      scale_x_continuous(breaks = c(1.25, 1.5, 2, 2.5, 3, 3.5, 4, 5), minor_breaks = NULL) + 
      scale_y_continuous(limits = c(0, 1))+
      xlab("Fold Rise") + 
      ylab("Percentage of Positives\nCorrectly Diagnosed by Assay Fold Rise") +
      theme_bw()+
      theme(legend.title = element_blank(),
            legend.position = "right",
            legend.direction = "vertical",
            legend.text = element_text (size = 12),
            title = element_text(size = 18),
            strip.text.y =  element_text (size = 14),
            strip.text.x =  element_text (size = 14),
            axis.title.x =  element_text (size = 14),
            axis.title.y =  element_text (size = 14),
            axis.text.x =  element_text (size = 12),
            axis.text.y =  element_text (size = 12))
  }
  
  return(p)
}

labf <- function(marker){
  case_when(marker == "bindSpike" ~ "bAb-IgG Spike Reference (AU/ml)",
            marker == "bindSpike_beta" ~ "bAb-IgG Spike Beta (AU/ml)",
            marker == "bindSpike_alpha" ~ "bAb-IgG Spike Alpha (AU/ml)",
            marker == "bindSpike_gamma" ~ "bAb-IgG Spike Gamma (AU/ml)",
            marker == "bindSpike_delta1" ~ "bAb-IgG Spike Delta-AY.4.2 (AU/ml)",
            marker == "bindSpike_delta2" ~ "bAb-IgG Spike Delta-B.1.617.2 (AU/ml)",
            marker == "bindSpike_delta3" ~ "bAb-IgG Spike Delta-B.1.617.2/AY.4 (AU/ml)",
            marker == "bindSpike_omicron" ~ "bAb-IgG Spike Omicron (AU/ml)",
            marker == "pseudoneutid50" ~ "nAb-ID50 Reference",
            marker == "pseudoneutid50_B.1.351" ~ "nAb-ID50 Beta",
            marker == "pseudoneutid50_BA.1" ~ "nAb-ID50 Omicron-BA.1",
            marker == "pseudoneutid50_BA.2" ~ "nAb-ID50 Omicron-BA.2",
            marker == "pseudoneutid50_BA.4.5" ~ "nAb-ID50 Omicron-BA.4/BA.5") 
} 

#LLOQf returns lower limit of quantification for IgG markers but LOD for titers
LLOQf <- function(marker){
  case_when(marker == "bindSpike" ~ log10(563.8),
            marker == "bindSpike_beta" ~ log10(656.0),
            marker == "bindSpike_alpha" ~ log10(573.4),
            marker == "bindSpike_gamma" ~ log10(503.7),
            marker == "bindSpike_delta1" ~ log10(614.2),
            marker == "bindSpike_delta2" ~ log10(559.4),
            marker == "bindSpike_delta3" ~ log10(730.6),
            marker == "bindSpike_omicron" ~ log10(649.2),
            marker == "pseudoneutid50" ~ log10(40), #LOD for titers
            marker == "pseudoneutid50_B.1.351" ~ log10(40), #LOD for titers
            marker == "pseudoneutid50_BA.1" ~ log10(40),#LOD for titers
            marker == "pseudoneutid50_BA.2" ~ log10(40),#LOD for titers
            marker == "pseudoneutid50_BA.4.5" ~ log10(40)) #LOD for titers 
}

ULOQf <- function(marker){
  case_when(marker == "bindSpike" ~ log10(1041179.4),
            marker == "bindSpike_beta" ~ log10(712445.3),
            marker == "bindSpike_alpha" ~ log10(760341.5),
            marker == "bindSpike_gamma" ~ log10(603712.1),
            marker == "bindSpike_delta1" ~ log10(66668.9),
            marker == "bindSpike_delta2" ~ log10(880573.8),
            marker == "bindSpike_delta3" ~ log10(998379.2),
            marker == "bindSpike_omicron" ~ log10(209419.9),
            marker == "pseudoneutid50" ~ log10(198904),
            marker == "pseudoneutid50_B.1.351" ~ log10(35879),
            marker == "pseudoneutid50_BA.1" ~ log10(27104),
            marker == "pseudoneutid50_BA.2" ~ log10(45333),
            marker == "pseudoneutid50_BA.4.5" ~ log10(58293)) 
}

