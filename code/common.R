#Shared group definition and contrasts
#Define groups and labels
groups <- tibble("Stage" = numeric(), "Treatment" = numeric(), "Baseline" = character(), "Group" = numeric(), "Label" = character())
groups <- add_row(.data = groups, "Stage" = 1, "Treatment" = 1, "Baseline" = "naive", "Group" = 1, "Label" = "1-V-N")
groups <- add_row(.data = groups, "Stage" = 1, "Treatment" = 1, "Baseline" = "nonnaiveD01SposOnly", "Group" = 2, "Label" = "1-V-NN1")
groups <- add_row(.data = groups, "Stage" = 1, "Treatment" = 0, "Baseline" = "nonnaiveD01SposOnly", "Group" = 3, "Label" = "1-P-NN1")
groups <- add_row(.data = groups, "Stage" = 2, "Treatment" = 1, "Baseline" = "naive" , "Group" = 4, "Label" = "2-V-N")
groups <- add_row(.data = groups, "Stage" = 2, "Treatment" = 1, "Baseline" = "nonnaiveD01SposOnly" , "Group" = 5, "Label" = "2-V-NN1")
groups <- add_row(.data = groups, "Stage" = 2, "Treatment" = 1, "Baseline" = "nonnaiveElse" , "Group" = 6, "Label" = "2-V-NN2")
groups <- add_row(.data = groups, "Stage" = 2, "Treatment" = 0, "Baseline" = "nonnaiveD01SposOnly" , "Group" = 7, "Label" = "2-P-NN1")
groups <- add_row(.data = groups, "Stage" = 2, "Treatment" = 0, "Baseline" = "nonnaiveElse", "Group" = 8, "Label" = "2-P-NN2")
groups <- add_row(.data = groups, "Stage" = 1, "Treatment" = 1, "Baseline" = "nonnaive", "Group" = 9, "Label" = "1-V-NN")
groups <- add_row(.data = groups, "Stage" = 1, "Treatment" = 0, "Baseline" = "nonnaive", "Group" = 10, "Label" = "1-P-NN")
groups <- add_row(.data = groups, "Stage" = 2, "Treatment" = 1, "Baseline" = "nonnaive", "Group" = 11, "Label" = "2-V-NN")
groups <- add_row(.data = groups, "Stage" = 2, "Treatment" = 0, "Baseline" = "nonnaive", "Group" = 12, "Label" = "2-P-NN")

contrastLabel <- function(GA, GB, groups){
  GAlabel <- groups$Label[groups$Group == GA]
  GBlabel <- groups$Label[groups$Group == GB]
  return(paste0(GAlabel, " vs. ", GBlabel))
}

#contrast of interest
contrasts <- tibble("GA" = numeric(), "GB" = numeric(), "label" = character())
contrasts <- add_row(contrasts, "GA" = 9, "GB" = 1, "label" = contrastLabel(9, 1, groups))
contrasts <- add_row(contrasts, "GA" = 5, "GB" = 6, "label" = contrastLabel(5, 6, groups))
contrasts <- add_row(contrasts, "GA" = 11, "GB" = 4, "label" = contrastLabel(11, 4, groups))
contrasts <- add_row(contrasts, "GA" = 7, "GB" = 8, "label" = contrastLabel(7, 8, groups))
contrasts <- add_row(contrasts, "GA" = 9, "GB" = 10, "label" = contrastLabel(9, 10, groups))
contrasts <- add_row(contrasts, "GA" = 11, "GB" = 12, "label" = contrastLabel(11, 12, groups))
contrasts <- add_row(contrasts, "GA" = 4, "GB" = 1, "label" = contrastLabel(4, 1, groups))
contrasts <- add_row(contrasts, "GA" = 5, "GB" = 2, "label" = contrastLabel(5, 2, groups))
contrasts <- add_row(contrasts, "GA" = 11, "GB" = 9, "label" = contrastLabel(11, 9, groups))
contrasts <- add_row(contrasts, "GA" = 12, "GB" = 10, "label" = contrastLabel(12, 10, groups))
contrasts <- add_row(contrasts, "GA" = 7, "GB" = 3, "label" = contrastLabel(7, 3, groups))


labf <- function(marker){
  case_when(marker == "bindSpike" ~ "bAb-IgG Spike Index",
            marker == "bindSpike_beta" ~ "bAb-IgG Spike Beta",
            marker == "bindSpike_alpha" ~ "bAb-IgG Spike Alpha",
            marker == "bindSpike_gamma" ~ "bAb-IgG Spike Gamma",
            marker == "bindSpike_delta1" ~ "bAb-IgG Spike Delta-AY.4",
            marker == "bindSpike_delta2" ~ "bAb-IgG Spike Delta-B.1.617.2",
            marker == "bindSpike_delta3" ~ "bAb-IgG Spike Delta-B.1.617.2/AY.4",
            marker == "bindSpike_omicron" ~ "bAb-IgG Spike Omicron",
            marker == "pseudoneutid50" ~ "nAb-ID50 Reference",
            marker == "pseudoneutid50_B.1.351" ~ "nAb-ID50 Beta",
            marker == "pseudoneutid50_BA.1" ~ "nAb-ID50 Omicron-BA.1",
            marker == "pseudoneutid50_BA.2" ~ "nAb-ID50 Omicron-BA.2",
            marker == "pseudoneutid50_BA.4.5" ~ "nAb-ID50 Omicron-BA.4/BA.5") 
}

plot_labf <- function(marker){
  case_when(marker == "bindSpike" ~ "bAb-IgG Spike Index (AU/ml)",
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
            marker == "pseudoneutid50_BA.4.5" ~ log10(40) #LOD for titers
            
  ) 
}

vars <- c("Ptid","Trt", "Trialstage", "Sex", "Age", "BMI", "Age", "Bserostatus", "EventTimeFirstInfectionD1", "EventIndFirstInfectionD1",
          "EventTimeFirstInfectionD43", "EventIndFirstInfectionD43", "FirstEnrollmentDate", "Perprotocol", 
          "D01_S_pos_only_in_non_naive_group","NumberdaysD1toD22", "NumberdaysD1toD43", "NumberdaysD1toD78", "NumberdaysD1toD134", "NumberdaysD1toD202", "NumberdaysD1toD292",
          "NumberdaysD1toD387", "D01_S_pos", "D01_N_pos", "D01_NAAT_pos", "riskScore", "time", "visitn",
          "BbindSpike", "BbindSpike_beta", "BbindSpike_alpha", "BbindSpike_gamma", "BbindSpike_delta1", "BbindSpike_delta2", "BbindSpike_delta3", "BbindSpike_omicron",
          "Bpseudoneutid50", "Bpseudoneutid50_B.1.351", "Bpseudoneutid50_BA.1", "Bpseudoneutid50_BA.2", "Bpseudoneutid50_BA.4.5")


