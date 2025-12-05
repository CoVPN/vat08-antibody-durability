rm(list=ls(all=TRUE))
here::i_am("README.md")
repoDir <- here::here()
datDir <- file.path(repoDir,"data")
codeDir <- file.path(repoDir, "code")
outputDir <- file.path(repoDir, "output")
figDir <- file.path(repoDir, "figures")
tabDir <- file.path(repoDir, "tables")

library(tidyverse)
library(plyr)
source(file.path(codeDir, "common.R"))




plotmarkers <- c("bindSpike", "bindSpike_beta", "bindSpike_alpha", "bindSpike_gamma", "bindSpike_delta1", "bindSpike_delta2", "bindSpike_delta3", "bindSpike_omicron",
                 "pseudoneutid50", "pseudoneutid50_B.1.351", "pseudoneutid50_BA.1", "pseudoneutid50_BA.2", "pseudoneutid50_BA.4.5")
plotmarkers1 <- c("bindSpike", "bindSpike_beta", "bindSpike_alpha", "bindSpike_gamma", "bindSpike_delta1", "bindSpike_delta2", "bindSpike_delta3", "bindSpike_omicron"
)


ppdat <- read.csv(file.path(datDir, "antibodyDurabilityWide_imputed.csv")) 

cases43 <- filter(ppdat, (EventIndFirstInfectionD1 == 1 & (EventTimeFirstInfectionD1 > NumberdaysD1toD43 & EventTimeFirstInfectionD1 < NumberdaysD1toD78))|
                          (EventIndSecondInfectionD1 == 1 & (EventTimeSecondInfectionD1 > NumberdaysD1toD43 & EventTimeSecondInfectionD1 < NumberdaysD1toD78)))
cases78 <- filter(ppdat, (EventIndFirstInfectionD1 == 1 & (EventTimeFirstInfectionD1 > NumberdaysD1toD78 & EventTimeFirstInfectionD1 < NumberdaysD1toD134))|
                          (EventIndSecondInfectionD1 == 1 & (EventTimeSecondInfectionD1 > NumberdaysD1toD78 & EventTimeSecondInfectionD1 < NumberdaysD1toD134)))
cases134 <- filter(ppdat, (EventIndFirstInfectionD1 == 1 & (EventTimeFirstInfectionD1 > NumberdaysD1toD134 & EventTimeFirstInfectionD1 < NumberdaysD1toD202))|
                     (EventIndSecondInfectionD1 == 1 & (EventTimeSecondInfectionD1 > NumberdaysD1toD134 & EventTimeSecondInfectionD1 < NumberdaysD1toD202)))
cases202 <- filter(ppdat, (EventIndFirstInfectionD1 == 1 & (EventTimeFirstInfectionD1 > NumberdaysD1toD202 & EventTimeFirstInfectionD1 < NumberdaysD1toD292))|
                     (EventIndSecondInfectionD1 == 1 & (EventTimeSecondInfectionD1 > NumberdaysD1toD202 & EventTimeSecondInfectionD1 < NumberdaysD1toD292)))
cases292 <- filter(ppdat, (EventIndFirstInfectionD1 == 1 & (EventTimeFirstInfectionD1 > NumberdaysD1toD292 & EventTimeFirstInfectionD1 < NumberdaysD1toD387))|
                     (EventIndSecondInfectionD1 == 1 & (EventTimeSecondInfectionD1 > NumberdaysD1toD292 & EventTimeSecondInfectionD1 < NumberdaysD1toD387)))

#plot antibody markers before and after confirmed infections at during each interval
df_diff <- data.frame()

times <- c(43, 78, 134, 202, 292, 387)
for(i in 1:5){
  time = times[i]
 if(time == 43){
   cases_pooled <- cases43
 }
  if(time == 78){
    cases_pooled <- cases78
  }
  if(time == 134){
    cases_pooled <- cases134
  }
  if(time == 202){
    cases_pooled <- cases202
  }
  if(time == 292){
    cases_pooled <- cases292
  }
  for(Bs in c(0, 1)){
    if(Bs == 0){
      cases <- filter(cases_pooled, Bserostatus == 0)
    }else{
      if(time == 43){
        cases <- filter(cases_pooled, Bserostatus == 1)
      }else{
        cases <- filter(cases_pooled, Bserostatus == 1)
      }
    }
    if(dim(cases)[1] > 0){
      for(marker in c("bindSpike_omicron", "bindSpike_delta3", "pseudoneutid50_BA.4.5")){
        if(i < 5){
          tmp <- cases[, (paste0("Day", times[(i + 1) : (i+2)], marker))]
          markerMax <- apply(tmp, 1, function(x)ifelse(sum(!is.na(x)) >0, max(x, na.rm = TRUE), NA ))
        }else{
          tmp <- cases[, paste0("Day", times[i+1], marker)]
          markerMax <- cases[, paste0("Day", times[i+1], marker)]
        }
        
        df1_diff <- data.frame("antibody_diff"= c(cases[, paste0("Day", times[i+1], marker)]-
                                               cases[, paste0("Day", times[i], marker)]), 
                               "antibody_after"= cases[, paste0("Day", times[i+1], marker)],
                               "antibody_before"= cases[, paste0("Day", times[i], marker)],
                               "antibody_after_max"= markerMax,
                               "antibody_diff_max"= markerMax - cases[, paste0("Day", times[i], marker)],
                               "Interval" = paste0("D", times[i], "-D", times[i+1]),
                               "Bserostatus" = rep(ifelse(Bs == 0, "Naive", "Non-naive"), dim(cases)[1]),
                               "status" = ifelse(Bs == 0 & cases$Trt == 0, "Naïve placebo", "Vaccine/Non-naïve placebo"),
                               "marker" = marker)
        
        df1_diff$antibodyBeforeLessLLOQ <- ifelse(is.na(cases[, paste0("Day", times[i], marker)] ), NA,  ifelse(cases[, paste0("Day", times[i], marker)] < LLOQf(marker), 1, 0))
        df1_diff$antibodyAfterGreaterLLOQ <- ifelse(is.na(markerMax), NA, ifelse(markerMax > LLOQf(marker), 1, 0))
        df1_diff$antibody_conversion <- ifelse(!is.na(df1_diff$antibodyBeforeLessLLOQ) & !is.na(df1_diff$antibodyAfterGreaterLLOQ), 
                                             ifelse((df1_diff$antibodyBeforeLessLLOQ == 1 & df1_diff$antibodyAfterGreaterLLOQ == 1), 1, 0), NA)
        
        df1_diff$antibody_response <- ifelse(!is.na(df1_diff$antibodyBeforeLessLLOQ) & !is.na(df1_diff$antibodyAfterGreaterLLOQ), 
                                             ifelse(df1_diff$antibodyBeforeLessLLOQ == 1, df1_diff$antibody_conversion,
                                                    1*(df1_diff$antibody_after_max - df1_diff$antibody_before >= ifelse(marker == "pseudoneutid50_BA.4.5",log10(4), log10(4))))
                                             , NA)
        
        df_diff <- rbind(df_diff, df1_diff) 
      }}}}


df_diff2  <- df_diff 
df_diff2$Interval <- "Pooled"
df_diff  <- rbind(df_diff , df_diff2 )
df_diff$Interval <- factor(df_diff$Interval, levels = c("D43-D78", "D78-D134", "D134-D202", "D202-D292", "D292-D387","Pooled"))
df_diff$marker2 <- labf(df_diff$marker)




tab1 <- ddply(filter(df_diff, Interval == "Pooled"), .(marker2, Interval), 
              function(x)c(sum(!is.na(x$antibody_conversion) & x$antibodyBeforeLessLLOQ == 1) , 
                           sum(x$antibody_conversion, na.rm = TRUE)/sum(!is.na(x$antibody_conversion) & x$antibodyBeforeLessLLOQ == 1)))


tab2 <- ddply(filter(df_diff, Interval == "Pooled"), .(marker2, Interval), 
              function(x){
                ind = x$antibodyBeforeLessLLOQ == 0
                ans = c(sum(!is.na(x$antibody_diff_max[ind])),10^mean(x$antibody_diff_max[ind], na.rm = TRUE), 
                        10^quantile(x$antibody_diff_max[ind], c(0.25, 0.5, 0.75), na.rm = TRUE))

                
              })



tab1[, -c(1:3)] <- round(tab1[, -c(1:3)], 2)
write.csv(tab1, file.path(tabDir, "cases_response.csv"))

tab2[, -c(1:3)] <- round(tab2[, -c(1:3)], 2)
write.csv(tab2, file.path(tabDir, "cases_foldrise.csv"))


library(scales)
scientific_10 <- function(x) {
  parse(text=gsub("1e\\+*", " 10^", scales::scientific_format()(x)))
}

#plot before and after infection level
data = filter(df_diff, !is.na(antibodyBeforeLessLLOQ) & 
                Interval == "Pooled" & marker2 == "bAb-IgG Spike Delta-B.1.617.2/AY.4")
data$g <- ifelse(data$antibodyBeforeLessLLOQ == 1, "Pre-infection < LLOQ", "Pre-infection >= LLOQ")
p <- ggplot(data = data) +
  geom_point( aes(x = antibody_before, y = antibody_after_max, color = g, shape = status), stroke = 1.2, size = 3) +
  xlab("Pre-infection Antibody Marker Level ") + 
  ylab("Post-infection Antibody Marker Level") +
  ggtitle("bAb-IgG Spike Delta-B.1.617.2/AY.4") +
  scale_color_manual(breaks = c("Pre-infection < LLOQ", "Pre-infection >= LLOQ"),
                     values = c("darkblue", "orange"))+
  scale_shape_manual(breaks = c("Naïve placebo", "Vaccine/Non-naïve placebo"),
                     values = c(1, 2))+
  
  scale_y_continuous(limits = c(1, 6), 
                     breaks = seq(0, 6, 1), 
                     labels = scientific_10(10^seq(0, 6, 1))) +
  scale_x_continuous(limits = c(1, 6), 
                     breaks = seq(0, 6, 1), 
                     labels = scientific_10(10^seq(0, 6, 1))) +
  theme_bw()+guides(color=guide_legend(ncol=2))+
  theme_bw()+guides(shape=guide_legend(ncol=2))+
  theme(legend.title = element_blank(),
        plot.margin = unit(c(1, 1, 1, 1), "lines"),
        legend.direction = "vertical",
        legend.key.size=unit(0.5, "cm"),
        legend.text = element_text (size = 12),
        legend.position = c(0.8, 0.2), 
        legend.justification = c(1, 1),
        legend.spacing = unit(-0.2, "cm"),
        legend.background = element_rect(fill = NA),
        title = element_text(vjust = -5, size = 16),
        strip.text.y =  element_text (size = 16),
        strip.text.x =  element_text (size = 16),
        axis.title.x =  element_text (size = 16, vjust = -0.8),
        axis.title.y =  element_text (size = 16),
        axis.text.x =  element_text (size = 14,  vjust = 0.5),
        axis.text.y =  element_text (size = 14))
ggsave(filename = "pre_vs_post_infection_IgG_delta.pdf", 
       plot = p, path = figDir, width = 7, height = 7, units = "in") 



data = filter(df_diff, !is.na(antibodyBeforeLessLLOQ) & 
                Interval == "Pooled" & marker2 == "bAb-IgG Spike Omicron-B.1.1.529")
data$g <- ifelse(data$antibodyBeforeLessLLOQ == 1, "Pre-infection < LLOQ", "Pre-infection >= LLOQ")
p <- ggplot(data = data) +
  geom_point( aes(x = antibody_before, y = antibody_after_max, color = g, shape = status), stroke = 1.2, size = 3) +
  xlab("Pre-infection Antibody Marker Level ") + 
  ylab("Post-infection Antibody Marker Level") +
  ggtitle("bAb-IgG Spike Omicron-B.1.1.529") +
  scale_color_manual(breaks = c("Pre-infection < LLOQ", "Pre-infection >= LLOQ"),
                     values = c("darkblue", "orange"))+
  scale_shape_manual(breaks = c("Naïve placebo", "Vaccine/Non-naïve placebo"),
                     values = c(1, 2))+
  scale_y_continuous(limits = c(1, 6), 
                     breaks = seq(0, 6, 1), 
                     labels = scientific_10(10^seq(0, 6, 1))) +
  scale_x_continuous(limits = c(1, 6), 
                     breaks = seq(0, 6, 1), 
                     labels = scientific_10(10^seq(0, 6, 1))) +
  theme_bw()+
  theme_bw()+guides(color=guide_legend(ncol=2))+
  theme_bw()+guides(shape=guide_legend(ncol=2))+
  theme(legend.title = element_blank(),
        plot.margin = unit(c(1, 1, 1, 1), "lines"),
        legend.direction = "vertical",
        legend.key.size=unit(0.5, "cm"),
        legend.text = element_text (size = 12),
        legend.position = c(0.8, 0.2), 
        legend.justification = c(1, 1),
        legend.spacing = unit(-0.2, "cm"),
        legend.background = element_rect(fill = NA),
        title = element_text(vjust = -5, size = 16),
        strip.text.y =  element_text (size = 16),
        strip.text.x =  element_text (size = 16),
        axis.title.x =  element_text (size = 16, vjust = -0.8),
        axis.title.y =  element_text (size = 16),
        axis.text.x =  element_text (size = 14,  vjust = 0.5),
        axis.text.y =  element_text (size = 14))
ggsave(filename = "pre_vs_post_infection_IgG_omicron.pdf", 
       plot = p, path = figDir, width = 7, height = 7, units = "in") 




data = filter(df_diff, !is.na(antibodyBeforeLessLLOQ) & 
                Interval == "Pooled" & marker2 == "nAb-ID50 Omicron-BA.4/BA.5")
data$g <- ifelse(data$antibodyBeforeLessLLOQ == 1, "Pre-infection < LOD", "Pre-infection >= LOD")
p <- ggplot(data = data) +
  geom_point( aes(x = antibody_before, y = antibody_after_max, color = g, shape = status), stroke = 1.2, size = 3) +
  xlab("Pre-infection Antibody Marker Level ") + 
  ylab("Post-infection Antibody Marker Level") +
  ggtitle("nAb-ID50 Omicron-BA.4/BA.5") +
  scale_color_manual(breaks = c("Pre-infection < LOD", "Pre-infection >= LOD"),
                     values = c("darkblue", "orange"))+
  scale_shape_manual(breaks = c("Naïve placebo", "Vaccine/Non-naïve placebo"),
                     values = c(1, 2))+
  scale_y_continuous(limits = c(1, 6), 
                     breaks = seq(0, 6, 1), 
                     labels = scientific_10(10^seq(0, 6, 1))) +
  scale_x_continuous(limits = c(1, 6), 
                     breaks = seq(0, 6, 1), 
                     labels = scientific_10(10^seq(0, 6, 1))) +
  theme_bw()+
  theme_bw()+guides(color=guide_legend(ncol=2))+
  theme_bw()+guides(shape=guide_legend(ncol=2))+
  theme(legend.title = element_blank(),
        plot.margin = unit(c(1, 1, 1, 1), "lines"),
        legend.direction = "vertical",
        legend.key.size=unit(0.5, "cm"),
        legend.text = element_text (size = 12),
        legend.position = c(0.9, 0.2), 
        legend.justification = c(1, 1),
        legend.spacing = unit(-0.2, "cm"),
        legend.background = element_rect(fill = NA),
        title = element_text(vjust = -5, size = 16),
        strip.text.y =  element_text (size = 16),
        strip.text.x =  element_text (size = 16),
        axis.title.x =  element_text (size = 16, vjust = -0.8),
        axis.title.y =  element_text (size = 16),
        axis.text.x =  element_text (size = 14,  vjust = 0.5),
        axis.text.y =  element_text (size = 14))
ggsave(filename = "pre_vs_post_infection_nAbID50_omicron.pdf", 
       plot = p, path = figDir, width = 7, height = 7, units = "in") 



