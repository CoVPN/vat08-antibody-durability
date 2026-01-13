
format.p2 <- function(x){
  y = NULL
  if(grepl("<",as.character(x))){
    y <- paste0("< ", gsub("<","", x))
  }else{
    y  = paste0("= ", x)
  }
  
  return(y)
}

reformatRatio <- function(x){
  y <- strsplit(x, split = " ")
  y1 <- y[[1]][1]
  y2 <- strsplit(strsplit(y[[1]][2], split = "\\(")[[1]][2], split = "\\,")[[1]][1]
  y3 <- strsplit(y[[1]][3], split = "\\)")[[1]][1]
  
  est <- format(round(as.numeric(y1),2), nsmall = 2, digits = 2)
  ll <- format(round(as.numeric(y2),2), nsmall = 2, digits = 2)
  ul <- format(round(as.numeric(y3),2), nsmall = 2, digits = 2)
  return(paste0(est, " (", ll, ", ", ul, ")", sep = ""))
}

plotSummary_D43 <- function(comparisons, labels.x = c("Monovalent\nVaccine", "Stage 1\nPlacebo", "Bivalent\nVaccine", "Stage 2\nPlacebo"), figureLabel,
                            IgG_ylim = c(2.1, 6.5), ID50_ylim = c(0.6, 5.1), reverseGroup = FALSE){
 
  for(markeri in c("bAb-IgG Spike Omicron", "bAb-IgG Spike Index", "nAb-ID50 Omicron-BA.4/BA.5",
                   "nAb-ID50 Reference")){
    if(markeri %in% c("bAb-IgG Spike Omicron", "bAb-IgG Spike Index")){
      rangeData <- data.frame(cbind(c(1, 5), IgG_ylim))
      range.y = IgG_ylim[2] - IgG_ylim[1]
    }else{
      rangeData <- data.frame(cbind(c(1, 5), ID50_ylim))
      range.y = ID50_ylim[2] - ID50_ylim[1]
    }
    
    p <- ggplot(data = rangeData)+
      scale_x_continuous(limits = c(0.55, 4.4), minor_breaks = NULL, breaks = c(1, 2, 3, 4), 
                         labels = labels.x)+
      scale_y_continuous(limits = c(rangeData[1, 2], rangeData[2, 2]), minor_breaks = NULL, breaks = seq(1, 7, 1), labels = scientific_10(10^seq(1, 7, 1)))+
      ylab("D43 Geometric Mean")+
      xlab("")+
      ggtitle(ifelse(markeri%in% c("bAb-IgG Spike Omicron", "bAb-IgG Spike Index"), paste0(markeri, " (AU/ml)"), markeri))+
      coord_cartesian( clip = "off")+
      theme_bw()+
      theme(plot.margin = margin(t = 0.5, r = 0.5, b = 0.5, l = 0.5, unit = "cm"),
            title = element_text(size = 16),
            axis.text.x = element_text(size = 14),
            axis.text.y = element_text(size = 14),
            axis.title.x = element_text(size = 16, margin = margin(t = 10)),
            axis.title.y = element_text(size = 16, margin = margin(r = 10)))
    if(reverseGroup){
      loci <- c(2, 1)
    }else{
      loci <- c(1, 2)
    }
    
    permuTests_bootCI_j <- filter(output, Comparison== comparisons[1] & Marker == markeri)
    
    maxy <- max(log10(permuTests_bootCI_j$D43_GA_CU), log10(permuTests_bootCI_j$D43_GB_CU))
    #browser()
    p <- p+geom_segment(aes(x = loci[1], xend = loci[1], y = log10(permuTests_bootCI_j$D43_GA_CL), yend = log10(permuTests_bootCI_j$D43_GA_CU)), linewidth = 1)+
      geom_segment(aes(x = loci[1]-0.05, xend = loci[1] + 0.05, y = log10(permuTests_bootCI_j$D43_GA_CL), yend = log10(permuTests_bootCI_j$D43_GA_CL)), linewidth = 1)+
      geom_segment(aes(x = loci[1]-0.05, xend = loci[1] + 0.05, y = log10(permuTests_bootCI_j$D43_GA_CU), yend = log10(permuTests_bootCI_j$D43_GA_CU)), linewidth = 1)+
      geom_point(aes(x = loci[1],  y = log10(permuTests_bootCI_j$D43_GA_est)), size = 1.7)+
      annotate (geom = "text", x = loci[1]-0.32, y = log10(permuTests_bootCI_j$D43_GA_est), label = comma(round(permuTests_bootCI_j$D43_GA_est)), size = 4.5)+
      annotate (geom = "text", x = loci[1]-0.32, y = log10(permuTests_bootCI_j$D43_GA_CL)-
                  ifelse((log10(permuTests_bootCI_j$D43_GA_est)-log10(permuTests_bootCI_j$D43_GA_CL))<range.y*0.035, range.y*0.035, 0), 
                label = comma(round(permuTests_bootCI_j$D43_GA_CL)), size = 4.5)+
      annotate (geom = "text", x = loci[1]-0.32, y = log10(permuTests_bootCI_j$D43_GA_CU)+
                  ifelse((log10(permuTests_bootCI_j$D43_GA_CU)-log10(permuTests_bootCI_j$D43_GA_est))<range.y*0.035, range.y*0.035, 0), 
                label = comma(round(permuTests_bootCI_j$D43_GA_CU)), size = 4.5)+
      
      geom_segment(aes(x = loci[2], xend = loci[2], y = log10(permuTests_bootCI_j$D43_GB_CL), yend = log10(permuTests_bootCI_j$D43_GB_CU)), linewidth = 1)+
      geom_segment(aes(x = loci[2]-0.05, xend = loci[2] + 0.05, y = log10(permuTests_bootCI_j$D43_GB_CL), yend = log10(permuTests_bootCI_j$D43_GB_CL)), linewidth = 1)+
      geom_segment(aes(x = loci[2]-0.05, xend = loci[2] + 0.05, y = log10(permuTests_bootCI_j$D43_GB_CU), yend = log10(permuTests_bootCI_j$D43_GB_CU)), linewidth = 1)+
      geom_point(aes(x = loci[2],  y = log10(permuTests_bootCI_j$D43_GB_est)), size = 1.7)+
      annotate (geom = "text", x = loci[2]-0.32, y = log10(permuTests_bootCI_j$D43_GB_est), label = comma(round(permuTests_bootCI_j$D43_GB_est)), size = 4.5)+
      annotate (geom = "text", x = loci[2]-0.32, y = log10(permuTests_bootCI_j$D43_GB_CL)-
                  ifelse((log10(permuTests_bootCI_j$D43_GB_est)-log10(permuTests_bootCI_j$D43_GB_CL))<range.y*0.035, range.y*0.035, 0), 
                label = comma(round(permuTests_bootCI_j$D43_GB_CL)), size = 4.5)+
      annotate (geom = "text", x = loci[2]-0.32, y = log10(permuTests_bootCI_j$D43_GB_CU)+
                  ifelse((log10(permuTests_bootCI_j$D43_GB_CU)-log10(permuTests_bootCI_j$D43_GB_est))<range.y*0.035, range.y*0.035, 0), 
                label = comma(round(permuTests_bootCI_j$D43_GB_CU)), size = 4.5)+
      geom_segment(aes(x = loci[1], xend = loci[1], y = maxy + range.y*0.05, yend = maxy + range.y*0.07))+
      geom_segment(aes(x = loci[2], xend = loci[2], y = maxy + range.y*0.05, yend = maxy + range.y*0.07))+
      geom_segment(aes(x = loci[1], xend = loci[2], y = maxy + range.y*0.07, yend = maxy + range.y*0.07))+
      annotate (geom = "text", x = (loci[1] + loci[2])/2, y = maxy + range.y*0.15, label = paste0("GM Ratio: ", reformatRatio(permuTests_bootCI_j$D43_diff), "\nUnadj. P ",
                                                                                         format.p2(permuTests_bootCI_j$D43_Pvalue), "; FWER P ",
                                                                                         format.p2(permuTests_bootCI_j$D43_FWER_Pvalue)), size = 4.5)+
      annotate (geom = "text", x = loci[1], y = rangeData[1, 2], label = paste0("N = ", permuTests_bootCI_j$nA), size = 4.5)+
      annotate (geom = "text", x = loci[2], y = rangeData[1, 2], label = paste0("N = ", permuTests_bootCI_j$nB), size = 4.5)
    
    loci2 <- loci + 2
    permuTests_bootCI_j2 <- filter(output, Comparison== comparisons[2] & Marker == markeri)
    maxy2 <- max(log10(permuTests_bootCI_j2$D43_GA_CU), log10(permuTests_bootCI_j2$D43_GB_CU))
    p <- p+geom_segment(aes(x = loci2[1], xend = loci2[1], y = log10(permuTests_bootCI_j2$D43_GA_CL), yend = log10(permuTests_bootCI_j2$D43_GA_CU)), linewidth = 1)+
      geom_segment(aes(x = loci2[1]-0.05, xend = loci2[1] + 0.05, y = log10(permuTests_bootCI_j2$D43_GA_CL), yend = log10(permuTests_bootCI_j2$D43_GA_CL)), linewidth = 1)+
      geom_segment(aes(x = loci2[1]-0.05, xend = loci2[1] + 0.05, y = log10(permuTests_bootCI_j2$D43_GA_CU), yend = log10(permuTests_bootCI_j2$D43_GA_CU)), linewidth = 1)+
      geom_point(aes(x = loci2[1],  y = log10(permuTests_bootCI_j2$D43_GA_est)), size = 1.7)+
      
      annotate (geom = "text", x = loci2[1]-0.32, y = log10(permuTests_bootCI_j2$D43_GA_est), label = comma(round(permuTests_bootCI_j2$D43_GA_est)), size = 4.5)+
      annotate (geom = "text", x = loci2[1]-0.32, y = log10(permuTests_bootCI_j2$D43_GA_CL)-
                  ifelse((log10(permuTests_bootCI_j2$D43_GA_est)-log10(permuTests_bootCI_j2$D43_GA_CL))<range.y*0.035, range.y*0.035, 0), 
                label = comma(round(permuTests_bootCI_j2$D43_GA_CL)), size = 4.5)+
      annotate (geom = "text", x = loci2[1]-0.32, y = log10(permuTests_bootCI_j2$D43_GA_CU)+
                  ifelse((log10(permuTests_bootCI_j2$D43_GA_CU)-log10(permuTests_bootCI_j2$D43_GA_est))<range.y*0.035, range.y*0.035, 0), 
                label = comma(round(permuTests_bootCI_j2$D43_GA_CU)), size = 4.5)+
      
      geom_segment(aes(x = loci2[2], xend = loci2[2], y = log10(permuTests_bootCI_j2$D43_GB_CL), yend = log10(permuTests_bootCI_j2$D43_GB_CU)), linewidth = 1)+
      geom_segment(aes(x = loci2[2]-0.05, xend = loci2[2] + 0.05, y = log10(permuTests_bootCI_j2$D43_GB_CL), yend = log10(permuTests_bootCI_j2$D43_GB_CL)), linewidth = 1)+
      geom_segment(aes(x = loci2[2]-0.05, xend = loci2[2] + 0.05, y = log10(permuTests_bootCI_j2$D43_GB_CU), yend = log10(permuTests_bootCI_j2$D43_GB_CU)), linewidth = 1)+
      geom_point(aes(x = loci2[2],  y = log10(permuTests_bootCI_j2$D43_GB_est)), size = 1.7)+
      
      annotate (geom = "text", x = loci2[2]-0.32, y = log10(permuTests_bootCI_j2$D43_GB_est), label = comma(round(permuTests_bootCI_j2$D43_GB_est)), size = 4.5)+
      annotate (geom = "text", x = loci2[2]-0.32, y = log10(permuTests_bootCI_j2$D43_GB_CL)-
                  ifelse((log10(permuTests_bootCI_j2$D43_GB_est)-log10(permuTests_bootCI_j2$D43_GB_CL))<range.y*0.035, range.y*0.035, 0), 
                label = comma(round(permuTests_bootCI_j2$D43_GB_CL)), size = 4.5)+
      annotate (geom = "text", x = loci2[2]-0.32, y = log10(permuTests_bootCI_j2$D43_GB_CU)+
                  ifelse((log10(permuTests_bootCI_j2$D43_GB_CU)-log10(permuTests_bootCI_j2$D43_GB_est))<range.y*0.035, range.y*0.035, 0), 
                label = comma(round(permuTests_bootCI_j2$D43_GB_CU)), size = 4.5)+
      
      geom_segment(aes(x = loci2[1], xend = loci2[1], y = maxy2 + range.y*0.05, yend = maxy2 + range.y*0.07))+
      geom_segment(aes(x = loci2[2], xend = loci2[2], y = maxy2 + range.y*0.05, yend = maxy2 + range.y*0.07))+
      geom_segment(aes(x = loci2[1], xend = loci2[2], y = maxy2 + range.y*0.07, yend = maxy2 + range.y*0.07))+
      annotate (geom = "text", x = (loci2[1] + loci2[2])/2, y = maxy2 + range.y*0.15, label = paste0("GM Ratio: ", reformatRatio(permuTests_bootCI_j2$D43_diff), "\nUnadj. P ",
                                                                                            format.p2(permuTests_bootCI_j2$D43_Pvalue), "; FWER P ",
                                                                                            format.p2(permuTests_bootCI_j2$D43_FWER_Pvalue)), size = 4.5)+
      annotate (geom = "text", x = loci2[1], y = rangeData[1, 2], label = paste0("N = ", permuTests_bootCI_j2$nA), size = 4.5)+
      annotate (geom = "text", x = loci2[2], y = rangeData[1, 2], label = paste0("N = ", permuTests_bootCI_j2$nB), size = 4.5)
    
    
    ggsave(filename = paste0("D43_", gsub("/","", gsub(" ", "",markeri)), figureLabel, "summary.pdf"), plot = p, path = figDir, width = 7.3, height = 6, units = "in") 
    
  }
}

plotSummary_AUC <- function(comparisons, labels.x = c("Monovalent\nVaccine", "Stage 1\nPlacebo", "Bivalent\nVaccine", "Stage 2\nPlacebo"), figureLabel,
                            IgG_ylim = c(2.1, 6.5), ID50_ylim = c(0.6, 5.1), reverseGroup = FALSE){
  #browser()
  for(markeri in c("bAb-IgG Spike Omicron", "bAb-IgG Spike Index", "nAb-ID50 Omicron-BA.4/BA.5",
                   "nAb-ID50 Reference")){
    if(markeri %in% c("bAb-IgG Spike Omicron", "bAb-IgG Spike Index")){
      rangeData <- data.frame(cbind(c(1, 5), IgG_ylim))
      range.y = IgG_ylim[2] - IgG_ylim[1]
    }else{
      rangeData <- data.frame(cbind(c(1, 5), ID50_ylim))
      range.y = ID50_ylim[2] - ID50_ylim[1]
    }
    
    p <- ggplot(data = rangeData)+
      scale_x_continuous(limits = c(0.55, 4.4), minor_breaks = NULL, breaks = c(1, 2, 3, 4), 
                         labels = labels.x)+
      scale_y_continuous(limits = c(rangeData[1, 2], rangeData[2, 2]), minor_breaks = NULL, breaks = seq(1, 7, 1), labels = scientific_10(10^seq(1, 7, 1)))+
      ylab("Durability (D43 to D202)")+
      xlab("")+
      ggtitle(ifelse(markeri%in% c("bAb-IgG Spike Omicron", "bAb-IgG Spike Index"), paste0(markeri, " (AU/ml)"), markeri))+
      coord_cartesian( clip = "off")+
      theme_bw()+
      theme(plot.margin = margin(t = 0.5, r = 0.5, b = 0.5, l = 0.5, unit = "cm"),
            title = element_text(size = 16),
            axis.text.x = element_text(size = 14),
            axis.text.y = element_text(size = 14),
            axis.title.x = element_text(size = 16, margin = margin(t = 10)),
            axis.title.y = element_text(size = 16, margin = margin(r = 10)))
    
    if(reverseGroup){
      loci <- c(2, 1)
    }else{
      loci <- c(1, 2)
    }
    permuTests_bootCI_j <- filter(output, Comparison== comparisons[1] & Marker == markeri)
    maxy <- max(permuTests_bootCI_j$AUC_GA_CU, permuTests_bootCI_j$AUC_GB_CU)
    p <- p+geom_segment(aes(x = loci[1], xend = loci[1], y = (permuTests_bootCI_j$AUC_GA_CL), yend = (permuTests_bootCI_j$AUC_GA_CU)), linewidth = 1)+
      geom_segment(aes(x = loci[1]-0.05, xend = loci[1] + 0.05, y = (permuTests_bootCI_j$AUC_GA_CL), yend = (permuTests_bootCI_j$AUC_GA_CL)), linewidth = 1)+
      geom_segment(aes(x = loci[1]-0.05, xend = loci[1] + 0.05, y = (permuTests_bootCI_j$AUC_GA_CU), yend = (permuTests_bootCI_j$AUC_GA_CU)), linewidth = 1)+
      geom_point(aes(x = loci[1],  y = (permuTests_bootCI_j$AUC_GA_est)), size = 1.7)+
      
      annotate (geom = "text", x = loci[1]-0.28, y = permuTests_bootCI_j$AUC_GA_est, label = comma(round(10^permuTests_bootCI_j$AUC_GA_est)), size = 4.5)+
      annotate (geom = "text", x = loci[1]-0.28, y = permuTests_bootCI_j$AUC_GA_CL-
                  ifelse((permuTests_bootCI_j$AUC_GA_est-permuTests_bootCI_j$AUC_GA_CL)<range.y*0.035, range.y*0.035, 0), 
                label = comma(round(10^permuTests_bootCI_j$AUC_GA_CL)), size = 4.5)+
      annotate (geom = "text", x = loci[1]-0.28, y = permuTests_bootCI_j$AUC_GA_CU+
                  ifelse((permuTests_bootCI_j$AUC_GA_CU-permuTests_bootCI_j$AUC_GA_est)<range.y*0.035, range.y*0.035, 0), 
                label = comma(round(10^permuTests_bootCI_j$AUC_GA_CU)), size = 4.5)+
      
     
      geom_segment(aes(x = loci[2], xend = loci[2], y = (permuTests_bootCI_j$AUC_GB_CL), yend = (permuTests_bootCI_j$AUC_GB_CU)), linewidth = 1)+
      geom_segment(aes(x = loci[2]-0.05, xend = loci[2] + 0.05, y = (permuTests_bootCI_j$AUC_GB_CL), yend = (permuTests_bootCI_j$AUC_GB_CL)), linewidth = 1)+
      geom_segment(aes(x = loci[2]-0.05, xend = loci[2] + 0.05, y = (permuTests_bootCI_j$AUC_GB_CU), yend = (permuTests_bootCI_j$AUC_GB_CU)), linewidth = 1)+
      geom_point(aes(x = loci[2],  y = (permuTests_bootCI_j$AUC_GB_est)), size = 1.7)+
      
      annotate (geom = "text", x = loci[2]-0.32, y = permuTests_bootCI_j$AUC_GB_est, label = comma(round(10^permuTests_bootCI_j$AUC_GB_est)), size = 4.5)+
      annotate (geom = "text", x = loci[2]-0.32, y = permuTests_bootCI_j$AUC_GB_CL-
                  ifelse((permuTests_bootCI_j$AUC_GB_est-permuTests_bootCI_j$AUC_GB_CL)<range.y*0.035, range.y*0.035, 0), 
                label = comma(round(10^permuTests_bootCI_j$AUC_GB_CL)), size = 4.5)+
      annotate (geom = "text", x = loci[2]-0.32, y = permuTests_bootCI_j$AUC_GB_CU+
                  ifelse((permuTests_bootCI_j$AUC_GB_CU-permuTests_bootCI_j$AUC_GB_est)<range.y*0.035, range.y*0.035, 0), 
                label = comma(round(10^permuTests_bootCI_j$AUC_GB_CU)), size = 4.5)+
      
      geom_segment(aes(x = loci[1], xend = loci[1], y = maxy + range.y*0.05, yend = maxy + range.y*0.07))+
      geom_segment(aes(x = loci[2], xend = loci[2], y = maxy + range.y*0.05, yend = maxy + range.y*0.07))+
      geom_segment(aes(x = loci[1], xend = loci[2], y = maxy + range.y*0.07, yend = maxy + range.y*0.07))+
      annotate (geom = "text", x = (loci[1] + loci[2])/2, y = maxy+ range.y*0.15, label = paste0("Ratio: ", reformatRatio(permuTests_bootCI_j$AUC_diff), "\nUnadj. P ",
                                                                                                 format.p2(permuTests_bootCI_j$Pvalue), "; FWER P ",
                                                                                                 format.p2(permuTests_bootCI_j$FWER_Pvalue)), size = 4.5)+
      annotate (geom = "text", x = loci[1], y = rangeData[1, 2], label = paste0("N = ", permuTests_bootCI_j$nA), size = 4.5)+
      annotate (geom = "text", x = loci[2], y = rangeData[1, 2], label = paste0("N = ", permuTests_bootCI_j$nB), size = 4.5)
    
    
    loci2 <- loci + 2
    permuTests_bootCI_j2 <- filter(output, Comparison== comparisons[2] & Marker == markeri)
    maxy2 <- max(permuTests_bootCI_j2$AUC_GA_CU, permuTests_bootCI_j2$AUC_GB_CU)
    p <- p+geom_segment(aes(x = loci2[1], xend = loci2[1], y = (permuTests_bootCI_j2$AUC_GA_CL), yend = (permuTests_bootCI_j2$AUC_GA_CU)), linewidth = 1)+
      geom_segment(aes(x = loci2[1]-0.05, xend = loci2[1] + 0.05, y = (permuTests_bootCI_j2$AUC_GA_CL), yend = (permuTests_bootCI_j2$AUC_GA_CL)), linewidth = 1)+
      geom_segment(aes(x = loci2[1]-0.05, xend = loci2[1] + 0.05, y = (permuTests_bootCI_j2$AUC_GA_CU), yend = (permuTests_bootCI_j2$AUC_GA_CU)), linewidth = 1)+
      geom_point(aes(x = loci2[1],  y = (permuTests_bootCI_j2$AUC_GA_est)), size = 1.7)+
     
      annotate (geom = "text", x = loci2[1]-0.28, y = permuTests_bootCI_j2$AUC_GA_est, label = comma(round(10^permuTests_bootCI_j2$AUC_GA_est)), size = 4.5)+
      annotate (geom = "text", x = loci2[1]-0.28, y = permuTests_bootCI_j2$AUC_GA_CL-
                  ifelse((permuTests_bootCI_j2$AUC_GA_est-permuTests_bootCI_j2$AUC_GA_CL)<range.y*0.035, range.y*0.035, 0), 
                label = comma(round(10^permuTests_bootCI_j2$AUC_GA_CL)), size = 4.5)+
      annotate (geom = "text", x = loci2[1]-0.28, y = permuTests_bootCI_j2$AUC_GA_CU+
                  ifelse((permuTests_bootCI_j2$AUC_GA_CU-permuTests_bootCI_j2$AUC_GA_est)<range.y*0.035, range.y*0.035, 0), 
                label = comma(round(10^permuTests_bootCI_j2$AUC_GA_CU)), size = 4.5)+
      
      geom_segment(aes(x = loci2[2], xend = loci2[2], y = (permuTests_bootCI_j2$AUC_GB_CL), yend = (permuTests_bootCI_j2$AUC_GB_CU)), linewidth = 1)+
      geom_segment(aes(x = loci2[2]-0.05, xend = loci2[2] + 0.05, y = (permuTests_bootCI_j2$AUC_GB_CL), yend = (permuTests_bootCI_j2$AUC_GB_CL)), linewidth = 1)+
      geom_segment(aes(x = loci2[2]-0.05, xend = loci2[2] + 0.05, y = (permuTests_bootCI_j2$AUC_GB_CU), yend = (permuTests_bootCI_j2$AUC_GB_CU)), linewidth = 1)+
      geom_point(aes(x = loci2[2],  y = (permuTests_bootCI_j2$AUC_GB_est)), size = 1.7)+
      
      annotate (geom = "text", x = loci2[2]-0.32, y = permuTests_bootCI_j2$AUC_GB_est, label = comma(round(10^permuTests_bootCI_j2$AUC_GB_est)), size = 4.5)+
      annotate (geom = "text", x = loci2[2]-0.32, y = permuTests_bootCI_j2$AUC_GB_CL-
                  ifelse((permuTests_bootCI_j2$AUC_GB_est-permuTests_bootCI_j2$AUC_GB_CL)<range.y*0.035, range.y*0.035, 0), 
                label = comma(round(10^permuTests_bootCI_j2$AUC_GB_CL)), size = 4.5)+
      annotate (geom = "text", x = loci2[2]-0.32, y = permuTests_bootCI_j2$AUC_GB_CU+
                  ifelse((permuTests_bootCI_j2$AUC_GB_CU-permuTests_bootCI_j2$AUC_GB_est)<range.y*0.035, range.y*0.035, 0), 
                label = comma(round(10^permuTests_bootCI_j2$AUC_GB_CU)), size = 4.5)+
      
      geom_segment(aes(x = loci2[1], xend = loci2[1], y = maxy2 + range.y*0.05, yend = maxy2 + range.y*0.07))+
      geom_segment(aes(x = loci2[2], xend = loci2[2], y = maxy2 + range.y*0.05, yend = maxy2 + range.y*0.07))+
      geom_segment(aes(x = loci2[1], xend = loci2[2], y = maxy2 + range.y*0.07, yend = maxy2 + range.y*0.07))+
      annotate (geom = "text", x = (loci2[1] + loci2[2])/2, y = maxy2 + range.y*0.15, label = paste0("Ratio: ", reformatRatio(permuTests_bootCI_j2$AUC_diff), "\nUnadj. P ",
                                                                                                   format.p2(permuTests_bootCI_j2$Pvalue), "; FWER P ",
                                                                                                   format.p2(permuTests_bootCI_j2$FWER_Pvalue)), size = 4.5)+
      annotate (geom = "text", x = loci2[1], y = rangeData[1, 2], label = paste0("N = ", permuTests_bootCI_j2$nA), size = 4.5)+
      annotate (geom = "text", x = loci2[2], y = rangeData[1, 2], label = paste0("N = ", permuTests_bootCI_j2$nB), size = 4.5)
    
    
    ggsave(filename = paste0("AUC_", gsub("/","", gsub(" ", "",markeri)),figureLabel, "summary.pdf"), plot = p, path = figDir, width = 7.3, height = 6, units = "in") 
    
  }
}


plotSummary_GMR <- function(comparisons, labels.x = c("Monovalent\nVaccine", "Stage 1\nPlacebo", "Bivalent\nVaccine", "Stage 2\nPlacebo"), figureLabel,
                            IgG_ylim = c(0, 1.2), ID50_ylim = c(0, 1.2), reverseGroup = FALSE){
  
  for(markeri in c("bAb-IgG Spike Omicron", "bAb-IgG Spike Index", "nAb-ID50 Omicron-BA.4/BA.5",
                   "nAb-ID50 Reference")){
    if(markeri %in% c("bAb-IgG Spike Omicron", "bAb-IgG Spike Index")){
      rangeData <- data.frame(cbind(c(1, 5), IgG_ylim))
      range.y = IgG_ylim[2] - IgG_ylim[1]
    }else{
      rangeData <- data.frame(cbind(c(1, 5), ID50_ylim))
      range.y = ID50_ylim[2] - ID50_ylim[1]
    }
    
    p <- ggplot(data = rangeData)+
      scale_x_continuous(limits = c(0.55, 4.4), minor_breaks = NULL, breaks = c(1, 2, 3, 4), 
                         labels = labels.x)+
      scale_y_continuous(limits = c(rangeData[1, 2], rangeData[2, 2]), minor_breaks = NULL, breaks = seq(0.2, 1, 0.2), 
                         labels = seq(0.2, 1, 0.2))+
      ylab("D202-to-D43 Geometric Mean Ratio")+
      xlab("")+
      ggtitle(ifelse(markeri%in% c("bAb-IgG Spike Omicron", "bAb-IgG Spike Index"), paste0(markeri, " (AU/ml)"), markeri))+
      coord_cartesian( clip = "off")+
      theme_bw()+
      theme(plot.margin = margin(t = 0.5, r = 0.5, b = 0.5, l = 0.5, unit = "cm"),
            title = element_text(size = 16),
            axis.text.x = element_text(size = 14),
            axis.text.y = element_text(size = 14),
            axis.title.x = element_text(size = 16, margin = margin(t = 10)),
            axis.title.y = element_text(size = 16, margin = margin(r = 10)))
    
    if(reverseGroup){
      loci <- c(2, 1)
    }else{
      loci <- c(1, 2)
    }
    permuTests_bootCI_j <- filter(output, Comparison== comparisons[1] & Marker == markeri)
    maxy <- max(permuTests_bootCI_j$rate_GA_CU, permuTests_bootCI_j$rate_GB_CU)
    p <- p+geom_segment(aes(x = loci[1], xend = loci[1], y = (permuTests_bootCI_j$rate_GA_CL), yend = (permuTests_bootCI_j$rate_GA_CU)), linewidth = 1)+
      geom_segment(aes(x = loci[1]-0.05, xend = loci[1] + 0.05, y = (permuTests_bootCI_j$rate_GA_CL), yend = (permuTests_bootCI_j$rate_GA_CL)), linewidth = 1)+
      geom_segment(aes(x = loci[1]-0.05, xend = loci[1] + 0.05, y = (permuTests_bootCI_j$rate_GA_CU), yend = (permuTests_bootCI_j$rate_GA_CU)), linewidth = 1)+
      geom_point(aes(x = loci[1],  y = (permuTests_bootCI_j$rate_GA_est)), size = 1.7)+
     
      annotate (geom = "text", x = loci[1]-0.20, y = permuTests_bootCI_j$rate_GA_est, label = format(round(permuTests_bootCI_j$rate_GA_est,2),nsmall = 2), size = 4.5)+
      annotate (geom = "text", x = loci[1]-0.20, y = permuTests_bootCI_j$rate_GA_CL-
                  ifelse((permuTests_bootCI_j$rate_GA_est-permuTests_bootCI_j$rate_GA_CL)<range.y*0.035, range.y*0.035, 0), 
                label = format(round(permuTests_bootCI_j$rate_GA_CL,2),nsmall = 2), size = 4.5)+
      annotate (geom = "text", x = loci[1]-0.20, y = permuTests_bootCI_j$rate_GA_CU+
                  ifelse((permuTests_bootCI_j$rate_GA_CU-permuTests_bootCI_j$rate_GA_est)<range.y*0.035, range.y*0.035, 0), 
                label = format(round(permuTests_bootCI_j$rate_GA_CU,2),nsmall = 2), size = 4.5)+
      
      
      
      geom_segment(aes(x = loci[2], xend = loci[2], y = (permuTests_bootCI_j$rate_GB_CL), yend = (permuTests_bootCI_j$rate_GB_CU)), linewidth = 1)+
      geom_segment(aes(x = loci[2]-0.05, xend = loci[2] + 0.05, y = (permuTests_bootCI_j$rate_GB_CL), yend = (permuTests_bootCI_j$rate_GB_CL)), linewidth = 1)+
      geom_segment(aes(x = loci[2]-0.05, xend = loci[2] + 0.05, y = (permuTests_bootCI_j$rate_GB_CU), yend = (permuTests_bootCI_j$rate_GB_CU)), linewidth = 1)+
      geom_point(aes(x = loci[2],  y = (permuTests_bootCI_j$rate_GB_est)), size = 1.5)+
      
      annotate (geom = "text", x = loci[2]-0.20, y = permuTests_bootCI_j$rate_GB_est, label = format(round(permuTests_bootCI_j$rate_GB_est,2),nsmall = 2), size = 4.5)+
      annotate (geom = "text", x = loci[2]-0.20, y = permuTests_bootCI_j$rate_GB_CL-
                  ifelse((permuTests_bootCI_j$rate_GB_est-permuTests_bootCI_j$rate_GB_CL)<range.y*0.035, range.y*0.035, 0), 
                label = format(round(permuTests_bootCI_j$rate_GB_CL,2),nsmall = 2), size = 4.5)+
      annotate (geom = "text", x = loci[2]-0.20, y = permuTests_bootCI_j$rate_GB_CU+
                  ifelse((permuTests_bootCI_j$rate_GB_CU-permuTests_bootCI_j$rate_GB_est)<range.y*0.035, range.y*0.035, 0), 
                label = format(round(permuTests_bootCI_j$rate_GB_CU,2),nsmall = 2), size = 4.5)+
      
      
      geom_segment(aes(x = loci[1], xend = loci[1], y = maxy + range.y*0.05, yend = maxy + range.y*0.07))+
      geom_segment(aes(x = loci[2], xend = loci[2], y = maxy + range.y*0.05, yend = maxy + range.y*0.07))+
      geom_segment(aes(x = loci[1], xend = loci[2], y = maxy + range.y*0.07, yend = maxy + range.y*0.07))+
      annotate (geom = "text", x = (loci[1] + loci[2])/2, y = maxy + range.y*0.15, label = paste0("GMR Ratio: ", reformatRatio(permuTests_bootCI_j$rate_diff), "\nUnadj. P ",
                                                                                                 format.p2(permuTests_bootCI_j$rate_Pvalue), "; FWER P ",
                                                                                                 format.p2(permuTests_bootCI_j$rate_FWER_Pvalue)), size = 4.5)+
      annotate (geom = "text", x = loci[1], y = rangeData[1, 2], label = paste0("N = ", permuTests_bootCI_j$nA), size = 4.5)+
      annotate (geom = "text", x = loci[2], y = rangeData[1, 2], label = paste0("N = ", permuTests_bootCI_j$nB), size = 4.5)
    
    
    loci2 <- loci + 2
    permuTests_bootCI_j2 <- filter(output, Comparison== comparisons[2] & Marker == markeri)
    maxy2 <- max(permuTests_bootCI_j2$rate_GA_CU, permuTests_bootCI_j2$rate_GB_CU)
    p <- p+geom_segment(aes(x = loci2[1], xend = loci2[1], y = (permuTests_bootCI_j2$rate_GA_CL), yend = (permuTests_bootCI_j2$rate_GA_CU)), linewidth = 1)+
      geom_segment(aes(x = loci2[1]-0.05, xend = loci2[1] + 0.05, y = (permuTests_bootCI_j2$rate_GA_CL), yend = (permuTests_bootCI_j2$rate_GA_CL)), linewidth = 1)+
      geom_segment(aes(x = loci2[1]-0.05, xend = loci2[1] + 0.05, y = (permuTests_bootCI_j2$rate_GA_CU), yend = (permuTests_bootCI_j2$rate_GA_CU)), linewidth = 1)+
      geom_point(aes(x = loci2[1],  y = (permuTests_bootCI_j2$rate_GA_est)), size = 1.7)+
      annotate (geom = "text", x = loci2[1]-0.20, y = permuTests_bootCI_j2$rate_GA_est, label = format(round(permuTests_bootCI_j2$rate_GA_est,2),nsmall = 2), size = 4.5)+
      annotate (geom = "text", x = loci2[1]-0.20, y = permuTests_bootCI_j2$rate_GA_CL-
                  ifelse((permuTests_bootCI_j2$rate_GA_est-permuTests_bootCI_j2$rate_GA_CL)<range.y*0.035, range.y*0.035, 0), 
                label = format(round(permuTests_bootCI_j2$rate_GA_CL,2),nsmall = 2), size = 4.5)+
      annotate (geom = "text", x = loci2[1]-0.20, y = permuTests_bootCI_j2$rate_GA_CU+
                  ifelse((permuTests_bootCI_j2$rate_GA_CU-permuTests_bootCI_j2$rate_GA_est)<range.y*0.035, range.y*0.035, 0), 
                label = format(round(permuTests_bootCI_j2$rate_GA_CU,2),nsmall = 2), size = 4.5)+
      
      geom_segment(aes(x = loci2[2], xend = loci2[2], y = (permuTests_bootCI_j2$rate_GB_CL), yend = (permuTests_bootCI_j2$rate_GB_CU)), linewidth = 1)+
      geom_segment(aes(x = loci2[2]-0.05, xend = loci2[2] + 0.05, y = (permuTests_bootCI_j2$rate_GB_CL), yend = (permuTests_bootCI_j2$rate_GB_CL)), linewidth = 1)+
      geom_segment(aes(x = loci2[2]-0.05, xend = loci2[2] + 0.05, y = (permuTests_bootCI_j2$rate_GB_CU), yend = (permuTests_bootCI_j2$rate_GB_CU)), linewidth = 1)+
      geom_point(aes(x = loci2[2],  y = (permuTests_bootCI_j2$rate_GB_est)), size = 1.7)+
     
      annotate (geom = "text", x = loci2[2]-0.20, y = permuTests_bootCI_j2$rate_GB_est, label = format(round(permuTests_bootCI_j2$rate_GB_est,2),nsmall = 2), size = 4.5)+
      annotate (geom = "text", x = loci2[2]-0.20, y = permuTests_bootCI_j2$rate_GB_CL-
                  ifelse((permuTests_bootCI_j2$rate_GB_est-permuTests_bootCI_j2$rate_GB_CL)<range.y*0.035, range.y*0.035, 0), 
                label = format(round(permuTests_bootCI_j2$rate_GB_CL,2),nsmall = 2), size = 4.5)+
      annotate (geom = "text", x = loci2[2]-0.20, y = permuTests_bootCI_j2$rate_GB_CU+
                  ifelse((permuTests_bootCI_j2$rate_GB_CU-permuTests_bootCI_j2$rate_GB_est)<range.y*0.035, range.y*0.035, 0), 
                label = format(round(permuTests_bootCI_j2$rate_GB_CU,2),nsmall = 2), size = 4.5)+
      
      geom_segment(aes(x = loci2[1], xend = loci2[1], y = maxy2 + range.y*0.05, yend = maxy2 + range.y*0.07))+
      geom_segment(aes(x = loci2[2], xend = loci2[2], y = maxy2 + range.y*0.05, yend = maxy2 + range.y*0.07))+
      geom_segment(aes(x = loci2[1], xend = loci2[2], y = maxy2 + range.y*0.07, yend = maxy2 + range.y*0.07))+
      annotate (geom = "text", x = (loci2[1] + loci2[2])/2, y = maxy2 + range.y*0.15, label = paste0("GMR Ratio: ", reformatRatio(permuTests_bootCI_j2$rate_diff), "\nUnadj. P ",
                                                                                                 format.p2(permuTests_bootCI_j2$rate_Pvalue), "; FWER P ",
                                                                                                 format.p2(permuTests_bootCI_j2$rate_FWER_Pvalue)), size = 4.5)+
      annotate (geom = "text", x = loci2[1], y = rangeData[1, 2], label = paste0("N = ", permuTests_bootCI_j2$nA), size = 4.5)+
      annotate (geom = "text", x = loci2[2], y = rangeData[1, 2], label = paste0("N = ", permuTests_bootCI_j2$nB), size = 4.5)
    
    
    
    ggsave(filename = paste0("rate_", gsub("/","", gsub(" ", "",markeri)), figureLabel, "summary.pdf"), plot = p, path = figDir, width = 7.3, height = 6, units = "in") 
    
  }
}