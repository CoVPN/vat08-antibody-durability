library(GGally)
library(scales)
library(dplyr)
scientific_10 <- function(x) {
  y <- NULL
  for(i in 1:length(x)){
    y[i] <- parse(text=gsub("1e\\+*", " 10^", scales::scientific_format()(x[i])))
  }
  return(y)
}
mycorrelations <- function(data,mapping,...){
  data2 = data
  data2$x = as.numeric(data[,as_label(mapping$x)])
  data2$y = as.numeric(data[,as_label(mapping$y)])
  data2$group = data[,as_label(mapping$colour)]
  
  
  g <- NULL
  estimate <- NULL
  i = 1
  for(gr in unique(data2$group)){
    g[i] <- gr
    df <- filter(data2, group == gr)
    estimate[i] = round(as.numeric(cor(df$x,df$y,method="pearson", use = "complete.obs")),2)
    i <- i + 1
  }
  correlation_df <- data.frame("group" = g, "estimate" = estimate)
  # correlation_df = data2 %>% 
  #   #bind_rows(data2 %>% mutate(group="Overall Corr")) %>%
  #   dplyr::group_by(group) %>% 
  #   #filter(sum(!is.na(x),na.rm=T)>1) %>%
  #   #filter(sum(!is.na(y),na.rm=T)>1) %>%
  #   summarize(estimate = round(as.numeric(cor(x,y,method="spearman", use = "complete.obs")),2)
  #             #,
  #             #pvalue = cor.test(x,y,method="spearman", use = "complete.obs")$p.value,
  #             #pvalue_star = as.character(symnum(pvalue, corr = FALSE, na = FALSE, 
  #             #                                  cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), 
  #             #                                  symbols = c("***", "**", "*", "'", " ")))
  #   )%>%group_by(group) 
  #%>% mutate(group = factor(group, levels=c(as.character(unique(sort(data[,as_label(mapping$colour)]))), "Overall Corr")))
  correlation_df$group <- factor(correlation_df$group, levels = c("Naïve Placebo", "Naïve Vaccine", "Non-naïve Placebo", "Non-naïve Vaccine"))
  
  # ggplot(data=correlation_df, aes(x=1,y=group,color=group))+geom_text(aes(label=paste0(group,": ",estimate,pvalue_star)), size = 3)
  ggplot(data=correlation_df, aes(x=1,y=group,color=group))+
    geom_text(aes(label=paste0(group,": ",ifelse(!is.na(estimate), estimate, "NE"))), size = 3)+
    scale_y_discrete(expand = expansion(mult = c(0.6, 0.6)))+
    theme(panel.grid = element_blank()           # remove all gridlines
          # panel.background = element_blank(),     # remove plot panel background
          # axis.line = element_line(color = "black")  # keep black axis/frame lines
    )
}


scatterplot_lower_titer <- function(data, mapping, ...) {
  ggplot(data = data, mapping = mapping, ...) +
    geom_point(...) +
    #geom_smooth(method = "loess", se = FALSE, linewidth = 0.5,...)+
    scale_y_continuous(limits = c(1.5, 5.5), breaks = c(2, 3, 4, 5), labels = scientific_10(10^c(2, 3, 4, 5)), minor_breaks = NULL) + # Set desired y-axis limits
    scale_x_continuous(limits = c(1.5, 5.5), breaks = c(2, 3, 4, 5), labels = scientific_10(10^c(2, 3, 4, 5)), minor_breaks = NULL)   # Set desired x-axis limits
}

density_diag_titer <- function(data, mapping, ...) {
  ggplot(data = data, mapping = mapping, ...) +
    geom_density(...) +
    #scale_y_continuous(limits = c(1.5, 5.5), breaks = c(2, 3, 4, 5)) + # Set desired y-axis limits
    scale_x_continuous(limits = c(1.5, 5.5), breaks = c(2, 3, 4, 5), labels = scientific_10(10^c(2, 3, 4, 5)), minor_breaks = NULL)  # Set desired x-axis limits
}

scatterplot_lower_IgG <- function(data, mapping, ...) {
  ggplot(data = data, mapping = mapping, ...) +
    geom_point(...) +
    #geom_smooth(method = "loess", se = FALSE,linewidth = 0.5, ...)+
    scale_y_continuous(limits = c(2.3, 6.2), breaks = c(2, 3, 4, 5,6), labels = scientific_10(10^c(2, 3, 4, 5, 6)), minor_breaks = NULL) + # Set desired y-axis limits
    scale_x_continuous(limits = c(2.3, 6.2), breaks = c(2, 3, 4, 5,6), labels = scientific_10(10^c(2, 3, 4, 5, 6)), minor_breaks = NULL)   # Set desired x-axis limits
}

density_diag_IgG <- function(data, mapping, ...) {
  ggplot(data = data, mapping = mapping, ...) +
    geom_density(...) +
    #scale_y_continuous(limits = c(1.5, 5.5), breaks = c(2, 3, 4, 5)) + # Set desired y-axis limits
    scale_x_continuous(limits = c(2.3, 6.2), breaks = c(2, 3, 4, 5,6), labels = scientific_10(10^c(2, 3, 4, 5, 6)), minor_breaks = NULL)  # Set desired x-axis limits
}