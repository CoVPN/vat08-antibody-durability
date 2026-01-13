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
