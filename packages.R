required = c("shiny", "shinyFiles", "markdown", "shinyBS", "validate", "tidyverse",
             "stringr", "ggpubr", "dendextend", "BiocManager", "vegan", "grid",
             "ggforce", "gridExtra", "ggsci", "scales", "viridis", "DT", "circlize",
             "ggrepel", "devtools", "compositions", "forcats","shinyjs", "plotly")
sapply(required, function(x){
  if(!require(x, character.only = TRUE)){
    install.packages(x); library(x,  character.only = TRUE)}
  else{library(x, character.only = TRUE)}
  }
)

if(!require("ComplexHeatmap", character.only = TRUE)){
  BiocManager::install("ComplexHeatmap"); library("ComplexHeatmap", character.only = TRUE)
}else{
  library("ComplexHeatmap", character.only = TRUE)
}

if(!require("pairwiseAdonis", character.only = TRUE)){
  library(devtools);install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis");library("pairwiseAdonis", character.only = TRUE)
}else {
  library("pairwiseAdonis", character.only = TRUE)
}

