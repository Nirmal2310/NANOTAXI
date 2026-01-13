required_cran = c("shiny", "shinyBS", "reticulate", "tidyverse", "shinyjs", "DT",
"plotly", "shinyFiles", "markdown", "validate", "stringr", "ggpubr",
"dendextend", "BiocManager", "vegan", "grid", "gridExtra", "ggsci",
"scales", "viridis", "circlize", "ggrepel", "devtools", "compositions",
"forcats", "formattable", "future", "promises", "ggtext", "FactoMineR",
"ggforce", "bslib")

required_bioc <- c("ComplexHeatmap", "ANCOMBC")

sapply(required_cran, function(x){
  if(!require(x, character.only = TRUE)){
    install.packages(x); suppressPackageStartupMessages(library(x,  character.only = TRUE))}
  else{suppressPackageStartupMessages(library(x, character.only = TRUE))}
  }
)

sapply(required_bioc, function(x){
  if(!require(x, character.only = TRUE)){
    BiocManager::install(x)
    suppressPackageStartupMessages(library(x, character.only = TRUE))
  }else{
    suppressPackageStartupMessages(library(x, character.only = TRUE))
  }
})

if(!require("pairwiseAdonis", character.only = TRUE)){
  suppressPackageStartupMessages(library(devtools))
  devtools::install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
  suppressPackageStartupMessages(library("pairwiseAdonis", character.only = TRUE))
} else{
  suppressPackageStartupMessages(library("pairwiseAdonis", character.only = TRUE))
}

work_dir <- getwd()

install_dir <- paste0(work_dir, "/Installation")

system(paste0('bash ', install_dir,'/realtime_install.sh'))