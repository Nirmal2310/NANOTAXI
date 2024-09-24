required = c("shiny", "shinyBS", "reticulate", "tidyverse", "shinyjs", "DT",
"plotly", "shinyFiles", "markdown", "validate", "stringr", "ggpubr",
"dendextend", "BiocManager", "vegan", "grid", "gridExtra", "ggsci",
"scales", "viridis", "circlize", "ggrepel", "devtools", "compositions", "forcats", "formattable")

sapply(required, function(x){
  if(!require(x, character.only = TRUE)){
    install.packages(x); suppressPackageStartupMessages(library(x,  character.only = TRUE))}
  else{suppressPackageStartupMessages(library(x, character.only = TRUE))}
  }
)

if(!require("ComplexHeatmap", character.only = TRUE)){
  BiocManager::install("ComplexHeatmap")
  suppressPackageStartupMessages(library("ComplexHeatmap", character.only = TRUE))
}else{
  suppressPackageStartupMessages(library("ComplexHeatmap", character.only = TRUE))
}

if(!require("pairwiseAdonis", character.only = TRUE)){
  suppressPackageStartupMessages(library(devtools))
  install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
  suppressPackageStartupMessages(library("pairwiseAdonis", character.only = TRUE))
} else{
  suppressPackageStartupMessages(library("pairwiseAdonis", character.only = TRUE))
}

work_dir <- getwd()

install_dir <- paste0(work_dir, "/Installation")

system(paste0('bash ', install_dir,'/kraken_ncbi_install.sh'))
