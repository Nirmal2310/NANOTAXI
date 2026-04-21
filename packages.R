required_cran = c("shiny", "shinyBS", "reticulate", "tidyverse", "shinyjs", "DT",
"plotly", "shinyFiles", "markdown", "validate", "stringr", "ggpubr",
"dendextend", "BiocManager", "vegan", "grid", "gridExtra", "ggsci",
"scales", "viridis", "circlize", "ggrepel", "devtools", "compositions",
"forcats", "formattable", "future", "promises", "ggtext", "FactoMineR",
"ggforce", "bslib", "cowplot")

required_bioc <- c("ComplexHeatmap", "ANCOMBC", "MicrobiomeProfiler", "enrichplot")

sapply(required_cran, function(x){
  if(!require(x, character.only = TRUE)){
    install.packages(x, ask=FALSE); suppressMessages(suppressPackageStartupMessages(library(x,  character.only = TRUE)))}
  else{suppressMessages(suppressPackageStartupMessages(library(x, character.only = TRUE)))}
  }
)

sapply(required_bioc, function(x){
  if(!require(x, character.only = TRUE)){
    BiocManager::install(x, ask=FALSE)
    suppressMessages(suppressPackageStartupMessages(library(x, character.only = TRUE)))
  }else{
    suppressMessages(suppressPackageStartupMessages(library(x, character.only = TRUE)))
  }
})

if(!require("pairwiseAdonis", character.only = TRUE)){
  suppressMessages(suppressPackageStartupMessages(library(devtools)))
  devtools::install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis", upgrade = FALSE)
  suppressMessages(suppressPackageStartupMessages(library("pairwiseAdonis", character.only = TRUE)))
} else{
  suppressMessages(suppressPackageStartupMessages(library("pairwiseAdonis", character.only = TRUE)))
}

work_dir <- getwd()

install_dir <- paste0(work_dir, "/Installation")

system(paste0('bash ', install_dir,'/realtime_install.sh'))