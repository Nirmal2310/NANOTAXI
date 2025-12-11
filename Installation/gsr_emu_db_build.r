library(tidyverse)

set.seed(231098)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  stop("No arguments supplied. Please provide the taxa file.", call. = FALSE)
}

input_taxa_file <- args[1]

tax_data <- read.table(input_taxa_file, header=FALSE, sep="\t")

colnames(tax_data) <- c("id", "superkingdom", "phylum", "class", "order", "family", "genus", "species")

unique_taxa <- tax_data %>% distinct(species)

taxa_taxid_df <- cbind(unique_taxa, sample(10000:70000, nrow(unique_taxa), replace = FALSE))

colnames(taxa_taxid_df) <- c("species", "tax_id")

seq_tax_df <- inner_join(tax_data, taxa_taxid_df) %>% select(c("id", "tax_id"))

taxon_data <- inner_join(tax_data, taxa_taxid_df) %>% select(c("tax_id", "species", "genus", "family", "order", "class", "phylum", "superkingdom"))

write.table(seq_tax_df, "seq2tax.map.tsv", col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

write.table(taxon_data, "taxonomy.tsv", row.names = FALSE, quote = FALSE, sep = "\t")

