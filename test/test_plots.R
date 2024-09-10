library(tidyverse)

hist_data <- read.table("sample_hist.txt", header = FALSE, sep = "\t")

colnames(hist_data) <- c("Bin", "Counts")

hist_data <- hist_data %>% group_by(Bin) %>% summarise(Counts = sum(Counts)) %>%
  mutate(Selected=ifelse(Bin >= 1400 & Bin <= 1800, "Yes", "No"))


ggplot(hist_data, aes(Bin, log10(Counts), fill = Selected)) +
  geom_bar(stat = "identity", color = "black") +
  geom_text(aes(label=Counts), position=position_dodge(width=0.9), vjust=-0.25, color = "black") +
  scale_fill_manual(values = c("#008280FF","#631879FF")) +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    strip.text.x = element_text(size = 14, face = "bold"),
    axis.text.y=element_text(size=14, face = "bold"),
    axis.text.x=element_text(size=14, face = "bold"),
    axis.line = element_line(colour = "black", linewidth = 0.5, linetype = "solid" ),
    legend.title=element_blank(),
    legend.text=element_blank(),
    plot.title = element_text(hjust = 0.5,size = rel(2)),
    legend.position = "none"
  ) +
  xlab("Read Length") +
  ylab("Log10(Read Counts)")

quality_data <- read.table("sample_quality.txt", header = FALSE, sep = "\t")

colnames(quality_data) <- c("Phred")  

quality_data <- quality_data %>%
  mutate(Bin = cut(Phred, breaks = c(seq(1,max(Phred),1),max(Phred)), include.lowest = TRUE, right = TRUE, labels = as.character(c(seq(1,max(Phred),1))))) %>%
  select(-Phred) %>%
  group_by(Bin) %>% summarise(Counts = n())

colnames(quality_data) <- c("Q_Score", "Frequency")

quality_data$Q_Score <- as.numeric(quality_data$Q_Score)

quality_data <- quality_data %>% mutate(Selected = ifelse(Q_Score>=10, "Yes", "No"))

quality_data %>% ggplot(aes(Q_Score,Frequency, fill = Selected)) +
  geom_bar(stat = "identity", color = "black") +
  theme_minimal() +
  xlab("Phred Score") +
  ylab("Read Counts") +
  theme(
    axis.title.x = element_text(size = 14, face = "bold", colour = "black"),
    axis.title.y = element_text(size = 14, face = "bold", colour = "black"),
    strip.text.x = element_text(size = 14, face = "bold", colour = "black"),
    axis.text.y= element_text(size=14, face = "bold", colour = "black"),
    axis.text.x= element_text(size=14, face = "bold", colour = "black"),
    axis.line = element_line(colour = "black", linewidth = 0.5, linetype = "solid" ),
    strip.background = element_blank(),
    legend.position = "none"
  ) +
  scale_fill_manual(values = c("#631879FF", "#008280FF"))

taxa_counts_df <- read.table("barcode01_final_kraken2_result.txt", header = FALSE, sep = "\t")

taxa_counts_df[taxa_counts_df==""] <- "Unclassified"

colnames(taxa_counts_df) <- c("Tax_ID", "Counts", "Superkingdom", "Phylum",
                              "Class", "Order", "Family", "Genus", "Species")

taxa <- "Species"

taxa_counts_df %>% dplyr::select(c(!!sym(taxa), Counts)) %>% 
  dplyr::filter(!!sym(taxa) !="Unclassified") %>% dplyr::summarise(Counts = sum(Counts))

taxa_counts_df <- taxa_counts_df %>% group_by(Species) %>% summarise(Counts = sum(Counts))

taxa_counts_df <- taxa_counts_df %>% mutate(Abundance = (Counts/sum(Counts))*100) %>% filter(Abundance >=0.5)

ggplot(taxa_counts_df, aes(Species, Counts)) +
  geom_bar(stat = "identity", fill = "#631879FF", color = "black" ) +
  theme_minimal() +
  xlab("Species") +
  ylab("Read Counts") +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 14, face = "bold", colour = "black"),
    axis.line = element_line(colour = "black", linewidth = 0.5, linetype = "solid" ),
    strip.text.x = element_text(size = 14, face = "bold", colour = "black"),
    axis.text.y= element_text(size=10, face = "bold", colour = "black"),
    axis.text.x= element_text(size=10, face = "bold", colour = "black", angle = 45, vjust = 1, hjust = 1),
    legend.position = "none",
    title = element_text(size = 10, face = "bold")
  ) +
  ggtitle("Species with atleast 0.5% Relative Abundance")



sample_list <- gsub("_final_kraken2_result.txt", "", list.files(getwd())[grep("\\_final_kraken2_result.txt$", 
                                                       list.files(getwd()))])

classification_data_list <- list()

quality_data_list <- list()

hist_data_list <- list()

for (i in 1:length(sample_list))
{
  classification_data_list[[i]] <- read.delim(file = paste0(sample_list[i],"_final_kraken2_result.txt"), header = FALSE, sep = "\t")
  
  colnames(classification_data_list[[i]]) <- c("TAX_ID", "Counts", "Superkingdom", "Phylum",
                                       "Class", "Order", "Family", "Genus", "Species")
  
  classification_data_list[[i]][classification_data_list[[i]]==""] <- "Unclassified"
  
  quality_data_list[[i]] <- read.delim(file = paste0(sample_list[i],"_quality.txt"), header = FALSE, sep = "\t")
  
  colnames(quality_data_list[[i]]) <- "Phred"
  
  hist_data_list[[i]] <- read.delim(file = paste0(sample_list[i],"_hist.txt"), sep = "\t", header = FALSE)
  
  colnames(hist_data_list[[i]]) <- c("Bin", "Frequency")
}

taxa <- "Phylum"

df <- classification_data_list[[1]]

df <- df %>% dplyr::select(c(!!sym(taxa), Counts)) %>% filter(!!sym(taxa)!="Unclassified") %>% group_by(!!sym(taxa)) %>% summarise(Counts = sum(Counts)) %>% mutate(Abundance = (Counts/sum(Counts))*100)

df <- df %>% filter(Abundance>=0.5)

df <- classification_data_list[[1]]

df %>% group_by(Species) %>% mutate(Abundance = (Counts/sum(Counts))*100) %>% filter(Abundance > 1)

ggplot(df, aes(!!(sym(taxa)), Counts)) +
  geom_bar(stat = "identity", fill = "#631879FF", color = "black" ) +
  theme_minimal() +
  ylab("Read Counts") +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 14, face = "bold", colour = "black"),
    axis.line = element_line(colour = "black", linewidth = 0.5, linetype = "solid" ),
    strip.text.x = element_text(size = 14, face = "bold", colour = "black"),
    axis.text.y= element_text(size=10, face = "bold", colour = "black"),
    axis.text.x= element_text(size=10, face = "bold", colour = "black", angle = 90, vjust = 0.5, hjust=1),
    legend.position = "none",
    title = element_text(size = 10, face = "bold")
  ) +
  ggtitle(paste0(taxa," with atleast 0.5% Relative Abundance"))

library(httr)

library(jsonlite)

res = GET("https://www.bv-brc.org/api/genome")

fromJSON(rawToChar(res$content)) %>% as.data.frame() %>% View()
