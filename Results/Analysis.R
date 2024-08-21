library(tidyverse)
library(stringr)
library(ggpubr)
library(dendextend)
library(ComplexHeatmap)
library(vegan)
library(ggrepel)
library(grid)
library(gridExtra)
library(ggforce)
library(ggsci)
library(scales)
library(ComplexHeatmap)
library(dendextend)
library(viridis)
library(circlize)
library(FactoMineR)

######################## EMU RESULTS ##############################################################################################################

emu_dir <- paste0(getwd(), "/EMU_results/")

file_list <- gsub(".txt", "", list.files(emu_dir)[grep("\\_final_emu_result.txt$",
                                                       list.files(emu_dir))])

samples_header <- gsub("_final_emu_result","",file_list)

sample_metadata <- read.table(file = "Sample_Information.csv", sep = ",", 
                              header = TRUE)

sample_metadata$Group <- as.factor(sample_metadata$Group)

sample_data_list <- list()

abundance_data_list <- list()

rel_abundance_data_list <- list()

lineage <- "Genus"

for (i in 1:length(file_list))
{
  sample_data_list[[i]] <- read.delim(file = paste0(emu_dir,"/",file_list[i],".txt"), header = FALSE)
  
  colnames(sample_data_list[[i]]) <- c("TAX_ID", "Counts", "Kingdom",
                                       "Phylum", "Class", "Order", "Family", "Genus",
                                       "Species")
  
  sample_data_list[[i]][sample_data_list[[i]]==""] <- "Unclassified"
  
  col_num <- which(colnames(sample_data_list[[i]])==lineage)
  
  abundance_data_list[[i]] <- sample_data_list[[i]] %>% group_by(sample_data_list[[i]][,col_num]) %>% 
    summarise(Counts = sum(Counts))
  
  rel_abundance_data_list[[i]] <- abundance_data_list[[i]] %>% mutate(Freq = (Counts/sum(Counts))*100)
  
  rel_abundance_data_list[[i]] <- rel_abundance_data_list[[i]][,-2]
  
  colnames(rel_abundance_data_list[[i]]) <- c(lineage, samples_header[i])
  
  colnames(abundance_data_list[[i]]) <- c(lineage, samples_header[i])
  
  
}

  rel_abundance_data <- rel_abundance_data_list %>% purrr::reduce(full_join, by=lineage)
  
  rel_abundance_data <- as.data.frame(rel_abundance_data)
  
  rel_abundance_data[is.na(rel_abundance_data)] <- 0
  
  abundance_data <- abundance_data_list %>% purrr::reduce(full_join, by="Species")
  
  abundance_data <- as.data.frame(abundance_data)
  
  abundance_data[is.na(abundance_data)] <- 0
  
  rel_abundance_matrix <- as.matrix(rel_abundance_data[,-1])
  
  rownames(rel_abundance_matrix) <- rel_abundance_data[,which(colnames(rel_abundance_data)==lineage)]
  
  abundance_matrix <- abundance_data[,-1]
  
  rownames(abundance_matrix) <- abundance_data[,which(colnames(abundance_data)==lineage)]

############################## Alpha Diversity #################################

stacked_df <- rel_abundance_data

stacked_df[,which(colnames(stacked_df)==lineage)] <- as.character(stacked_df[,which(colnames(stacked_df)==lineage)])

stacked_df <- stacked_df %>% 
  pivot_longer(cols = -(which(colnames(stacked_df)==lineage)), names_to = "Sample_Id", values_to = "Abundance") %>% 
  filter(Abundance > 0) %>% 
  mutate(Name = ifelse(Abundance>2, stacked_df[,which(colnames(stacked_df)==lineage)], "< 2% Abundance"))

unique_names <- unique(stacked_df$Name)

fill_colors <- c("< 2% Abundance" = "black", setNames(viridis_pal(option = "D")(length(unique_names)-1), unique_names[
  unique_names != "< 2% Abundance"]))

ggplot(stacked_df, aes(x = Sample_Id, y = Abundance, fill = Name)) +
  geom_bar(stat = "identity", position = "stack") + 
  coord_flip() +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 10, face = "bold", angle = 45, hjust = 1),
    axis.text.y = element_text(size = 10, face = "bold"),
    axis.title.x = element_text(size = 10, face = "bold"),
    axis.title.y = element_text(size = 10, face = "bold"),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(face = "bold"),
    legend.position = "bottom"
  ) +
  guides(fill = guide_legend(ncol = 10)) +
  labs(x = "BARCODE", y="% RELATIVE ABUNDANCE", fill = lineage) + 
  scale_fill_manual(values = fill_colors) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 105))

sample_annotation <- data.frame(sample_metadata$Group)

colnames(sample_annotation) <- "Group"

sample_annotation$Group <- factor(sample_annotation$Group)

alpha_diversity_data <- rel_abundance_data

alpha_diversity_data[,which(colnames(alpha_diversity_data)==lineage)] <- as.character(
  alpha_diversity_data[,which(colnames(alpha_diversity_data)==lineage)])

alpha_diversity_data <- alpha_diversity_data %>% pivot_longer(cols = -(which(colnames(alpha_diversity_data)==lineage)), names_to = "Sample_Id", values_to = "Abundance")


alpha_diversity_data <- alpha_diversity_data %>%
  group_by(Sample_Id) %>% 
  summarise(Shannon = diversity(Abundance, index = "shannon"),
            Simpson = diversity(Abundance, index = "simpson"))


alpha_diversity_data <- inner_join(alpha_diversity_data, sample_metadata)

alpha_diversity_data$Group <- factor(alpha_diversity_data$Group)

alpha_diversity_data <- pivot_longer(alpha_diversity_data, cols = c(Shannon, Simpson), 
                                     names_to = "Diversity", values_to = "Value")

alpha_div_p <- compare_means(Value~Group, data = alpha_diversity_data, method = "wilcox",
                             p.adjust.method = "holm", group.by = "Diversity")

shannon_plot <- alpha_diversity_data %>% filter(Diversity == "Shannon") %>%
  ggplot(aes(x=Group, y=Value, color=Group)) +
  geom_boxplot() +
  scale_color_manual(values = pal_aaas("default")(length(levels(sample_annotation$Group)))) +
  theme_classic() +
  labs(y= "Alpha Diversity", x = "") +
  stat_pvalue_manual(subset(alpha_div_p, Diversity=="Shannon"), label = "p.signif", y.position = max(
    subset(alpha_diversity_data, Diversity=="Shannon")$Value) + 0.1, hide.ns = "p", step.increase = 0.1,
    tip.length = 0.02, bracket.size = 0.8, size = 8) +
  guides(color = guide_legend(title = "Shannon", title.position = "top")) +
  theme(
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    strip.text.x = element_text(size = 14, face = "bold"),
    axis.text.y=element_text(size=14, face = "bold"),
    axis.text.x= element_blank(),
    axis.ticks.x = element_blank(),
    legend.title=element_text(colour="black",size=10, face = "bold"),
    legend.text=element_text(colour="black", size=10, face = "bold"),
    plot.title = element_text(hjust = 0.5,size = rel(2))
  )

simpson_plot <- alpha_diversity_data %>% filter(Diversity == "Simpson") %>%
  ggplot(aes(x=Group, y=Value, color=Group)) +
  geom_boxplot() +
  scale_color_manual(values = pal_aaas("default")(length(levels(sample_annotation$Group)))) + 
  theme_classic() +
  labs(y= "Alpha diversity", x = "") +
  stat_pvalue_manual(subset(alpha_div_p, Diversity=="Simpson"), y.position = max(
    subset(alpha_diversity_data, Diversity=="Simpson")$Value) + 0.01, label = "p.signif",
    hide.ns = "p", step.increase = 0.1, tip.length = 0.02, bracket.size = 0.8, size = 8) +
  guides(color = guide_legend(title = "Simpson", title.position = "top")) +
  theme(
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    strip.text.x = element_text(size = 14, face = "bold"),
    axis.text.y=element_text(size=14, face = "bold"),
    axis.text.x= element_blank(),
    axis.ticks.x = element_blank(),
    legend.title=element_text(colour="black",size=10, face = "bold"),
    legend.text=element_text(colour="black", size=10, face = "bold"),
    plot.title = element_text(hjust = 0.5,size = rel(2))
  )


as_ggplot(grid.grabExpr(grid.arrange(shannon_plot, simpson_plot, ncol=2)))

###################### PCA Plot ################################################

pca <- prcomp(t(rel_abundance_matrix), scale = TRUE)

princomp(t(rel_abundance_matrix), cor=TRUE)

pca.var <- pca$sdev^2

pca.var.per <- round(pca.var/sum(pca.var)*100, 1)

pca_df <- data.frame(Sample_Id = rownames(pca$x), X=pca$x[,1], Y=pca$x[,2])

pca_df <- inner_join(pca_df, sample_metadata)

pca_df$Group <- factor(pca_df$Group)

ggplot(data = pca_df, aes(x = X, y = Y, color = Group)) +
  geom_point(size=5) +
  xlab(paste0("PC1 (", pca.var.per[1], "%", ")")) +
  ylab(paste0("PC2 (", pca.var.per[2], "%", ")")) +
  scale_color_manual(values = pal_aaas("default")(length(levels(sample_metadata$Group)))) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 15, face = "bold"),
    axis.text.y = element_text(size = 15, face = "bold"),
    legend.text = element_text(size = 15, face = "bold"),
    axis.title.x = element_text(size = 15, face = "bold"),
    axis.title.y = element_text(size = 15, face = "bold"),
    legend.title = element_text(size = 15, face = "bold")
  ) +
  geom_vline(xintercept = 0, color = "black", linetype = "dashed", linewidth = 1) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed", linewidth = 1)


nmds_dist <- metaMDS(t(abundance_matrix), distance = "bray", trymax = 100)

nmds_df <- as.data.frame(nmds_dist$points)

nmds_df$Group <- sample_metadata$Group

ggplot(data = nmds_df, aes(x = MDS1, y = MDS2, color = Group)) +
  geom_point(size=5) +
  scale_color_manual(values = pal_aaas("default")(length(levels(sample_annotation$Group)))) +
  theme_classic() +
  geom_vline(xintercept = 0, color = "black", linetype = "dashed", linewidth = 1) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed", linewidth = 1) +
  theme(
    axis.text.x = element_text(size = 25, face = "bold"),
    axis.text.y = element_text(size = 25, face = "bold"),
    legend.text = element_text(size = 25, face = "bold"),
    axis.title.x = element_text(size = 25, face = "bold"),
    axis.title.y = element_text(size = 25, face = "bold"),
    legend.title = element_text(size = 25, face = "bold")
  )

########################### PCoA Plot ###########################################

pcoa_df <- wcmdscale(vegdist(t(abundance_matrix), method = "aitchison", pseudocount=1), k=2) %>% as.data.frame()

pcoa_df$Group <- sample_metadata$Group

pcoa_df$Group <- factor(pcoa_df$Group)

colnames(pcoa_df) <- c("Axis.1", "Axis.2", "Group")

pal_aaas("default")(length(levels(pcoa_df$Group)))

ggplot(data = pcoa_df, aes(x = Axis.1, y = Axis.2, color = Group)) +
  geom_point(size=5) +
  scale_color_manual(values = pal_aaas("default")(length(levels(pcoa_df$Group)))) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 25, face = "bold"),
    axis.text.y = element_text(size = 25, face = "bold"),
    legend.text = element_text(size = 25, face = "bold"),
    axis.title.x = element_text(size = 25, face = "bold"),
    axis.title.y = element_text(size = 25, face = "bold"),
    legend.title = element_text(size = 25, face = "bold")
  ) +
  geom_vline(xintercept = 0, color = "black", linetype = "dashed", linewidth = 1) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed", linewidth = 1)
  

############################ Heatmap ###########################################

heatmap_df <- as.data.frame(log10(rel_abundance_matrix+0.00001))

sample_annotation <- data.frame(sample_metadata$Group)

colnames(sample_annotation) <- "Group"

sample_annotation$Group <- factor(sample_annotation$Group)

col_list <- list(Group = setNames(pal_aaas("default")(length(levels(sample_metadata$Group))) ,levels(sample_metadata$Group)))

row_annotate <- rowAnnotation(
  df = sample_annotation,
  col = col_list, show_annotation_name = FALSE)

row_dend <-  hclust(dist(t(heatmap_df)), method = "complete")

column_dend <- hclust(dist(heatmap_df), method = "complete")

tmp <- heatmap_df %>% pivot_longer(cols = 1:length(colnames(heatmap_df)), names_to = "Sample",
                                   values_to = "Abundance")

col_fun <- colorRamp2(c(0, as.vector(quantile(tmp$Abundance)[4]), 30), c("#181C7D", "white", "#AC373D"))

heatmap_plot <- Heatmap(t(heatmap_df), name = "Relative Abundance",
                        row_names_gp = gpar(fontsize = 10, fontface = "bold"),
                        cluster_rows = color_branches(row_dend),
                        cluster_columns = color_branches(column_dend),
                        show_column_names = FALSE,
                        show_row_names = TRUE,
                        right_annotation = row_annotate,
                        col = col_fun(seq(0,30)))

heatmap_ggplot <- as_ggplot(grid.grabExpr(print(heatmap_plot)))

print(heatmap_ggplot)


################################### Kraken2 Data ###############################

kraken_dir <- paste0(getwd(),"/Kraken2_Results")

file_list <- gsub(".txt", "", list.files(kraken_dir)[grep("\\_final_kraken2_result.txt$",
                                                       list.files(kraken_dir))])

samples_header <- gsub("_final_kraken2_result","",file_list)

sample_metadata <- read.table(file = paste0(kraken_dir,"/Sample_Information.csv"), sep = ",", 
                              header = TRUE)

sample_metadata$Group <- as.factor(sample_metadata$Group)

sample_data_list <- list()

abundance_data_list <- list()

rel_abundance_data_list <- list()

lineage <- "Genus"

for (i in 1:length(file_list))
{
  sample_data_list[[i]] <- read.delim(file = paste0(kraken_dir,"/",file_list[i],".txt"), header = FALSE)
  
  colnames(sample_data_list[[i]]) <- c("TAX_ID", "Counts", "Kingdom",
                                       "Phylum", "Class", "Order", "Family", "Genus",
                                       "Species")
  
  sample_data_list[[i]][sample_data_list[[i]]==""] <- "Unclassified"
  
  col_num <- which(colnames(sample_data_list[[i]])==lineage)
  
  abundance_data_list[[i]] <- sample_data_list[[i]] %>% group_by(sample_data_list[[i]][,col_num]) %>% 
    summarise(Counts = sum(Counts))
  
  rel_abundance_data_list[[i]] <- abundance_data_list[[i]] %>% mutate(Freq = (Counts/sum(Counts))*100)
  
  rel_abundance_data_list[[i]] <- rel_abundance_data_list[[i]][,-2]
  
  colnames(rel_abundance_data_list[[i]]) <- c(lineage, samples_header[i])
  
  colnames(abundance_data_list[[i]]) <- c(lineage, samples_header[i])
  
  
}

rel_abundance_data <- rel_abundance_data_list %>% purrr::reduce(full_join, by=lineage)

rel_abundance_data <- as.data.frame(rel_abundance_data)

rel_abundance_data[is.na(rel_abundance_data)] <- 0

abundance_data <- abundance_data_list %>% purrr::reduce(full_join, by=lineage)

abundance_data <- as.data.frame(abundance_data)

abundance_data[is.na(abundance_data)] <- 0

rel_abundance_matrix <- as.matrix(rel_abundance_data[,-1])

rownames(rel_abundance_matrix) <- rel_abundance_data[,which(colnames(rel_abundance_data)==lineage)]

abundance_matrix <- abundance_data[,-1]

rownames(abundance_matrix) <- abundance_data[,which(colnames(abundance_data)==lineage)]

rel_abundance_matrix <- rel_abundance_matrix[apply(rel_abundance_matrix, 1, function(row) any(mean(row) >0.1 )), ]

abundance_matrix <- abundance_matrix[which(rownames(abundance_matrix) %in% rownames(rel_abundance_matrix)),]

rel_abundance_data <- rel_abundance_data[which(rel_abundance_data[,1] %in% rownames(rel_abundance_matrix)),]

############################## Alpha Diversity #################################

stacked_df <- rel_abundance_data

required_col <- which(colnames(stacked_df)==lineage)

stacked_df[,required_col] <- as.character(stacked_df[,required_col])

stacked_df <- stacked_df %>% 
  pivot_longer(cols = -(required_col), names_to = "Sample_Id", values_to = "Abundance") %>% 
  filter(Abundance > 0)

stacked_df <- stacked_df %>% group_by(Sample_Id) %>%
  arrange(desc(Abundance)) %>% slice(1:ifelse(n()<10,n(),10))

stacked_df$Sample_Id <- gsub("barcode","", stacked_df$Sample_Id)

temp <- stacked_df %>% group_by(Sample_Id) %>% summarise(Abundance = 100-sum(Abundance)) %>% mutate(Species="Others") %>% 
  dplyr::select(Species, Sample_Id, Abundance)

colnames(temp) <- colnames(stacked_df)

stacked_df <- rbind(stacked_df, temp) %>% arrange(Sample_Id)

stacked_df[[required_col]] <- factor(stacked_df[[required_col]])

stacked_df[[required_col]] <- factor(stacked_df[[required_col]], 
                                     levels = c("Others", setdiff(levels(stacked_df[[required_col]]), "Others")))

fill_colors <- c("Others" = "#D3D3D3",setNames(viridis_pal(option = "D")(
  length(levels(stacked_df[[required_col]]))-1), 
  levels(stacked_df[[required_col]])[levels(stacked_df[[required_col]]) != "Others"]))

ggplot(stacked_df, aes(x = Sample_Id, y = Abundance, fill = stacked_df[[required_col]])) +
  geom_bar(stat = "identity", position = "stack", color = "black", linewidth = 0.2) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 10, face = "bold", angle = 45, hjust = 1),
    axis.text.y = element_text(size = 10, face = "bold"),
    axis.title.x = element_text(size = 10, face = "bold"),
    axis.title.y = element_text(size = 10, face = "bold"),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(face = "bold"),
    legend.position = "bottom",
    title = element_text(size = 10, face = "bold")
  ) +
  guides(fill = guide_legend(ncol = 5)) +
  labs(x = "BARCODE", y="% RELATIVE ABUNDANCE", fill = lineage) + 
  scale_fill_manual(values = fill_colors) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,102)) +
  ggtitle(label = paste0("Stacked Barplot Showing Top 10 ",lineage," Across All Samples"))

sample_annotation <- data.frame(sample_metadata$Group)

colnames(sample_annotation) <- "Group"

sample_annotation$Group <- factor(sample_annotation$Group)

alpha_diversity_data <- rel_abundance_data

alpha_diversity_data[,which(colnames(alpha_diversity_data)==lineage)] <- as.character(
  alpha_diversity_data[,which(colnames(alpha_diversity_data)==lineage)])

alpha_diversity_data <- alpha_diversity_data %>% pivot_longer(cols = -(which(colnames(alpha_diversity_data)==lineage)), names_to = "Sample_Id", values_to = "Abundance")


alpha_diversity_data <- alpha_diversity_data %>%
  group_by(Sample_Id) %>% 
  summarise(Shannon = diversity(Abundance, index = "shannon"),
            Simpson = diversity(Abundance, index = "simpson"))


alpha_diversity_data <- inner_join(alpha_diversity_data, sample_metadata)

alpha_diversity_data$Group <- factor(alpha_diversity_data$Group)

alpha_diversity_data <- pivot_longer(alpha_diversity_data, cols = c(Shannon, Simpson), 
                                     names_to = "Diversity", values_to = "Value")

alpha_div_p <- compare_means(Value~Group, data = alpha_diversity_data, method = "wilcox",
                             p.adjust.method = "holm", group.by = "Diversity")

shannon_plot <- alpha_diversity_data %>% filter(Diversity == "Shannon") %>%
  ggplot(aes(x=Group, y=Value, color=Group)) +
  geom_boxplot() +
  scale_color_manual(values = pal_aaas("default")(length(levels(sample_annotation$Group)))) +
  theme_classic() +
  labs(y= "Alpha Diversity", x = "") +
  stat_pvalue_manual(subset(alpha_div_p, Diversity=="Shannon"), label = "p.signif", y.position = max(
    subset(alpha_diversity_data, Diversity=="Shannon")$Value) + 0.1, hide.ns = "p", step.increase = 0.1,
    tip.length = 0.02, bracket.size = 0.8, size = 8) +
  guides(color = guide_legend(title = "Shannon", title.position = "top")) +
  theme(
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    strip.text.x = element_text(size = 14, face = "bold"),
    axis.text.y=element_text(size=14, face = "bold"),
    axis.text.x= element_blank(),
    axis.ticks.x = element_blank(),
    legend.title=element_text(colour="black",size=10, face = "bold"),
    legend.text=element_text(colour="black", size=10, face = "bold"),
    plot.title = element_text(hjust = 0.5,size = rel(2))
  )

simpson_plot <- alpha_diversity_data %>% filter(Diversity == "Simpson") %>%
  ggplot(aes(x=Group, y=Value, color=Group)) +
  geom_boxplot() +
  scale_color_manual(values = pal_aaas("default")(length(levels(sample_annotation$Group)))) + 
  theme_classic() +
  labs(y= "Alpha diversity", x = "") +
  stat_pvalue_manual(subset(alpha_div_p, Diversity=="Simpson"), y.position = max(
    subset(alpha_diversity_data, Diversity=="Simpson")$Value) + 0.01, label = "p.signif",
    hide.ns = "p", step.increase = 0.1, tip.length = 0.02, bracket.size = 0.8, size = 8) +
  guides(color = guide_legend(title = "Simpson", title.position = "top")) +
  theme(
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    strip.text.x = element_text(size = 14, face = "bold"),
    axis.text.y=element_text(size=14, face = "bold"),
    axis.text.x= element_blank(),
    axis.ticks.x = element_blank(),
    legend.title=element_text(colour="black",size=10, face = "bold"),
    legend.text=element_text(colour="black", size=10, face = "bold"),
    plot.title = element_text(hjust = 0.5,size = rel(2))
  )


as_ggplot(grid.grabExpr(grid.arrange(shannon_plot, simpson_plot, ncol=2)))

###################### PCA Plot ################################################

pca <- PCA(rel_abundance_matrix, ncp = 2, graph = FALSE)

pca_df <- pca$var$coord %>% as.data.frame()

pca_df$Sample_Id <- rownames(pca_df)

pca_df <- inner_join(pca_df, sample_metadata)

pca_df %>% ggplot(aes(x=Dim.1, y=Dim.2, color = Group)) + geom_point(size = 5) + xlab(paste0("PC1 (",round(test$eig[1,3], 2),")")) + ylab(paste0("PC2 (",round(test$eig[2,3], 2),")")) + scale_color_manual(values = pal_aaas("default")(length(levels(sample_annotation$Group)))) +
  theme_classic() +
  geom_vline(xintercept = 0, color = "black", linetype = "dashed", linewidth = 1) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed", linewidth = 1) +
  theme(
    axis.text.x = element_text(size = 25, face = "bold"),
    axis.text.y = element_text(size = 25, face = "bold"),
    legend.text = element_text(size = 25, face = "bold"),
    axis.title.x = element_text(size = 25, face = "bold"),
    axis.title.y = element_text(size = 25, face = "bold"),
    legend.title = element_text(size = 25, face = "bold")
  )


nmds_dist <- metaMDS(t(abundance_matrix), distance = "bray", trymax = 100)

nmds_df <- as.data.frame(nmds_dist$points)

nmds_df$Group <- sample_metadata$Group

ggplot(data = nmds_df, aes(x = MDS1, y = MDS2, color = Group)) +
  geom_point(size=5) +
  scale_color_manual(values = pal_aaas("default")(length(levels(sample_annotation$Group)))) +
  theme_classic() +
  geom_vline(xintercept = 0, color = "black", linetype = "dashed", linewidth = 1) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed", linewidth = 1) +
  theme(
    axis.text.x = element_text(size = 25, face = "bold"),
    axis.text.y = element_text(size = 25, face = "bold"),
    legend.text = element_text(size = 25, face = "bold"),
    axis.title.x = element_text(size = 25, face = "bold"),
    axis.title.y = element_text(size = 25, face = "bold"),
    legend.title = element_text(size = 25, face = "bold")
  )

########################### PCoA Plot ###########################################

pcoa_df <- wcmdscale(vegdist(t(abundance_matrix), method = "aitchison", pseudocount=1), k=2) %>% as.data.frame()

pcoa_df$Group <- sample_metadata$Group

pcoa_df$Group <- factor(pcoa_df$Group)

colnames(pcoa_df) <- c("Axis.1", "Axis.2", "Group")

pal_aaas("default")(length(levels(pcoa_df$Group)))

ggplot(data = pcoa_df, aes(x = Axis.1, y = Axis.2, color = Group)) +
  geom_point(size=5) +
  scale_color_manual(values = pal_aaas("default")(length(levels(pcoa_df$Group)))) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 25, face = "bold"),
    axis.text.y = element_text(size = 25, face = "bold"),
    legend.text = element_text(size = 25, face = "bold"),
    axis.title.x = element_text(size = 25, face = "bold"),
    axis.title.y = element_text(size = 25, face = "bold"),
    legend.title = element_text(size = 25, face = "bold")
  ) +
  geom_vline(xintercept = 0, color = "black", linetype = "dashed", linewidth = 1) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed", linewidth = 1)


############################ Heatmap ###########################################

heatmap_df <- as.data.frame(log10(rel_abundance_matrix+0.00001))

sample_metadata$Group <- factor(sample_metadata$Group)

col_list <- list(Group = setNames(pal_aaas("default")(length(levels(sample_metadata$Group))) ,levels(sample_metadata$Group)))

row_annotate <- rowAnnotation(
  df = sample_annotation,
  col = col_list, show_annotation_name = FALSE)

row_dend <-  hclust(dist(t(heatmap_df)), method = "complete")

column_dend <- hclust(dist(heatmap_df), method = "complete")

tmp <- heatmap_df %>% pivot_longer(cols = 1:length(colnames(heatmap_df)), names_to = "Sample",
                                   values_to = "Abundance")

col_fun <- colorRamp2(c(as.vector(min(tmp$Abundance)),
                        as.vector(quantile(tmp$Abundance)[3]),
                        as.vector(max(tmp$Abundance))), c("#1010fe", "white", "#FF1212"))

heatmap_plot <- Heatmap(t(heatmap_df), heatmap_legend_param = list(title = expression("Log"[10]*"Relative Abundance"),
                                                                   title_gp = gpar(fontsize = 10)),
                        row_names_gp = gpar(fontsize = 10, fontface = "bold"),
                        cluster_rows = color_branches(row_dend),
                        cluster_columns = color_branches(column_dend),
                        show_column_names = FALSE,
                        show_row_names = TRUE,
                        right_annotation = row_annotate,
                        col = col_fun(seq(min(tmp$Abundance), max(tmp$Abundance))))

heatmap_ggplot <- as_ggplot(grid.grabExpr(print(heatmap_plot)))

print(heatmap_ggplot)
