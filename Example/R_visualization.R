library(tidyverse)
library(stringr)
library(ggpubr)
library(dendextend)
library(ComplexHeatmap)
library(vegan)
library(ggrepel)
library(grid)
library(gridExtra)
library(phyloseq)
library(ggforce)
library(ggsci)
library(scales)
library(ComplexHeatmap)
library(dendextend)
library(viridis)
library(circlize)
library(compositions)
library(pairwiseAdonis)
library(ANCOMBC)
library(phyloseq)
library(lme4)

file_list <- gsub(".txt", "", list.files(getwd())[grep("\\_final_blast_result.txt$", 
                                                       list.files(getwd()))])
samples_header <- gsub("_final_blast_result","",file_list)

sample_metadata <- read.table(file = "Sample_Information.csv", sep = ",", 
                              header = TRUE)

rownames(sample_metadata) <- sample_metadata$Sample_Id

sample_annotation <- data.frame(sample_metadata$Group)

colnames(sample_annotation) <- "Group"

sample_annotation$Group <- factor(sample_annotation$Group)

rownames(sample_annotation) <- sample_metadata$Sample_Id

sample_data_list <- list()

rel_abundance_data_list <- list()

abundance_data_list <- list()

lineage <- "Species"

phyloseq_df <- NULL

for (i in 1:length(file_list))
{
  sample_data_list[[i]] <- read.delim(file = paste0(file_list[i],".txt"), header = TRUE)
  
  colnames(sample_data_list[[i]]) <- c("TAX_ID", "Counts", "Superkingdom", "Phylum",
                                       "Class", "Order", "Family", "Genus", "Species")
  
  sample_data_list[[i]][sample_data_list[[i]]==""] <- "Unclassified"
  
  sample_data_list[[i]]$Sample_Id <- samples_header[i]
  
  phyloseq_df <- rbind(phyloseq_df, sample_data_list[[i]])
  
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

rel_abundance_matrix <- as.matrix(rel_abundance_data[,-1])

rownames(rel_abundance_matrix) <- rel_abundance_data[,which(colnames(rel_abundance_data)==lineage)]

rel_abundance_filtered_matrix <- rel_abundance_matrix[apply(rel_abundance_matrix, 1, function(row) any(mean(row) >0.1 )), ]

abundance_data <- abundance_data_list %>% purrr::reduce(full_join, by=lineage)

abundance_data <- as.data.frame(abundance_data)

abundance_data[is.na(abundance_data)] <- 0

abundance_matrix <- abundance_data[,-1]

rownames(abundance_matrix) <- abundance_data[,which(colnames(abundance_data)==lineage)]

rel_abundance_filtered_data <- rel_abundance_data[which(rel_abundance_data[,1] %in% rownames(rel_abundance_filtered_matrix)),]

# abundance_matrix <- abundance_matrix[which(rownames(abundance_matrix) %in% rownames(rel_abundance_matrix)),]

############# Stacked Bar-Plot #################################################

stacked_df <- rel_abundance_data

required_col <- which(colnames(stacked_df)==lineage)

stacked_df[,required_col] <- as.character(stacked_df[,required_col])

stacked_df <- stacked_df %>% 
  pivot_longer(cols = -(all_of(required_col)), names_to = "Sample_Id", values_to = "Abundance") %>% 
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

#################### Alpha Diversity Plot ######################################

alpha_diversity_data <- rel_abundance_filtered_data

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
  stat_pvalue_manual(subset(alpha_div_p, Diversity=="Shannon"), label = "p.adj", y.position = max(
    subset(alpha_diversity_data, Diversity=="Shannon")$Value) + 0.1, hide.ns = "p.adj", step.increase = 0.1,
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
    subset(alpha_diversity_data, Diversity=="Simpson")$Value) + 0.01, label = "p.adj",
    hide.ns = "p.adj", step.increase = 0.1, tip.length = 0.02, bracket.size = 0.8, size = 8) +
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

#################### PCA Plot ##################################################

clr_transformed_abundance_matrix <- clr(abundance_matrix+0.5) %>% as.data.frame()

pca <- prcomp(t(clr_transformed_abundance_matrix), scale. = FALSE, center = TRUE,
              retx = TRUE)

pca.var <- pca$sdev^2

pca.var.per <- round(pca.var/sum(pca.var)*100, 1)

pca_data <- data.frame(Sample_Id = rownames(pca$x), X=pca$x[,1], Y=pca$x[,2])

pca_data <- inner_join(pca_data, sample_metadata)

pca_data$Group <- factor(pca_data$Group)

ggplot(data = pca_data, aes(x = X, y = Y, color = Group)) +
  geom_point(size=5) +
  stat_ellipse(aes(colour = Group, fill = Group), level = 0.95, alpha = 0.25, geom = "polygon") +
  xlab(paste0("PC1 (", pca.var.per[1], "%", ")")) +
  ylab(paste0("PC2 (", pca.var.per[2], "%", ")")) +
  scale_color_manual(values = pal_aaas("default")(length(levels(sample_annotation$Group)))) +
  scale_fill_manual(values = pal_aaas("default")(length(levels(pcoa_df$Group)))) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 15, face = "bold"),
    axis.text.y = element_text(size = 15, face = "bold"),
    legend.text = element_text(size = 15, face = "bold"),
    axis.title.x = element_text(size = 15, face = "bold"),
    axis.title.y = element_text(size = 15, face = "bold"),
    legend.title = element_text(size = 15, face = "bold"),
    title = element_text(size = 15, face = "bold")
  ) +
  geom_vline(xintercept = 0, color = "black") +
  geom_hline(yintercept = 0, color = "black") +
  ggtitle("Principal Component Analysis (PCA) Plot on CLR Transformed Data")

######################### NMDS Plot #############################################

nmds_dist <- metaMDS(t(rel_abundance_matrix), distance = "bray", trymax = 100)

nmds_df <- as.data.frame(nmds_dist$points)

nmds_df$Group <- sample_metadata$Group

ggplot(data = nmds_df, aes(x = MDS1, y = MDS2, color = Group)) +
  geom_point(size=5) +
  stat_ellipse(aes(colour = Group, fill = Group), level = 0.95, alpha = 0.25, geom = "polygon") +
  scale_color_manual(values = pal_aaas("default")(length(levels(sample_annotation$Group)))) +
  scale_fill_manual(values = pal_aaas("default")(length(levels(pcoa_df$Group)))) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 15, face = "bold"),
    axis.text.y = element_text(size = 15, face = "bold"),
    legend.text = element_text(size = 15, face = "bold"),
    axis.title.x = element_text(size = 15, face = "bold"),
    axis.title.y = element_text(size = 15, face = "bold"),
    legend.title = element_text(size = 15, face = "bold"),
    title = element_text(size = 15, face = "bold")
  ) +
  geom_vline(xintercept = 0, color = "black", linetype = "dashed") +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  ggtitle("Non-metric MultiDimensional Scaling (NMDS) ordination plot with Bray-Curtis Distance")

#################### PCoA Plot #####################################################

pcoa_dist <- wcmdscale(vegdist(t(abundance_matrix), method = "aitchison", pseudocount=1), k=2, eig = TRUE)

pcoa_df <- pcoa_dist$points[,1:2] %>% as.data.frame()

pcoa_eigenvalues <- pcoa_dist$eig

pcoa.var <- round(pcoa_eigenvalues/sum(pcoa_eigenvalues)*100, 1)

pcoa_df$Group <- sample_metadata$Group

pcoa_df$Group <- factor(pcoa_df$Group)

colnames(pcoa_df) <- c("Axis.1", "Axis.2", "Group")

ggplot(data = pcoa_df, aes(x = Axis.1, y = Axis.2, color = Group)) +
  geom_point(size=5) +
  stat_ellipse(aes(colour = Group, fill = Group), level = 0.95, alpha = 0.25, geom = "polygon") +
  scale_color_manual(values = pal_aaas("default")(length(levels(pcoa_df$Group)))) +
  scale_fill_manual(values = pal_aaas("default")(length(levels(pcoa_df$Group)))) +
  theme_minimal() +
  xlab(paste0("PC1 (", pcoa.var[1], "%", ")")) +
  ylab(paste0("PC2 (", pcoa.var[2], "%", ")")) +
  theme(
    axis.text.x = element_text(size = 15, face = "bold"),
    axis.text.y = element_text(size = 15, face = "bold"),
    legend.text = element_text(size = 15, face = "bold"),
    axis.title.x = element_text(size = 15, face = "bold"),
    axis.title.y = element_text(size = 15, face = "bold"),
    legend.title = element_text(size = 15, face = "bold"),
    title = element_text(size = 15, face = "bold")
  ) +
  geom_vline(xintercept = 0, color = "black") +
  geom_hline(yintercept = 0, color = "black") +
  ggtitle("Principal Coordination Analysis (PCoA) ordination plot with Aitchison Distance")


##################### Permanova Analysis #######################################

perm_dist <- vegdist(t(abundance_matrix), method = "aitchison", pseudocount=1)

permanova_res <- pairwise.adonis(perm_dist, as.factor(sample_metadata$Group))

permanova_res <- permanova_res %>% select(c(pairs, R2, p.value, p.adjusted, sig))

############### HeatMap Plot ###################################################

# heatmap_data <- as.data.frame(rel_abundance_matrix)

# sample_annotation <- data.frame(sample_metadata$Group)
# 
# colnames(sample_annotation) <- "Group"
# 
# sample_annotation$Group <- factor(sample_annotation$Group)
# 
# col_list <- list(Group = setNames(pal_aaas("default")(length(levels(pcoa_df$Group))) ,levels(pcoa_df$Group)))
# 
# row_annotate <- rowAnnotation(
#   df = sample_annotation,
#   col = col_list, show_annotation_name = FALSE)
# 
# row_dend <-  hclust(dist(t(heatmap_data)), method = "complete")
# 
# column_dend <- hclust(dist(heatmap_data), method = "complete")
# 
# heatmap_plot <- Heatmap(t(heatmap_data), name = "Relative Abundance",
#                         row_names_gp = gpar(fontsize = 10, fontface = "bold"),
#                         cluster_rows = color_branches(row_dend),
#                         cluster_columns = color_branches(column_dend),
#                         show_column_names = FALSE,
#                         show_row_names = TRUE,
#                         right_annotation = row_annotate,
#                         col = col_fun(seq(0,30)))
# 
# print(heatmap_plot)

heatmap_df <- as.data.frame(log10(rel_abundance_filtered_matrix+0.00001))

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


# pheatmap::pheatmap(t(heatmap_df),
#          col = col_fun(seq(min(tmp$Abundance),max(tmp$Abundance))),
#          cluster_rows = T, cluster_cols = T,
#          clustering_distance_cols = 'euclidean',
#          clustering_distance_rows = 'euclidean',
#          clustering_method = 'complete',
#          show_rownames = T,
#          show_colnames = F,
#          fontsize_col = 10,
#          annotation_row = sample_annotation,
#          annotation_colors = col_list)

heatmap_plot <- Heatmap(t(heatmap_df), heatmap_legend_param = list(title = expression("Log"[10]*"Relative Abundance"),
                                                                   title_gp = gpar(fontsize = 10)),
                        row_names_gp = gpar(fontsize = 10, fontface = "bold"),
                        cluster_rows = color_branches(row_dend),
                        cluster_columns = color_branches(column_dend),
                        show_column_names = FALSE,
                        show_row_names = TRUE,
                        right_annotation = row_annotate,
                        col = col_fun(seq(min(tmp$Abundance), max(tmp$Abundance)))
                        )

heatmap_ggplot <- as_ggplot(grid.grabExpr(print(heatmap_plot)))

print(heatmap_ggplot)

################## Miscellaneous ###############################################

perm_dist <- vegdist(t(abundance_matrix), method = "bray")


permanova_test <- adonis2(perm_dist~as.factor(sample_metadata$Group), data = perm_dist, permutations = 9999)


tax_df <- phyloseq_df %>% dplyr::select(-c(TAX_ID, Counts, Sample_Id)) %>% unique()

library(phyloseq)

tax_df <- tax_df %>% as.matrix()

rownames(tax_df) <- tax_df[,7]

abundance_phyloseq <- phyloseq(otu_table(abundance_matrix, taxa_are_rows = TRUE), tax_table(tax_df), sample_data(sample_metadata))

abundance_phyloseq@sam_data$Group <- factor(abundance_phyloseq@sam_data$Group)

ancombc2(data = abundance_phyloseq, assay_name = "otu_table", tax_level = "Species",
         fix_formula = NULL,
         rand_formula = NULL,
         p_adj_method = "holm", pseudo_sens = TRUE,
         prv_cut = 0.10, lib_cut = 1000, s0_perc = 0.05,
         group = "Group", struc_zero = TRUE, neg_lb = TRUE,
         alpha = 0.05, n_cl = 1, verbose = TRUE,
         global = TRUE, pairwise = TRUE, dunnet = TRUE, trend = TRUE,
         iter_control = list(tol = 1e-2, max_iter = 1, verbose = TRUE),
         em_control = list(tol = 1e-5, max_iter = 1),
         lme_control = lme4::lmerControl(),
         mdfdr_control = list(fwer_ctrl_method = "holm", B = 1),
         trend_control = list(contrast =
                                list(matrix(c(1, 0, -1, 1),
                                            nrow = 2,
                                            byrow = TRUE)),
                              node = list(2),
                              solver = "ECOS",
                              B = 1))

abundance_phyloseq@sam_data$Group <- factor(abundance_phyloseq@sam_data$Group)

ancom_res <- ancombc2(data=abundance_phyloseq, rank = "Species", fix_formula = "Group",
         rand_formula = NULL, pseudo_sens = TRUE, prv_cut = 0.10, lib_cut = 1000,
         s0_perc = 0.05, group = "Group", struc_zero = TRUE, neg_lb = TRUE,
         alpha = 0.05, n_cl = 1, verbose = TRUE, global = TRUE, pairwise = TRUE,
         dunnet = TRUE, trend = FALSE, iter_control = list(tol = 1e-2, max_iter = 100, verbose = TRUE),
         em_control = list(tol = 1e-5, max_iter = 100, verbose = TRUE),
         lme_control = lme4::lmerControl(),
         mdfdr_control = list(fwer_ctrl_method = "holm", B = 100),
         trend_control = NULL)


ancom_res$res_pair %>% View()


res_pair <- ancom_res$res_pair

df_pair1 <- res_pair %>% dplyr::filter(`diff_GroupJoint aspirate` == 1 |
                             `diff_GroupPeritoneal dialysis fluid` == 1 | 
                             `diff_GroupPeritoneal fluid` == 1) %>% 
  dplyr::mutate(lfc1 = ifelse(`diff_GroupJoint aspirate`==1, round(`lfc_GroupJoint aspirate`,2),0),
                lfc2 = ifelse(`diff_GroupPeritoneal dialysis fluid`==1, round(`lfc_GroupPeritoneal dialysis fluid`,2),0),
                lfc3 = ifelse(`diff_GroupPeritoneal dialysis fluid_GroupJoint aspirate`==1, round(`lfc_GroupPeritoneal dialysis fluid_GroupJoint aspirate`,2),0)) %>% 
  tidyr::pivot_longer(cols = lfc1:lfc3,
                      names_to = "Group", values_to = "value") %>% 
  dplyr::arrange(taxon)

df_pair2 <- res_pair %>% dplyr::filter(`diff_GroupPeritoneal fluid` == 1 |
                                         `diff_GroupPeritoneal dialysis fluid` == 1 | 
                                         `diff_GroupPeritoneal fluid` == 1) %>% 
  dplyr::mutate(lfc1 = ifelse(`passed_ss_GroupJoint aspirate`==1 & `diff_GroupJoint aspirate`==1, "aquamarine3", "black"),
                lfc2 = ifelse(`passed_ss_GroupPeritoneal dialysis fluid`==1 & `diff_GroupPeritoneal dialysis fluid`==1, "aquamarine3", "black"),
                lfc3 = ifelse(`passed_ss_GroupPeritoneal fluid`==1 & `diff_GroupPeritoneal fluid`==1, "aquamarine3", "black")) %>% 
  tidyr::pivot_longer(cols = lfc1:lfc3,
                      names_to = "Group", values_to = "color") %>% 
  dplyr::arrange(taxon)


df_fig_pair <- df_pair1 %>%
  dplyr::left_join(df_pair2, by = c("taxon", "Group"))

df_fig_pair$Group = recode(df_fig_pair$Group, 
                           `lfc1` = "Overweight - Obese",
                           `lfc2` = "Lean - Obese",
                           `lfc3` = "Lean - Overweight")

df_fig_pair$Group = factor(df_fig_pair$Group, 
                           levels = c("Overweight - Obese",
                                      "Lean - Obese", 
                                      "Lean - Overweight"))

lo = floor(min(df_fig_pair$value))
up = ceiling(max(df_fig_pair$value))
mid = (lo + up)/2
df_fig_pair %>%
  ggplot(aes(x = Group, y = taxon, fill = value)) + 
  geom_tile(color = "black") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       na.value = "white", midpoint = mid, limit = c(lo, up),
                       name = NULL) +
  geom_text(aes(Group, taxon, label = value, color = color), size = 4) +
  scale_color_identity(guide = FALSE) +
  labs(x = NULL, y = NULL, title = "Log fold changes as compared to obese subjects") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))


df_fig_pair %>% dplyr::select(c(Group, taxon, value, color))

tmp <- read.table(file="qual.txt", header = FALSE)

colnames(tmp) <- "Phred"

table(tmp$Phred)

tmp <- tmp %>%
  mutate(Bin = cut(Phred, breaks = c(seq(1,max(Phred),1),max(Phred)), include.lowest = TRUE, right = TRUE, labels = as.character(c(seq(1,max(Phred),1))))) %>%
  select(-Phred) %>%
  group_by(Bin) %>% summarise(Counts = n())

tmp %>% ggplot(aes(Bin,Counts)) +
  geom_bar(stat = "identity", fill = "#00005E") +
  theme_classic() +
  xlab("Phred Score") +
  ylab("Frequency") +
  theme(
    axis.title.x = element_text(size = 14, face = "bold", colour = "black"),
    axis.title.y = element_text(size = 14, face = "bold", colour = "black"),
    strip.text.x = element_text(size = 14, face = "bold", colour = "black"),
    axis.text.y= element_text(size=14, face = "bold", colour = "black"),
    axis.text.x= element_text(size=14, face = "bold", colour = "black"),
    legend.title= element_text(colour="black",size=14, face = "bold"),
    legend.text= element_text(colour="black", size=14, face = "bold"),
    axis.line = element_line(colour = "black", linewidth = 0.5, linetype = "solid" ),
    strip.background = element_blank(),
    legend.position = "right"
  )
  
