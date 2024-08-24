cohort_analysis_plots <- reactive({
    
    sample_metadata <- analyze_data_reactive()$sample_metadata

    rel_abundance_data <- analyze_data_reactive()$rel_abundance_data

    abundance_matrix <- analyze_data_reactive()$abundance_matrix

    rel_abundance_matrix <- analyze_data_reactive()$rel_abundance_matrix

    lineage <- analyze_data_reactive()$lineage
    
    sample_metadata$Group <- factor(sample_metadata$Group)
    
    rownames(sample_metadata) <- sample_metadata$Sample_Id

    rel_abundance_filtered_data <- analyze_data_reactive()$rel_abundance_filtered_data

    rel_abundance_filtered_matrix <- analyze_data_reactive()$rel_abundance_filtered_matrix
    
    fun_stacked_bar_plot <- function(data, sample_metadata, lineage) {     
      
      stacked_df <- data

      required_col <- which(colnames(stacked_df)==lineage)

      stacked_df[,required_col] <- as.character(stacked_df[,required_col])

      stacked_df <- stacked_df %>% 
        pivot_longer(cols = -all_of(required_col), names_to = "Sample_Id", values_to = "Abundance") %>% 
        filter(Abundance > 0)

      stacked_df <- stacked_df %>% group_by(Sample_Id) %>%
        arrange(desc(Abundance)) %>% dplyr::slice(1:ifelse(n()<10,n(),10))

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

      stacked_barplot <- ggplot(stacked_df, aes(x = Sample_Id, y = Abundance, fill = stacked_df[[required_col]])) +
        geom_bar(stat = "identity", position = "stack", color = "black", linewidth=0.2) +
        theme_linedraw() +
        theme(
          axis.text.x = element_text(size = 10, face = "bold"),
          axis.text.y = element_text(size = 10, face = "bold"),
          axis.title.x = element_text(size = 10, face = "bold"),
          axis.title.y = element_text(size = 10, face = "bold"),
          legend.title = element_text(face = "bold"),
          legend.text = element_text(face = "bold"),
          legend.position = "right",
          title = element_text(size = 10, face = "bold"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()
        ) +
        guides(fill = guide_legend(ncol = 4)) +
        labs(x = "BARCODE", y="% RELATIVE ABUNDANCE", fill = lineage) + 
        scale_fill_manual(values = fill_colors) +
        scale_y_continuous(expand = c(0, 0), limits = c(0,102)) +
        ggtitle(label = paste0("Stacked Barplot Showing Top 10 ",lineage," Across All Samples"))
    }

    fun_alpha_diversity_plot <- function(data, sample_metadata, lineage) {
      
      alpha_diversity_data <- data
        
      required_col <- which(colnames(alpha_diversity_data)==lineage)
        
      alpha_diversity_data[,required_col] <- as.character(alpha_diversity_data[,required_col])
        
      alpha_diversity_data <- alpha_diversity_data %>%
        pivot_longer(cols = -(required_col), names_to = "Sample_Id", values_to = "Abundance") %>%
        filter(Abundance > 0)
      
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
        scale_color_manual(values = pal_aaas("default")(length(levels(sample_metadata$Group)))) +
        theme_linedraw() +
        labs(y= "Alpha Diversity", x = "") +
        stat_pvalue_manual(subset(alpha_div_p, Diversity=="Shannon"), label = "p.signif", y.position = max(
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
          legend.title=element_text(colour="black",size=14, face = "bold"),
          legend.text=element_text(colour="black", size=14, face = "bold"),
          plot.title = element_text(hjust = 0.5,size = rel(2)),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()
        )
      
      simpson_plot <- alpha_diversity_data %>% filter(Diversity == "Simpson") %>%
        ggplot(aes(x=Group, y=Value, color=Group)) +
        geom_boxplot() +
        scale_color_manual(values = pal_aaas("default")(length(levels(sample_metadata$Group)))) + 
        theme_linedraw() +
        labs(y= "Alpha diversity", x = "") +
        stat_pvalue_manual(subset(alpha_div_p, Diversity=="Simpson"), y.position = max(
          subset(alpha_diversity_data, Diversity=="Simpson")$Value) + 0.01, label = "p.signif",
          hide.ns = "p.adj", step.increase = 0.1, tip.length = 0.02, bracket.size = 0.8, size = 8) +
        guides(color = guide_legend(title = "Simpson", title.position = "top")) +
        theme(
          axis.title.x = element_text(size = 14, face = "bold"),
          axis.title.y = element_text(size = 14, face = "bold"),
          strip.text.x = element_text(size = 14, face = "bold"),
          axis.text.y=element_text(size=14, face = "bold"),
          axis.text.x= element_blank(),
          axis.ticks.x = element_blank(),
          legend.title=element_text(colour="black",size=14, face = "bold"),
          legend.text=element_text(colour="black", size=14, face = "bold"),
          plot.title = element_text(hjust = 0.5,size = rel(2)),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()
        )
      
      
      diversity_facet <- as_ggplot(grid.grabExpr(grid.arrange(shannon_plot,
                                          simpson_plot, ncol=2)))
      
      return(diversity_facet)
    }


    fun_pcoa_plot <- function(matrix) {

      pcoa_dist <- wcmdscale(vegdist(t(matrix), method = "aitchison", pseudocount=0.5), k=2, eig = TRUE)

      pcoa_df <- pcoa_dist$points[,1:2] %>% as.data.frame()
      
      pcoa_eigenvalues <- pcoa_dist$eig

      pcoa.var <- round(pcoa_eigenvalues/sum(pcoa_eigenvalues)*100, 1)  

      pcoa_df$Group <- sample_metadata$Group

      pcoa_df$Group <- factor(pcoa_df$Group)
      
      colnames(pcoa_df) <- c("Axis.1", "Axis.2", "Group")

      pcoa_plot <- ggplot(data = pcoa_df, aes(x = Axis.1, y = Axis.2, color = Group)) +
        geom_point(size=5) +
        stat_ellipse(aes(colour = Group, fill = Group), level = 0.95, alpha = 0.25, geom = "polygon") +
        scale_color_manual(values = pal_aaas("default")(length(levels(sample_metadata$Group)))) +
        scale_fill_manual(values = pal_aaas("default")(length(levels(sample_metadata$Group)))) +
        theme_linedraw() +
        xlab(paste0("PC1 (", pcoa.var[1], "%", ")")) +
        ylab(paste0("PC2 (", pcoa.var[2], "%", ")")) +
        geom_vline(xintercept = 0, color = "black", linetype = "dashed") +
        geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
        theme(
          axis.text.x = element_text(size = 15, face = "bold"),
          axis.text.y = element_text(size = 15, face = "bold"),
          legend.text = element_text(size = 15, face = "bold"),
          axis.title.x = element_text(size = 15, face = "bold"),
          axis.title.y = element_text(size = 15, face = "bold"),
          legend.title = element_text(size = 15, face = "bold"),
          title = element_text(size = 15, face = "bold")
        ) +
        ggtitle("Principal Coordination Analysis (PCoA) ordination plot using Aitchison Distance")

      return(pcoa_plot)
    }

    fun_pca_plot <- function(matrix) {
      
      clr_transformed_abundance_matrix <- clr(matrix+0.5) %>% as.data.frame()

      pca <- prcomp(t(clr_transformed_abundance_matrix), scale. = FALSE, center = TRUE,
              retx = TRUE)

      pca.var <- pca$sdev^2

      pca.var.per <- round(pca.var/sum(pca.var)*100, 1)

      pca_df <- data.frame(Sample_Id = rownames(pca$x), X=pca$x[,1], Y=pca$x[,2])

      pca_df <- inner_join(pca_df, sample_metadata)

      pca_df$Group <- factor(pca_df$Group)

      pca_plot <- ggplot(data = pca_df, aes(x = X, y = Y, color = Group)) +
        geom_point(size=5) +
        stat_ellipse(aes(colour = Group, fill = Group), level = 0.95, alpha = 0.25, geom = "polygon") +
        xlab(paste0("PC1 (", pca.var.per[1], "%", ")")) +
        ylab(paste0("PC2 (", pca.var.per[2], "%", ")")) +
        scale_color_manual(values = pal_aaas("default")(length(levels(sample_metadata$Group)))) +
        scale_fill_manual(values = pal_aaas("default")(length(levels(sample_metadata$Group)))) +
        theme_linedraw() +
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
        ggtitle("Principal Component Analysis (PCA) Plot using CLR Transformed Data")

      return(pca_plot)
    }

    fun_nmds_plot <- function(matrix) {
      
      nmds_dist <- metaMDS(t(matrix), distance = "bray", trymax = 100)

      nmds_df <- as.data.frame(nmds_dist$points)

      nmds_df$Group <- sample_metadata$Group

      nmds_plot <- ggplot(data = nmds_df, aes(x = MDS1, y = MDS2, color = Group)) +
        geom_point(size=5) +
        stat_ellipse(aes(colour = Group, fill = Group), level = 0.95, alpha = 0.25, geom = "polygon") +
        scale_color_manual(values = pal_aaas("default")(length(levels(sample_metadata$Group)))) +
        scale_fill_manual(values = pal_aaas("default")(length(levels(sample_metadata$Group)))) +
        theme_linedraw() +
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
      
      return(nmds_plot)
    }

    fun_permanova <- function(matrix)
    {
      perm_dist <- vegdist(t(matrix), method = "aitchison", pseudocount=0.5)

      permanova_res <- pairwise.adonis(perm_dist, as.factor(sample_metadata$Group), p.adjust.m = "BH")

      permanova_res <- permanova_res %>% dplyr::select(c(pairs, R2, p.value, p.adjusted))

      colnames(permanova_res) <- c("Pair", "R2", "P", "Padj")

      return(permanova_res)
    
    }

    fun_heatmap_plot <- function(matrix)
    {
      
        heatmap_df <- as.data.frame(log10(matrix+0.00001))

        sample_annotation <- data.frame(sample_metadata$Group)

        sample_metadata$Group <- factor(sample_metadata$Group)

        colnames(sample_annotation) <- "Group"

        sample_annotation$Group <- factor(sample_annotation$Group)

        col_list <- list(Group = setNames(pal_aaas("default")(length(levels(sample_metadata$Group))), 
        levels(sample_metadata$Group)))

        row_annotate <- rowAnnotation(
        df = sample_annotation,
        col = col_list, show_annotation_name = FALSE)

        row_dend <-  hclust(dist(t(heatmap_df)), method = "complete")

        column_dend <- hclust(dist(heatmap_df), method = "complete")

        tmp <- heatmap_df %>% pivot_longer(cols = 1:length(colnames(heatmap_df)), names_to = "Sample",
                                          values_to = "Abundance")

        col_fun <- colorRamp2(c(as.vector(quantile(tmp$Abundance)[1]),
                              as.vector(quantile(tmp$Abundance)[3]),
                              as.vector(quantile(tmp$Abundance)[5])), c("#1010fe", "#FFD700", "#FF1212"))

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

      return(heatmap_ggplot)
    }
    
    return(list(stacked_barplot = fun_stacked_bar_plot(rel_abundance_data, sample_metadata, lineage),
                alpha_diversity_plot = fun_alpha_diversity_plot(rel_abundance_filtered_data, sample_metadata, lineage),
                pcoa_plot = fun_pcoa_plot(abundance_matrix),
                pca_plot = fun_pca_plot(abundance_matrix),
                nmds_plot = fun_nmds_plot(rel_abundance_matrix),
                permanova_res = fun_permanova(abundance_matrix),
                heatmap_plot = fun_heatmap_plot(rel_abundance_filtered_matrix)
                ))
})

observeEvent(c(input$upload_data, input$taxa), ignoreNULL=TRUE, ignoreInit=TRUE, {
    
    
    select_taxon <- reactive({
      
      return(list('taxon'=input$taxa))
      
    })
  
    lineage <- select_taxon()$taxon
  
    plots_data <- cohort_analysis_plots()

    output$plot_stacked_barplot <- renderPlot({
        plots_data$stacked_barplot}, height = 700, width = 1500
        )
    
    output$plot_boxplot <- renderPlot({
        plots_data$alpha_diversity_plot}, height = 500
        )
    output$plot_pcoa <- renderPlot({
      plots_data$pcoa_plot}, height = 500
    )


    output$plot_nmds <- renderPlot({
      plots_data$nmds_plot}, height = 500
    )

    output$plot_pca <- renderPlot({
        plots_data$pca_plot}, height = 500
        )

    output$permanova_data <- renderDataTable({
  
      print("PERMANOVA Results")
      
      plots_data$permanova_res
      
    })
    
    output$plot_heatmap <- renderPlot({
      plots_data$heatmap_plot}, height = 500
      )

    output$download_stacked_barplot <- downloadHandler(
      filename = function() {
        paste("Stacked_bar_plot_", lineage,"_", Sys.Date(), ".png", sep = "")
      },
      content = function(file) {
        ggsave(file, plots_data$stacked_barplot,
               width = 15, height = 10, units = "in", dpi = "retina", bg = "white")
      }
    )
    
    output$download_boxplot <- downloadHandler(
      filename = function() {
        paste("Alpha_Diversity_plot_", lineage, "_", Sys.Date(), ".png", sep = "")
      },
      content = function(file) {
        ggsave(file, plots_data$alpha_diversity_plot,
               width = 13.69, height = 8.27, units = "in", dpi = "retina", bg = "white")
      }
    )
    
    output$download_pcoa_plot <- downloadHandler(
      filename = function() {
        paste("PCoA_plot_", lineage, "_", Sys.Date(), ".png", sep = "")
      },
      content = function(file) {
        ggsave(file, plots_data$pcoa_plot,
               width = 13.69, height = 8.27, units = "in", dpi = "retina", bg = "white")
      }
    )
    
    output$download_nmds_plot <- downloadHandler(
      filename = function() {
        paste("NMDS_plot_", lineage, "_", Sys.Date(), ".png", sep = "")
      },
      content = function(file) {
        ggsave(file, plots_data$nmds_plot,
               width = 13.69, height = 8.27, units = "in", dpi = "retina", bg = "white")
      }
    )
    
    output$download_pca_plot <- downloadHandler(
      filename = function() {
        paste("PCA_plot_", lineage, "_", Sys.Date(), ".png", sep = "")
      },
      content = function(file) {
        ggsave(file, plots_data$pca_plot,
               width = 13.69, height = 8.27, units = "in", dpi = "retina", bg = "white")
      }
    )
    
    output$download_heatmap <- downloadHandler(
      filename = function() {
        paste("HeatMap_", lineage, "_", Sys.Date(), ".png", sep = "")
      },
      content = function(file) {
        ggsave(file, plots_data$heatmap_plot,
               width = 13.69, height = 8.27, units = "in", dpi = "retina", bg = "white")
      }
    )

    output$download_permanova_csv <- downloadHandler(
      filename = paste0("PERMANOVA_Result_", Sys.Date(), ".csv"),
      
      content = function(file) {
        
        write.csv(plots_data$permanova_res,
                  
                  file, row.names = FALSE)
      }
    )
})