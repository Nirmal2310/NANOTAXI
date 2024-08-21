beta_diversity_analysis <- reactive({
    
    data <- analyze_data_reactive()$countsmetadata
    
    sample_metadata <- analyze_data_reactive()$sample_metadata
    
    fun_heatmap_plot <- function(data, sample_metadata) {
        data_4 <- data %>%
          group_by(Sample_Id, ARO_term) %>%
          summarise(Counts = sum(Normalized_counts)) %>%
          mutate(Counts = log2(Counts))
        
        data_4 <- data_4 %>%
          pivot_wider(names_from = Sample_Id, values_from = Counts)
        
        data_4 <- as.data.frame(data_4)
        
        rownames(data_4) <- data_4$ARO_term
        
        data_4 <- data_4 %>%
          subset(. , select = -ARO_term) %>%
          replace(. , is.na(.), 0)
        
        sample_annotation <- data.frame(sample_metadata$Group)
        
        colnames(sample_annotation) <- "Group"
        
        sample_annotation$Group <- factor(sample_annotation$Group, levels = c("Control", "Case"))
        
        col_list <- list(Group = 
                           setNames(rep(c("#746AB0", "#ffbf00"),
                                        c(length(sample_annotation$Group[sample_annotation$Group=="Case"]),
                                          length(sample_annotation$Group[sample_annotation$Group=="Control"]))),sample_annotation$Group))
        
        row_annotate <- rowAnnotation(
          df = sample_annotation,
          col = col_list, show_annotation_name = FALSE)
        
        row_dend <-  hclust(dist(t(data_4)), method = "complete")
        
        column_dend <- hclust(dist(data_4), method = "complete")
        
        heatmap_plot <- Heatmap(t(data_4), name = "Log2Counts",
                                row_names_gp = gpar(fontsize = 10, fontface = "bold"),
                                cluster_rows = color_branches(row_dend),
                                cluster_columns = color_branches(column_dend),
                                show_column_names = FALSE,
                                show_row_names = TRUE,
                                right_annotation = row_annotate)
        
        heatmap_ggplot <- as_ggplot(grid.grabExpr(print(heatmap_plot)))
        
        return(heatmap_ggplot)
    }
    fun_pca_plot <- function(data) {
        data_4 <- data %>%
          group_by(Sample_Id, ARO_term) %>%
          summarise(Counts = sum(Normalized_counts)) %>%
          mutate(Counts = log2(Counts))
        data_4 <- data_4 %>%
          pivot_wider(names_from = Sample_Id, values_from = Counts)
        data_4 <- as.data.frame(data_4)
        rownames(data_4) <- data_4$ARO_term
        data_4 <- data_4 %>%
          subset(. , select = -ARO_term) %>%
          replace(. , is.na(.), 0)
        pca <- prcomp(t(data_4), scale. = TRUE, center = TRUE)
        pca.var <- pca$sdev^2
        pca.var.per <- round(pca.var / sum(pca.var) * 100, 1)
        pca.data <- data.frame(Sample_Id = rownames(pca$x), X=pca$x[,1],Y=pca$x[,2])
        pca.data <- inner_join(pca.data, sample_metadata, by = "Sample_Id")
        pca.data$Group <- factor(pca.data$Group, levels = c("Control", "Case"))
        pca_plot <- ggplot(data = pca.data, aes(x = X, y = Y, color = Group)) +
          geom_point() +
          scale_color_manual(values = c( "#ffbf00", "#746AB0")) +
          geom_mark_ellipse(aes(color = Group), expand = unit(0.5, "mm"), 
                            linetype = 2) +
          xlab(paste0("PC1 (", pca.var.per[1], "%", ")")) +
          ylab(paste0("PC2 (", pca.var.per[2], "%", ")")) +
          theme_bw() +
          theme(
            axis.text.x = element_text(size = 25, face = "bold"),
            axis.text.y = element_text(size = 25, face = "bold"),
            legend.text = element_text(size = 25, face = "bold"),
            axis.title.x = element_text(size = 25, face = "bold"),
            axis.title.y = element_text(size = 25, face = "bold"),
            legend.title = element_text(size = 25, face = "bold"),
            panel.border = element_rect(color = "black", size = 1, linetype = "solid")
          )
        return(pca_plot)
    }
    return(list(heatmap_plot = fun_heatmap_plot(data, sample_metadata),
    pca_plot = fun_pca_plot(data)))
})
observeEvent(input$upload_data, {
    beta_diversity_plots <- beta_diversity_analysis()

    output$plot_heatmap <- renderPlot({
      beta_diversity_plots$heatmap_plot}, height = 500)

    output$plot_pca <- renderPlot({
      beta_diversity_plots$pca_plot}, height = 500)

    output$download_heatmap <- downloadHandler(
        filename = function() {
            paste("HeatMap_", Sys.Date(),".png", sep = "")
        },
        content = function(file) {
            ggsave(file, beta_diversity_plots$heatmap_plot,
                   width = 14.69, height = 8.27, units = "in", dpi = "retina")
        }
    )

    output$download_pca <- downloadHandler(
        filename = function() {
            paste("PCA_plot_", Sys.Date(), ".png", sep = "")
        },
        content = function(file) {
            ggsave(file, beta_diversity_plots$pca_plot,
                   width = 11.69, height = 8.27, units = "in", dpi = "retina")
        }
    )
})