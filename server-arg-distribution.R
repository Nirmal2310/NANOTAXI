arg_distribution_analysis <- reactive({
    data <- analyze_data_reactive()$countsmetadata
    
    sample_metadata <- analyze_data_reactive()$sample_metadata
    
    data$Group <- factor(data$Group, levels = c("Control", "Case"))
    
    data$Family <- gsub("unclassified", "Unclassified", data$Family)
    
    data$Classification <- gsub("unclassified", "Unclassified", data$Classification)
    
    fun_richness_plot <- function(data, treatment) {
      data_2 <- data %>% filter(Group == treatment)
      
      data_2$Family <- gsub("unclassified", "Unclassified", data_2$Family)
      
      data_2$Classification <- gsub("unclassified", "Unclassified", data_2$Classification)
      
      data_2 <- data_2 %>%
                    group_by(Family, Classification) %>%
                    summarise(ARG_Richness = n_distinct(ARO_term))

      data_2$Family <- as.factor(data_2$Family)

      data_2 <- data_2 %>% arrange(Family, ARG_Richness)

      empty_bar <- 4
      
      to_add <- data.frame(matrix(NA, empty_bar * nlevels(data_2$Family),
                                ncol(data_2)))
      colnames(to_add) <- colnames(data_2)
      
      to_add$Family <- rep(levels(data_2$Family), each = empty_bar)
      
      data_2 <- rbind(data_2, to_add)
      
      data_2 <- data_2 %>% arrange(as.factor(Family))
      
      data_2$id <- seq(1, nrow(data_2))

      label_data <- data_2
      
      number_of_bar <- nrow(label_data)
      
      angle <- 90 - 360 * (label_data$id - 0.5)/ number_of_bar
      
      label_data$hjust <- ifelse(angle < -90, 1, 0)
      
      label_data$angle <- ifelse(angle < -90, angle + 180, angle)

      base_data <- data_2 %>%
            group_by(Family) %>%
            summarize(start = min(id), end = max(id) - empty_bar) %>%
            rowwise() %>%
            mutate(title = mean(c(start, end)))

      grid_data <- base_data
      
      grid_data$end <- grid_data$end[ c( nrow(grid_data), 1:nrow(grid_data)-1)] + 1
      
      grid_data$start <- grid_data$start - 1
        
      grid_data <- grid_data[-1, ]

      y_data_richness <- floor(seq(5, floor(max(data_2$ARG_Richness %>%
            replace(., is.na(.), 0))), length.out = 5))

      richness_plot <- ggplot(data_2, aes(x = as.factor(id), y = ARG_Richness, fill = Family)) +
            geom_bar(aes(x = as.factor(id), y = ARG_Richness, fill = Family),
                        stat = "identity", alpha = 0.3) +
        geom_segment(data = grid_data, aes(x = end, y = y_data_richness[5],
                xend = start, yend = y_data_richness[5]),
                colour = "grey",  size = 0.3, inherit.aes = FALSE) +
        geom_segment(data = grid_data, aes(x = end, y = y_data_richness[4],
                    xend = start, yend = y_data_richness[4]),
                colour = "grey",  size = 0.3, inherit.aes = FALSE) +
        geom_segment(data=grid_data, aes(x = end, y = y_data_richness[3],
                    xend = start, yend = y_data_richness[3]),
                colour = "grey",  size = 0.3, inherit.aes = FALSE) +
        geom_segment(data=grid_data, aes(x = end, y = y_data_richness[2],
                    xend = start, yend = y_data_richness[2]),
                colour = "grey",  size = 0.3, inherit.aes = FALSE) +
        geom_segment(data=grid_data, aes(x = end, y = y_data_richness[1],
                    xend = start, yend = y_data_richness[1]),
                colour = "grey",  size = 0.3, inherit.aes = FALSE) +
        annotate("text", x = rep(max(data_2$id)-1,5), y = c(y_data_richness[1:5]),
                 label = c(y_data_richness[1:5]), color = "black", size = 3, angle = 0,
                 fontface = "bold", hjust = 1) +
        geom_bar(aes(x = as.factor(id), y = ARG_Richness, fill = Family),
            stat = "identity", alpha = 0.3) +
        ylim(-max(data_2$ARG_Richness + 10, na.rm = TRUE), max(data_2$ARG_Richness + 10, na.rm = TRUE)) +
        theme_minimal() +
            theme(
                    axis.text = element_blank(),
                    axis.title = element_blank(),
                    panel.grid = element_blank(),
                    legend.text = element_text(size = 15, face = "bold"),
                    legend.title = element_text(size = 15, face = "bold"),
                    plot.title = element_text(family = "sans", size = 18,
                                              face = "bold", hjust = 0.5,
                                              margin = margin(0,0,-10,0))) +
        coord_polar() +
        geom_text(data = label_data, aes(x = id, y = ARG_Richness + 5, label = Classification, 
                    hjust = hjust), color = "black", fontface = "bold",
                    alpha = 0.6, size = 3.5, angle = label_data$angle, inherit.aes = FALSE ) +
        geom_text(x = max(data_2$id+0.45), y = y_data_richness[5] + 10, label = "ARG Richness", 
                color = "black", size = 3.5, angle = 0, fontface = "bold", hjust = 1) +
        ggtitle(paste0("ARGs Richness Plot for ", treatment, " Samples"))
        
        return(richness_plot)
    }
    
    fun_abundance_plot <- function(data, treatment)
      {
        data_3 <- data %>% filter(Group == treatment)
        
        data_3 <- data_3 %>%
                    group_by(Classification, AMR_Gene_Family) %>%
                            summarise(Counts = sum(Normalized_counts)) %>%
                            mutate(Abundance = Counts / sum(Counts) * 100) %>% 
          select(Classification, AMR_Gene_Family, Abundance)
        
        data_3$Classification <- as.factor(data_3$Classification)
        
        data_3 <- data_3 %>% pivot_wider(names_from = Classification, values_from = Abundance)
        
        data_3 <- as.data.frame(data_3)
        
        rownames(data_3) <- data_3$AMR_Gene_Family
        
        data_3 <- data_3 %>% subset(., select = -AMR_Gene_Family) %>% replace(., is.na(.), 0)
        
        arg_row_dend <- hclust(dist(t(data_3)), method = "complete")
        
        arg_col_dend <- hclust(dist(data_3), method = "complete")
        
        abundance_plot <- Heatmap(t(data_3), name = "Abundance",
                row_names_gp = gpar(fontsize = 10, fontface = "bold"),
                cluster_rows = arg_row_dend,
                cluster_columns = arg_col_dend,
                show_column_names = TRUE,
                show_row_names = TRUE,
                column_names_rot = -45,
                column_names_gp = gpar(fontsize = 7, fontface = "bold"))
        
        return(as_ggplot(grid.grabExpr(print(abundance_plot))))
    }
    
    return(list(control_richness_plot = fun_richness_plot(data, "Control"),
                case_richness_plot = fun_richness_plot(data, "Case"), 
                control_abundance_plot = fun_abundance_plot(data, "Control"),
                case_abundance_plot = fun_abundance_plot(data, "Case"))
    )
})

observeEvent(input$upload_data, {
    
    plots_data <- arg_distribution_analysis()
    
    output$plot_control_circular_richness_plot <- renderPlot({
        plots_data$control_richness_plot}, height = 1000
        )
    output$plot_case_circular_richness_plot <- renderPlot({
      plots_data$case_richness_plot}, height = 1000
      )
    
    output$plot_control_abundance_heatmap <- renderPlot({
        plots_data$control_abundance_plot}, height = 500
        )
    output$plot_case_abundance_heatmap <- renderPlot({
      plots_data$case_abundance_plot}, height = 500
      )
    
    output$control_download_circular_richness_plot <- downloadHandler(
        filename = function() {
            paste("Control_Circular_Richness_plot_",
            Sys.Date(), ".png", sep = "")
        },
        content = function(file) {
            ggsave(file, plots_data$control_richness_plot, width = 12.27, height = 11.69,
                   units = "in", dpi = "retina", bg = "white")
        }
    )
    
    output$case_download_circular_richness_plot <- downloadHandler(
      filename = function() {
        paste("Case_Circular_Richness_plot_",
              Sys.Date(), ".png", sep = "")
      },
      content = function(file) {
        ggsave(file, plots_data$case_richness_plot, width = 12.27, height = 11.69,
               units = "in", dpi = "retina", bg = "white")
      }
    )
    
    output$download_control_abundance_heatmap <- downloadHandler(
        filename = function() {
            paste("Control_Abundance_HeatMap_",
            Sys.Date(), ".png", sep = "")
        },
        content = function(file) {
            ggsave(file, plots_data$control_abundance_plot, width = 20, height = 10,
                   units = "in", dpi = "retina", bg = "white")
        }
    )
    
    output$download_case_abundance_heatmap <- downloadHandler(
      filename = function() {
        paste("Case_Abundance_HeatMap_",
              Sys.Date(), ".png", sep = "")
      },
      content = function(file) {
        ggsave(file, plots_data$case_abundance_plot, width = 30, height = 15,
               units = "in", dpi = "retina", bg = "white")
      }
    )
})