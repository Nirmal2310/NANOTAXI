alpha_diversity_analysis <- reactive({
    data <- analyze_data_reactive()$countsmetadata
    sample_metadata <- analyze_data_reactive()$sample_metadata
    fun_diversity_plot <- function(data)
    {
        diversity_data <- data %>% group_by(Sample_Id,ARO_term) %>%
          summarise(Counts = sum(Normalized_counts)) %>%
          mutate(Abundance = Counts/sum(Counts)*100)
        diversity_data <- inner_join(diversity_data, sample_metadata, by="Sample_Id")
        diversity_data <- diversity_data %>% select(Sample_Id, Group, everything())
        diversity_data <- diversity_data %>% group_by(Sample_Id) %>% summarise(Shannon = diversity(Abundance, index = "shannon"),
        Simpson = diversity(Abundance, index = "simpson"))
        diversity_data <- inner_join(diversity_data, sample_metadata, by="Sample_Id")
        diversity_data <- diversity_data %>% select(Sample_Id, Group, everything())
        diversity_data$Group <- factor(diversity_data$Group,
                                levels = c("Control", "Case"))
        diversity_data <- pivot_longer(diversity_data, cols = c(Shannon, Simpson),
                                        names_to = "Diversity", values_to = "Value")
        alpha_div_p <- compare_means(Value ~ Group, data = diversity_data,
                                     method = "wilcox", p.adjust.method = "holm",
                                     group.by = "Diversity")

        shannon_plot <- diversity_data %>% filter(Diversity == "Shannon") %>% 
          ggplot(aes(x = Group, y = Value, color = Group)) +
          geom_boxplot() +
          geom_jitter(shape = 16, position = position_jitter(0.2)) +
          scale_color_manual(values = c("#ffbf00", "#746AB0")) +
          theme_light() +
          labs(y = "ARGs diversity", x = "") +
          stat_pvalue_manual(subset(alpha_div_p, Diversity=="Shannon"), y.position = max(
            subset(diversity_data, Diversity=="Shannon")$Value) + 0.1, label = "p.adj") +
          guides(color = guide_legend(title = "Shannon", title.position = "top")) +
          theme(
            axis.title.x = element_text(size = 14, face = "bold"),
            axis.title.y = element_text(size = 14, face = "bold"),
            strip.text.x = element_text(size = 14, face = "bold"),
            axis.text.y = element_text(size = 14, face = "bold"),
            axis.text.x = element_text(size = 14, face = "bold"),
            legend.title = element_text(colour = "black",size = 10, face = "bold"),
            legend.text = element_text(colour = "black", size = 10, face = "bold"),
            plot.title = element_text(hjust = 0.5, size = rel(2)))
        
        simpson_plot <- diversity_data %>% filter(Diversity == "Simpson") %>% 
          ggplot(aes(x = Group, y = Value, color = Group)) +
          geom_boxplot() +
          geom_jitter(shape = 16, position = position_jitter(0.2)) +
          scale_color_manual(values = c("#ffbf00", "#746AB0")) +
          theme_light() +
          labs(y = "ARGs diversity", x = "") +
          stat_pvalue_manual(subset(alpha_div_p, Diversity=="Simpson"), y.position = max(
            subset(diversity_data, Diversity=="Simpson")$Value) + 0.01, label = "p.adj") +
          guides(color = guide_legend(title = "Simpson", title.position = "top")) +
          theme(
            axis.title.x = element_text(size = 14, face = "bold"),
            axis.title.y = element_text(size = 14, face = "bold"),
            strip.text.x = element_text(size = 14, face = "bold"),
            axis.text.y = element_text(size = 14, face = "bold"),
            axis.text.x = element_text(size = 14, face = "bold"),
            legend.title = element_text(colour = "black",size = 10, face = "bold"),
            legend.text = element_text(colour = "black", size = 10, face = "bold"),
            plot.title = element_text(hjust = 0.5, size = rel(2)))
        
        diversity_facet <- as_ggplot(grid.grabExpr(grid.arrange(shannon_plot,
                                                                      simpson_plot,
                                                                      ncol = 2)))
        
        return(diversity_facet)
    }
    fun_abundance_plot <- function(data)
    {
        abundance_data <- data %>% 
          group_by(Sample_Id, ARO_term) %>%
          summarise(Counts = sum(Normalized_counts))
        abundance_data <- inner_join(abundance_data, sample_metadata,
                                    by = "Sample_Id")
        abundance_data <- abundance_data %>% select(Sample_Id, Group, everything())
        abundance_data$Group <- factor(abundance_data$Group,
                                    levels = c("Control", "Case"))
        abundance_data <- subset(abundance_data, select = -Sample_Id)
        abundance_p_value <- compare_means(Counts~Group, data = abundance_data, method = "wilcox",
                                           p.adjust.method = "holm")
        abundance_plot <- ggplot(abundance_data, aes(x = Group, y = log2(Counts), color = Group)) +
          geom_boxplot() +
          stat_pvalue_manual(abundance_p_value, y.position = max(log2(abundance_data$Counts)) + 0.5,
                             label = "p.signif") +
          scale_color_manual(values = c("#ffbf00", "#746AB0")) +
          geom_jitter(shape = 16, position = position_jitter(0.2), alpha = 0.5) +
          theme_light() +
          xlab("Group") +
          ylab("Log2(Gene GCPM Counts)") +
          theme(
            axis.title.x = element_text(size = 14, face = "bold"),
            axis.title.y = element_text(size = 14, face = "bold"),
            strip.text.x = element_text(size = 14, face = "bold"),
            axis.text.y = element_text(size = 14, face = "bold"),
            axis.text.x = element_text(size = 14, face = "bold"),
            legend.title = element_text(colour = "black",size = 10, face = "bold"),
            legend.text = element_text(colour = "black", size = 10, face = "bold"),
            plot.title = element_text(hjust = 0.5, size = rel(2))
          )
        association_data <- data %>% group_by(Sample_Id, ARO_term) %>% 
          summarise(Summed_Counts = sum(Counts), length = max(ARG_length)/1000)
        
        association_data <- unique(association_data)
        
        association_data <- association_data %>% group_by(Sample_Id) %>% 
          mutate(RPKM_Counts = (Summed_Counts*1000000/(sum(Summed_Counts)*length)))
        
        association_data <- association_data %>% group_by(Sample_Id) %>% 
          summarise(RPKM_Abundance = sum(RPKM_Counts))
        
        association_data <- inner_join(association_data, sample_metadata,
                                     by = "Sample_Id")
        
        association_data$Group <- factor(association_data$Group,
                                       levels = c("Control", "Case"))
        
        association_p_value <- compare_means(RPKM_Abundance~Group, data = association_data, method = "wilcox",
                                           p.adjust.method = "holm")
        association_plot <- ggplot(association_data, aes(x = Group, y = log2(RPKM_Abundance), color = Group)) +
          geom_boxplot() +
          stat_pvalue_manual(association_p_value, y.position = max(log2(association_data$RPKM_Abundance)) + 0.1,
                             label = "p.signif") +
          geom_jitter(shape = 16, position = position_jitter(0.2)) +
          scale_color_manual(values = c("#ffbf00", "#746AB0")) +
          theme_light() +
          xlab("Group") +
          ylab("Log2(RPKM Total ARG Abundace)") +
          theme(
            axis.title.x = element_text(size = 14, face = "bold"),
            axis.title.y = element_text(size = 14, face = "bold"),
            strip.text.x = element_text(size = 14, face = "bold"),
            axis.text.y = element_text(size = 14, face = "bold"),
            axis.text.x = element_text(size = 14, face = "bold"),
            legend.title = element_text(colour = "black",size = 10, face = "bold"),
            legend.text = element_text(colour = "black", size = 10, face = "bold"),
            plot.title = element_text(hjust = 0.5, size = rel(2))
          )
        
        abundance_facet <- as_ggplot(grid.grabExpr(grid.arrange(association_plot,
                                                                abundance_plot,
                                                                ncol = 2)))
        
        return(abundance_facet)
    }
    return(list(diversity_plot = fun_diversity_plot(data),
        abundance_plot = fun_abundance_plot(data)))
})
observeEvent(input$upload_data, {
    plots_data <- alpha_diversity_analysis()
    output$plot_alpha_diversity <- renderPlot({
      plots_data$diversity_plot},height = 500
    )
    
  output$plot_abundance <- renderPlot({
        plots_data$abundance_plot}, height = 500)

    output$download_alpha_diversity_plot <- downloadHandler(
        filename = function() {
            paste("Alpha_diversity_plot_", Sys.Date(), ".png", sep = "")
        },
        content = function(file) {
            ggsave(file, plots_data$diversity_plot, width = 10.07, height = 3.96, units = "in", dpi = "retina")
        }
    )

    output$download_abundance_plot <- downloadHandler(
        filename = function() {
            paste("Abundance_plot_", Sys.Date(), ".png", sep = "")
        },
        content = function(file) {
            ggsave(file, plots_data$abundance_plot, width = 10.07, height = 3.96, units = "in", dpi = "retina")
        }
    )
})