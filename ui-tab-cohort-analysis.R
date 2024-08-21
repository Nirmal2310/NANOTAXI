tabPanel(
  "Cohort Analysis",
  fluidRow(
    column(
      width = 2,
      selectInput("taxa", "Taxon Level", choices = list("Species", "Genus",
                                                          "Family", "Order"),
                  selected = "Species")
    ),
    column(
      width = 10,
      tabsetPanel(id="cohort_analysis_plots",
                  tabPanel(title = div("Taxa Stacked Barplot", style="font-size: 15px; font-weight: bold; color: #000000;
                                                    font-family: serif"),
                           downloadButton(outputId = "download_stacked_barplot",
                                          label = "Download Stacked Bar Plot (PNG)"),
                           plotOutput("plot_stacked_barplot", height = "100%")
                  ),
                  tabPanel(title = div("Alpha Diversity Plot", style="font-size: 15px; font-weight: bold; color: #000000;
                                                    font-family: serif"),
                           downloadButton(outputId = "download_boxplot",
                                          label = "Download Boxplot (PNG)"),
                           plotOutput("plot_boxplot", height = "100%")
                           
                  ),
                  tabPanel(title = div("PCoA Plot", style="font-size: 15px; font-weight: bold; color: #000000;
                                                    font-family: serif"),
                           downloadButton(outputId = "download_pcoa_plot",
                                          label = "Download PCoA Plot (PNG)"),
                           plotOutput("plot_pcoa", height = "100%")
                  ),
                  tabPanel(title = div("NMDS Plot", style="font-size: 15px; font-weight: bold; color: #000000;
                                                    font-family: serif"),
                           downloadButton(outputId = "download_nmds_plot",
                                          label = "Download NMDS Plot (PNG)"),
                           plotOutput("plot_nmds", height = "100%")),
                  tabPanel(title = div("PCA Plot", style="font-size: 15px; font-weight: bold; color: #000000;
                                                    font-family: serif"),
                           downloadButton(outputId = "download_pca_plot",
                                          label = "Download PCA Plot (PNG)"),
                           plotOutput("plot_pca", height = "100%")
                  ),
                  tabPanel(title = div("HeatMap", style="font-size: 15px; font-weight: bold; color: #000000;
                                                    font-family: serif"),
                           downloadButton(outputId = "download_heatmap",
                                          label = "Download HeatMap (PNG)"),
                           plotOutput("plot_heatmap", height = "100%")
                  )
                  
      )
      
    )
  )
)