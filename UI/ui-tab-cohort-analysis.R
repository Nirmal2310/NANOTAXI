tabPanel(
  title = div("COHORT", style="font-size: 12px; font-weight: bold; color: #007B9A;
             font-family: serif"),
  value = "cohort_tab",
    fluidRow(
      column(
        width = 2,
        div(selectInput("taxa", "Taxon Level", choices = list("Species", "Genus",
                                                          "Family", "Order", "Class", "Phylum"),
                    selected = "Species"), style="color: #00607D;"),
        div(sliderInput("top_taxa", "Top N", min = 1, max = 25, value = 5), style="color: #00607D;"),
        div(selectizeInput("toi", "Select Taxon Name", choices = NULL, selected = NULL, multiple = TRUE, options = list(searchConjunction = 'and')), style="color: #00607D;"),
        div(sliderInput("prevalence_cutoff", "Prevalence Cutoff (%)", min = 5, max = 100, value = 10), style="color: #00607D;"),
        div(numericInput("abundance_cutoff", "Relative Abundance Cutoff (%)", value = 0.1), style="color: #00607D;"),
        div(numericInput("counts_cutoff", "Absolute Counts Cutoff", value = 10), style="color: #00607D;"),
        div(bslib::input_switch("show_names", "Show Barcodes Name"), style="font-size: 15px; font-weight: bold; color: #00607D"),
        div(sliderInput("biplot_taxa", "Biplot Top Taxa", min = 5, max = 25, value = 10), style="color: #00607D;"),
        div(selectInput("fun_cat", "Functional Pathway Group", choices = list("Kegg Orthologs (KO)", "MetaCyc Pathway"),
                        selected = "Kegg Ontology (KO)"), style="color: #00607D;"),
        div(sliderInput("path_count", "DA Functional Categories", min = 6, max = 20, value = 10, step = 2), style="color: #00607D;")
      ),
      column(
        width = 10,
        tabsetPanel(id="cohort_analysis_plots",
                    tabPanel(title = div("Taxa Stacked Barplot", style="font-size: 12px; font-weight: bold; color: #00607D;
                                                    font-family: serif"),
                              value = "cohort_stacked_plot",
                             div(downloadButton(outputId = "download_stacked_barplot",
                                            label = "Download Stacked Bar Plot (PDF)"), style="margin-top: 10px;"),
                             plotlyOutput("plot_stacked_barplot", height = "100%", width = "100%")
                    ),
                    tabPanel(title = div("Alpha Diversity Plot", style="font-size: 12px; font-weight: bold; color: #00607D;
                                                    font-family: serif"),
                             div(downloadButton(outputId = "download_boxplot",
                                            label = "Download Boxplot (PDF)"), style="margin-top: 10px;"),
                             plotOutput("plot_boxplot", height = "100%")
                             
                    ),
                    tabPanel(title = div("NMDS Plot", style="font-size: 12px; font-weight: bold; color: #00607D;
                                                    font-family: serif"),
                             div(downloadButton(outputId = "download_nmds_plot",
                                            label = "Download NMDS Plot (PDF)"), style="margin-top: 10px;"),
                              plotOutput("plot_nmds", height = "100%")
                    ),
                    tabPanel(title = div("PCoA Plot", style="font-size: 12px; font-weight: bold; color: #00607D;
                                                    font-family: serif"),
                             div(downloadButton(outputId = "download_pcoa_plot",
                                            label = "Download PCoA Plot (PDF)"), style="margin-top: 10px;"),
                             plotOutput("plot_pcoa", height = "100%")
                    ),
                    tabPanel(title = div("PCA Plot", style="font-size: 12px; font-weight: bold; color: #00607D;
                                                    font-family: serif"),
                             div(downloadButton(outputId = "download_pca_plot",
                                            label = "Download PCA Plot (PDF)"), style="margin-top: 10px;"),
                             plotOutput("plot_pca", height = "100%")
                    ),
                    tabPanel(title = div("PERMANOVA", style="font-size: 12px; font-weight: bold; color: #00607D;
                                                    font-family: serif"),
                             div(downloadButton(outputId = "download_permanova_csv",
                                            label = "Download Result (TSV)"), style="margin-top: 10px;"),
                             DT::dataTableOutput("permanova_data"), height = "100%", width = "50%"),
                    tabPanel(title = div("HeatMap", style="font-size: 12px; font-weight: bold; color: #00607D;
                                                    font-family: serif"),
                             div(downloadButton(outputId = "download_heatmap",
                                            label = "Download HeatMap (PDF)"), style="margin-top: 10px;"),
                             div(plotOutput("plot_heatmap", height = "100%"), 
                                style = "margin-bottom: 10px;")
                    ),
                    tabPanel(title = div("DAA", style="font-size: 12px; font-weight: bold; color: #00607D;
                                                    font-family: serif"),
                                      div(div(selectizeInput("daa_comp", "", choices = NULL, selected = NULL, multiple = FALSE, width = '20%'),
                                          style="color: #00607D; display: flex; justify-content: flex-end;"),
                                          div(downloadButton(outputId = "download_daa",
                                          label = "Download Volcano Plot (PDF)"),
                                          downloadButton(outputId = "download_daa_csv",
                                          label = "Download Result (TSV)"), style="margin-top: -50px;")
                                        ),
                                      div(plotOutput("plot_volcano", height = "100%"), style="margin: 15px 0 0 0;")
                    ),
                    tabPanel(title =  div("PCA(Functional)", style="font-size: 12px; font-weight: bold; color: #00607D;
                                                    font-family: serif"),
                                      div(style="margin-top: 10px;",
                                      downloadButton(outputId = "download_functional_pca",
                                      label = "Download Pathway PCA Plot (PDF)")),
                                      plotOutput("plot_functional_pca", height = "100%")
                    ),
                    tabPanel(title =  div("DAA(Functional)", style="font-size: 12px; font-weight: bold; color: #00607D;
                                                    font-family: serif"),
                                      div(div(selectizeInput("daa_fun_com", "", choices = NULL, selected = NULL, multiple = FALSE, width = '20%'), 
                                          style="color: #00607D; display: flex; justify-content: flex-end;"),
                                          div(downloadButton(outputId = "download_functional_DAA",
                                          label = "Download Enrichment Plot (PDF)"),
                                          downloadButton(outputId = "download_functional_daa_csv",
                                          label = "Download Result (TSV)"), style="margin-top: -50px;")
                                      ),
                                      div(plotOutput("plot_functional_daa", height = "100%"), style="margin: 15px 0 0 0;")
                    )
                    
        )
        
      )
    )
)