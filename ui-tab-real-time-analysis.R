tabPanel(
  title = div("REAL-TIME", style="font-size: 12px; font-weight: bold; color: #007B9A;
             font-family: serif"),
  value = "realtime_tab",
        fluidRow(id="realtime_analysis",
                 column(id="realtime_one", 6,
                        column(id="realtime_barcode", 3,
                               fluidRow(
                                div(selectInput("barcode_select", "BARCODE", choices = paste0("barcode",sprintf("%02d", 1:24)),
                                selected = "barcode01"),
                                style="padding-left: 2%; text-align: center; width: 95%; font-size: 13px; font-weight: bold; color: #00607D;"),
                                tags$head(
                                  tags$style(
                                    HTML(
                                      "#barcode_select {color: #00607D;}"
                                    )
                                  )
                                )
                              ),
                              tags$head(tags$style(HTML(".selectize-input {height: 40px}")))
                              ),
                              column(id="realtime_taxa", 3,
                              div(selectInput("taxon_select", "TAXON", choices = list("Species", "Genus", "Family",
                              "Order", "Class", "Phylum", "Kingdom"), selected = "Species"),
                              style="padding-left: 2%; height: 40px; text-align: center; font-size: 13px; font-weight: bold; color: #00607D;"),
                              tags$head(
                                tags$style(
                                  HTML(
                                    "#taxon_select {color: #00607D;}"
                                  )
                                )
                              )
                            ),
                            column(id="realtime_status", 3,
                            div(tags$b("STATE"), style="padding-bottom: 4px; text-align: center; color: #00607D;"),
                            verbatimTextOutput("state", placeholder = T),
                            tags$head(
                              tags$style(
                                HTML(
                                  "#state {height: 40px; font-size: 13px; text-align: center; font-weight: bold}")
                                  )
                                  )
                            ),
                            column(id="realtime_pores", 3,
                            div(tags$b("PORES"), style="padding-bottom: 4px; text-align: center; color: #00607D;"),
                            verbatimTextOutput("pores", placeholder = T),
                            tags$head(
                              tags$style(
                                HTML("#pores {height: 40px; font-size: 13px; text-align: center; font-weight: bold}"
                                )
                                )
                                )
                              )
                          ),
                          column(id="realtime_two", 6,
                          column(id="realtime_reads", 3,
                          div(tags$b("READS"), style="padding-bottom: 4px; text-align: center; color: #00607D;"),
                          verbatimTextOutput("passed_reads", placeholder = T),
                          tags$head(
                            tags$style(
                              HTML(
                                "#passed_reads {height: 40px; font-size: 13px; text-align: center; font-weight: bold}"
                                )
                                )
                                )
                              ),
                              column(id="realtime_mean", 3,
                              div(tags$b("MEAN LENGTH"), style="padding-bottom: 4px; text-align: center; color: #00607D;"),
                              verbatimTextOutput("mean_length", placeholder = T),
                              tags$head(
                                tags$style(
                                  HTML(
                                    "#mean_length {height: 40px; font-size: 13px; text-align: center; font-weight: bold;}"
                                    )
                                    )
                                    )
                                ),
                                column(id="realtime_classified", 3,
                                div(tags$b("CLASSIFIED"), style="padding-bottom: 4px; text-align: center; color: #00607D;"),
                                verbatimTextOutput("classified_reads", placeholder = T),
                                tags$head(
                                  tags$style(
                                    HTML(
                                     "#classified_reads {height: 40px; font-size: 13px; text-align: center; font-weight: bold;}"
                                     )
                                     )
                                     )
                                ),
                                column(id="realtime_species", 3,
                                div(tags$b("TAXON COUNT"), style="padding-bottom: 4px; text-align: center; color: #00607D;"),
                                verbatimTextOutput("taxon_count", placeholder = T),
                                tags$head(
                                  tags$style(
                                    HTML("#taxon_count {height: 40px; font-size: 13px; text-align: center; font-weight: bold; color: #00607D;}")))))),
                                    br(),
                                    fluidRow(id="realtime_readlength_analysis",
                                    column(width = 12,
                                    tabsetPanel(id="realtime_plots",
                                    tabPanel(title = div("TAXON TABLE", style="font-size: 15px; font-weight: bold; color: #00607D; font-family: serif"),
                                    div(downloadButton(outputId = "download_taxa_table",
                                    label = "Download Table (TSV)"), style = "padding-top: 10px;"),
                                    DT::dataTableOutput("taxa_table", height = "100%", width = "50%")),
                                    tabPanel(title = div("READ LENGTH", style="font-size: 15px; font-weight: bold; color: #00607D; font-family: serif"),
                                    div(downloadButton(outputId = "download_readlength_plot",
                                    label = "Download Plot (PNG)"), style = "padding-top: 10px;"),
                                    plotOutput("plot_read_length", height = "100%", width = "100%")
                                    ),
                                    tabPanel(title = div("Q-SCORE", style="font-size: 15px; font-weight: bold; color: #00607D; font-family: serif"),
                                    div(downloadButton(outputId = "download_qscore_plot",
                                    label = "Download Plot (PNG)"), style = "padding-top: 10px;"),
                                    plotOutput("plot_q_score_plot", height = "100%", width = "75%")
                                    ),
                                    tabPanel(title = div("CLASSIFICATION", style="font-size: 15px; font-weight: bold; color: #00607D; font-family: serif"),
                                    fluidRow(column(width = 2,
                                    sliderInput("top_n", "Top N", min = 1, max = 25, value = 10)),
                                    column(width=10,
                                    div(downloadButton(outputId = "download_classification_plot",
                                    label = "Download Plot (PNG)"), style = "float: right; padding-right: 15px; padding-top: 10px;"),
                                    plotOutput("plot_classification_plot", height = "100%")))
                                    ))))
    )
