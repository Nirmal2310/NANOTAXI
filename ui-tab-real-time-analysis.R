tabPanel(
  div("Real-Time Analysis", style="font-size: 20px; font-weight: bold; color: #000000;
             font-family: serif"),
  fluidRow(id="realtime_analysis",
    column(id="realtime_data", 2,
           fluidRow(
             div(selectInput("select", "BARCODE", choices = list("barcode1"=1,   "barcode2"=2,   "barcode3"=3,   "barcode4"=4,   "barcode5"=5,   "barcode6"=6,
                                                                    "barcode7"=7,   "barcode8"=8,   "barcode9"=9,   "barcode10"=10, "barcode11"=11, "barcode12"=12,
                                                                    "barcode13"=13, "barcode14"=14, "barcode15"=15, "barcode16"=16, "barcode17"=17, "barcode18"=18,
                                                                    "barcode19"=19, "barcode20"=20, "barcode21"=21, "barcode22"=22, "barcode23"=23, "barcode24"=24),
                         selected = 1),
                 style="padding-left: 2%; text-align: center; width: 95%")
             ),
           tags$head(tags$style(HTML(".selectize-input {height: 40px}")))
           ),
    column(id="realtime_taxa", 2,
           div(selectInput("select", "TAXON", choices = list("Species"=1, "Genus"=2, "Family"=3,
                                                             "Order"=4), selected = 1),
               style="padding-left: 2%; height: 40px; text-align: center")
    ),
    column(id="realtime_status", 2,
           div(tags$b("SEQUENCING STATE"), style="padding-bottom: 4px; text-align: center"),
           verbatimTextOutput("state", placeholder = T),
           tags$head(
             tags$style(
               HTML(
               "#state {
        height: 40px;
        }"
               )))),
    column(id="realtime_status", 2,
           div(tags$b("N50 READ LENGTH"), style="padding-bottom: 4px; text-align: center"),
           verbatimTextOutput("n50", placeholder = T),
           tags$head(
             tags$style(
               HTML(
                 "#n50 {
        height: 40px;
        }")))),
    column(id="realtime_status", 2,
           div(tags$b("PASSED READS"), style="padding-bottom: 4px; text-align: center"),
           verbatimTextOutput("reads", placeholder = T),
           tags$head(
             tags$style(
               HTML(
                 "#reads {
        height: 40px;
        }")))),
    column(id="realtime_status", 2,
           div(tags$b("NUMBER OF SPECIES"), style="padding-bottom: 4px; text-align: center"),
           verbatimTextOutput("species", placeholder = T),
           tags$head(
             tags$style(
               HTML(
                 "#species {
        height: 40px;
        }"))))
  ),
  fluidRow(id="realtime_readlength_analysis",
           br(),
           column(id="realtime_readlength_plot", 6,
                  tags$b("READ LENGTH DISTRIBUTION"),
                           downloadButton(outputId = "download_readlength_plot",
                                          label = "Download Read Length Plot (PNG)"),
                           plotOutput("plot_read_classification_plot", height = "100%"))
           ,
           column(id="realtime_plot", 6,
                  tags$b("Q-SCORE HISTOGRAM"),
                                 downloadButton(outputId = "download_qscore_plot",
                                                label = "Download Q score Plot (PNG)"),
                                 plotOutput("plot_q_score_plot", height = "100%"))),
  br(),
  fluidRow(id="realtime_read_classification",
           column(id="realtime_readclassification_plot", 6,
                  tags$b("READS CLASSIFICATION"),
                           downloadButton(outputId = "download_classification_plot",
                                          label = "Download Reads Classification Plot (PNG)"),
                           plotOutput("plot_classification_plot",height = "100%", width = "100%")),
           column(id="realtime_sankey_plot", 6,
                  tags$b("SANKEY CLASSIFICATION"),
                  downloadButton(outputId = "download_sankey_plot",
                                 label = "Download Sankey Classification Plot (PNG)"),
                  plotOutput("plot_sankey_plot",height = "100%", width = "100%")))
)