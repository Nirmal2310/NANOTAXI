tabPanel(
  div("INPUT", style="font-size: 12px; font-weight: bold; color: #007B9A;
             font-family: serif"),
  fluidRow(
    column(
      width = 3,
      wellPanel(
        radioButtons(
          "data_file_type",
          "Use Example or Upload List",
          c("Example Data" = "examplelist", "Offline Analysis" = "precomputed",
            "Real-Time Analysis" = "upload"),
          selected = "examplelist"
        ),
        conditionalPanel(
          condition = "input.data_file_type == 'upload'",
          fileInput("inputfile", "Select the sample Information File", accept = ".csv", multiple = FALSE),
          textInput("realtime_control", "Select Control Group", value = ""),
          selectInput("realtime_pipeline", "Select Analysis Tool", choices = list("Kraken2", "Minimap2", "BLASTn", "EMU"),
                      selected = "Minimap2"),
          selectInput("realtime_database", "Select Database", choices = list("GTDB", "MIMT", "GSR", "REFSEQ", "EMUDB"),
                      selected = "REFSEQ"),
          selectInput("realtime_tax", "Taxon Level", choices = list("Species", "Genus",
                                                              "Family", "Order"), selected = "Species"),
          numericInput("realtime_threads", tags$span("Number of Threads", icon("info-circle", id = "info_threads",
                        style = "cursor: pointer; color: #337ab7;")), value = 24),
          bsTooltip("info_threads", "Maximum number of threads available. NANOTAXI will distribute the computational threads based on the barcodes available for processing.", placement = "right", trigger = "hover"),
          numericInput("realtime_min", "Minimum Length", value = 1400),
          numericInput("realtime_max", "Maximum Length", value = 1800),
          numericInput("realtime_q_score", "Q-Score", value = 10),
          sliderInput("realtime_conf", "Kraken Confidence Score", min = 0.00, max = 1.00, value = 0.00),
          numericInput("realtime_iden", "Percent Identity", value = 85),
          numericInput("realtime_cov", "Percent Coverage", value = 85),
          sliderInput("chunk_size", tags$span("Batch Size", icon("info-circle", id = "info_chunk", style = "cursor: pointer; color: #337ab7;")), 
                      min = 100, max = 4000, value = 500, step = 100, ticks = FALSE),
          bsTooltip("info_chunk", "Number of reads processed per barcode in each iteration.", placement = "right", trigger = "hover"),
          selectInput("refresh_rate", tags$span("Update Interval", icon("info-circle", id = "info_interval", style = "cursor: pointer; color: #337ab7;")), 
              choices = c("1 second" = 1, "10 seconds" = 10, "30 seconds" = 30, "60 seconds" = 120), selected = 10),
          bsTooltip("info_interval", "Buffer time between each iteration. ", placement = "right", trigger = "hover"),
          actionButton("start_analysis", "Start Analysis")
        ),
        conditionalPanel(
          condition = "input.data_file_type == 'precomputed'",
          textInput("fastqdir", "Enter the Directory Containing Raw Data", value = ""),
          fileInput("metafile", "Select the sample Information File", accept = ".csv", multiple = FALSE),
          textInput("control", "Select Control Group", value = ""),
          selectInput("kitname", "Kit Name", choices = list("SQK-16S024", "SQK-16S114-24"), selected = "SQK-16S114-24"),
          selectInput("pipeline", "Select Analysis Tool", choices = list("EMU", "BLASTn", "Kraken2", "Minimap2"),
                      selected = "EMU"),
          selectInput("database", "Select Database", choices = list("GTDB", "MIMT", "GSR", "REFSEQ", "EMUDB"),
                      selected = "EMUDB"),
          selectInput("tax", "Taxon Level", choices = list("Species", "Genus",
                                                              "Family", "Order"),
                      selected = "Species"),
          numericInput("min", "Minimum Length", value = 1400),
          numericInput("max", "Maximum Length", value = 1800),
          numericInput("q_score", "Q-Score", value = 10),
          sliderInput("conf", "Kraken Confidence Score", min = 0.00, max = 1.00, value = 0.00),
          numericInput("threads", "Number of Threads", value = 16),
          numericInput("iden", "Percent Identity", value = 85),
          numericInput("cov", "Percent Coverage", value = 85),
          checkboxInput("setup", "Setup"),
          actionButton("upload_data", "Start Analysis")
        ),
        conditionalPanel(
          condition = "input.data_file_type == 'examplelist'",
          textInput("control", "Select Control Group", value = ""),
          actionButton("example_run", "Use Example Data")
        )
    )),

    column(
      width = 9,
      bsCollapse(id = "input_collapse_panel",open = "data_panel",multiple = FALSE,
        bsCollapsePanel(title = "Data Contents: Check Before `Submit`",
          value = "data_panel",
          DTOutput("sampleinfo")
        ),
        bsCollapsePanel(title = "Analysis Results: Ready to View Other Tabs",
                        value = "analysis_panel",
                        downloadButton("download_results_csv","Save Results as CSV File"),
                        DTOutput("analysisoutput")
                        )
      )
    )
  ),
  p(hr(), p(("ShinyApp created by Nirmal Singh Mahar and Ishaan Gupta*"), align = "center", width=2),
    p(("Copyright (C) 2024, code licensed under GPLv3"), align="center", width=2),
    p(("Code available on GitHub:"), a("https://github.com/Nirmal2310/NANOTAXI",
                                       href="https://github.com/Nirmal2310/NANOTAXI"),
      align="center",width=2)
  )
)