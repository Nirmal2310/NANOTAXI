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
          selectInput("pipeline", "Select Analysis Tool", choices = list("KRAKEN2 + 16S NCBI DB"),
                      selected = "KRAKEN2 + 16S NCBI DB"),
          numericInput("min", "Minimum Length", value = 1400),
          numericInput("max", "Maximum Length", value = 1800),
          actionButton("start_analysis", "Start Analysis")
        ),
        conditionalPanel(
          condition = "input.data_file_type == 'precomputed'",
          textInput("fastqdir", "Enter the Directory Containing Raw Data", value = ""),
          fileInput("metafile", "Select the sample Information File", accept = ".csv", multiple = FALSE),
          selectInput("pipeline", "Select Analysis Tool", choices = list("BLASTn + 16s DB", "Kraken2 + Greengenes",
                                                             "Qiime2", "Centrifuge",
                                                             "EMU + Standard DB"),
                      selected = "BLASTn + 16s DB"),
          selectInput("tax", "Taxon Level", choices = list("Species", "Genus",
                                                              "Family", "Order"),
                      selected = "Species"),
          numericInput("min", "Minimum Length", value = 1400),
          numericInput("max", "Maximum Length", value = 1800),
          numericInput("threads", "Number of Threads", value = 16),
          numericInput("iden", "Percent Identity", value = 75),
          checkboxInput("setup", "Setup"),
          actionButton("upload_data", "Start Analysis")
        ),
        conditionalPanel(
          condition = "input.data_file_type == 'examplelist'",
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
