tabPanel(
  # div("Input Data", style="font-size: 20px; font-weight: bold; color: #000000;
  #            font-family: serif")
  "Input Data",
  fluidRow(
    column(
      width = 3,
      wellPanel(
        radioButtons(
          "data_file_type",
          "Use example data or upload the list",
          c("Example Data" = "examplelist","Cohort Analysis" = "precomputed",
            "Real-Time Analysis" = "upload"),
          selected = "examplelist"
        ),
        conditionalPanel(
          condition = "input.data_file_type == 'upload'",
          fileInput("inputfile", "Select the sample Information File", accept = ".csv", multiple = FALSE),
          selectInput("pipeline", "Select Analysis Tool", choices = list("KRAKEN2 + 16S NCBI DB"),
                      selected = "KRAKEN2 + 16S NCBI DB"),
          numericInput("num", "Minimum Length", value = 1400),
          numericInput("num", "Maximum Length", value = 1800),
          checkboxInput("setup", "Setup"),
          actionButton("real_time","Start Analysis")
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
                                                              "Family", "Order", "Class", "Phylum"),
                      selected = "Species"),
          numericInput("min", "Minimum Length", value = 1400),
          numericInput("max", "Maximum Length", value = 1800),
          numericInput("threads", "Number of Threads", value = 16),
          
          checkboxInput("setup", "Setup"),
          actionButton("upload_data","Start Analysis")
        ),
        conditionalPanel(
          condition = "input.data_file_type == 'examplelist'",
          actionButton("upload_data","Use Example Data")
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
                        downloadButton("download_results_CSV","Save Results as CSV File"),
                        DTOutput("analysisoutput")
                        )
      )
    )
  ),
  p(hr(), p(("ShinyApp created by Nirmal Singh Mahar and Ishaan Gupta*"), align = "center", width=2),
    p(("Copyrigth (C) 2023, code licensed under GPLv3"), align="center", width=2),
    p(("Code available on Github:"), a("https://github.com/Nirmal2310/NanoTAXI",
                                       href="https://github.com/Nirmal2310/NanoTAXI"),
      align="center",width=2)
  )
)


