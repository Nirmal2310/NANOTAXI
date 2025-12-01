tabPanel(
  div("HELP", style="font-size: 12px; font-weight: bold; color: #007B9A;
             font-family: serif"),
         fluidRow(
           column(2, wellPanel(
             h4("Additional Information"),
             a("Installation", href="#installation"),br(),
             a("Input Data", href="#inputdata"),br(),
             a("Pipeline", href="#algorithm"),br(),
             a("Data Format", href="#dataformat"),br(),
             a("Output Data", href="#outputdata"),br(),
             a("Visualization", href="#vis"),br(),
             a("Real-time Analysis", href="#realtime"),br(),
             a("Cohort Analysis", href="#cohortanalysis"),br(),
           )
           ),
           column(10, includeMarkdown("Tabs/additional_information.md"))
         )
         
)
