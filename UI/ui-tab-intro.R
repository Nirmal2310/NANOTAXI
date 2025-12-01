tabPanel(
  div("INTRODUCTION", style="font-size: 12px; font-weight: bold; color: #007B9A;
             font-family: serif"),
  fluidRow(
    column(2, wellPanel(
      h4("Introduction"),
      a("Features", href="#Features"),br(),
      a("Data Visualization", href="#Visualization"),br(),
      a("Data Format", href="#DataFormat"),br(),
      a("Input Data", href="#InputData"),br(),
      a("Output Data", href="#Output"),br(),
      a("Additional Information", href="#Help"),br(),
      a("Author Details", href="#Developer"),br(),
      a("Citation Details", href="#Citation")
    )
    ),
    column(10, includeMarkdown("Tabs/Instructions.md"))
  )
  
)
