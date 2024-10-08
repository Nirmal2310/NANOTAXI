tabPanel(
  div("HELP", style="font-size: 12px; font-weight: bold; color: #007B9A;
             font-family: serif"),
         fluidRow(
           column(2, wellPanel(
             h4("Additonal Information"),
             a("Input Data", href="#inputdata"),br(),
             a("Data Format", href="#dataformat"),br(),
             a("Output Metadata", href="#outputdata"),br(),
             a("Visulization", href="#vis"),br(),
             a("ARG Cohort Analysis", href="#cohortanalysis"),br(),
             a("ARG Distribution", href="#distributionplot"),br(),
             a("Alpha Diversity", href="#alphadiversity"),br(),
             a("Beta Diversity", href="#betadiversity"),br(),
           )
           ),
           column(8, includeMarkdown("Tabs/additional_information.md"))
         )
         
)
