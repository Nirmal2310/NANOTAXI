tabPanel(
div("INTRODUCTION", style="font-size: 12px; font-weight: bold; color: #007B9A;
             font-family: serif"),
         fluidRow(
           column(2, wellPanel(
             h4("Introduction"),
             a("Features", href="#Features"),br(),
             a("Data Formats", href="#DataFormat"),br(),
             a("Additional Information", href="#help"),br()
           )
           ),
           column(8, includeMarkdown("Tabs/Instructions.md"))
         ),
         p(hr(), p(("ShinyApp created by Nirmal Singh Mahar and Ishaan Gupta*"), align = "center", width=2),
           p(("Copyrigth (C) 2024, code licensed under GPLv3"), align="center", width=2),
           p(("Code available on Github:"), a("https://github.com/Nirmal2310/NanoTAXI",
           href="https://github.com/Nirmal2310/NanoTaxi"), align="center",width=2)
         )
  
)
