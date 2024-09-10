# library(shiny)
# library(reticulate)  # For running Python scripts
# library(reactable)   # For interactive tables
# library(plotly)      # For interactive plots
# 
# # UI
# ui <- fluidPage(
#   tabPanel(
#         "Real-Time Analysis",
#         fluidRow(id="realtime_analysis",
#                  column(id="realtime_one", 6,
#                         column(id="realtime_barcode", 3,
#                                fluidRow(
#                                 div(selectInput("barcode_select", "BARCODE", choices = paste0("barcode",sprintf("%02d", 1:24)),
#                                 selected = "barcode01"),
#                                 style="padding-left: 2%; text-align: center; width: 95%")
#                               ),
#                               tags$head(tags$style(HTML(".selectize-input {height: 40px}")))
#                               ),
#                               column(id="realtime_taxa", 3,
#                               div(selectInput("select", "TAXON", choices = list("Species", "Genus", "Family",
#                               "Order"), selected = "Species"),
#                               style="padding-left: 2%; height: 40px; text-align: center")
#                             ),
#                             column(id="realtime_status", 3,
#                             div(tags$b("STATE"), style="padding-bottom: 4px; text-align: center"),
#                             verbatimTextOutput("state", placeholder = T),
#                             tags$head(
#                               tags$style(
#                                 HTML(
#                                   "#state {height: 40px; font-size: 13px; text-align: center; font-weight: bold;}")
#                                   )
#                                   )
#                             ),
#                             column(id="realtime_pores", 3,
#                             div(tags$b("PORES"), style="padding-bottom: 4px; text-align: center"),
#                             verbatimTextOutput("pores", placeholder = T),
#                             tags$head(
#                               tags$style(
#                                 HTML("#pores {height: 40px; font-size: 13px; text-align: center; font-weight: bold;}"
#                                 )
#                                 )
#                                 )
#                               )
#                           ),
#                           column(id="realtime_two", 6,
#                           column(id="realtime_reads", 3,
#                           div(tags$b("PASSED READS"), style="padding-bottom: 4px; text-align: center"),
#                           verbatimTextOutput("passed_reads", placeholder = T),
#                           tags$head(
#                             tags$style(
#                               HTML(
#                                 "#passed_reads {height: 40px; font-size: 13px; text-align: center; font-weight: bold;}"
#                                 )
#                                 )
#                                 )
#                               ),
#                               column(id="realtime_mean", 3,
#                               div(tags$b("MEAN LENGTH"), style="padding-bottom: 4px; text-align: center"),
#                               verbatimTextOutput("mean_length", placeholder = T),
#                               tags$head(
#                                 tags$style(
#                                   HTML(
#                                     "#mean_length {height: 40px; font-size: 13px; text-align: center; font-weight: bold;}"
#                                     )
#                                     )
#                                     )
#                                 ),
#                                 column(id="realtime_classified", 3,
#                                 div(tags$b("CLASSIFIED"), style="padding-bottom: 4px; text-align: center"),
#                                 verbatimTextOutput("classified_reads", placeholder = T),
#                                 tags$head(
#                                   tags$style(
#                                     HTML(
#                                      "#classified_reads {height: 40px; font-size: 13px; text-align: center; font-weight: bold;}"
#                                      )
#                                      )
#                                      )
#                                 ),
#                                 column(id="realtime_species", 3,
#                                 div(tags$b("# SPECIES"), style="padding-bottom: 4px; text-align: center"),
#                                 verbatimTextOutput("species", placeholder = T),
#                                 tags$head(
#                                   tags$style(
#                                     HTML("#species {height: 40px; font-size: 13px; text-align: center; font-weight: bold;}"
#                                     )
#                                     )
#                                     )
#                                 )
#                         )
#         )
#         # ,
#         # fluidRow(id="realtime_readlength_analysis",
#         #          br(),
#         #          column(id="realtime_readlength_plot", 6,
#         #                 tags$b("READ LENGTH DISTRIBUTION"),
#         #                 downloadButton(outputId = "download_readlength_plot",
#         #                                label = "Download Read Length Plot (PNG)"),
#         #                 plotOutput("plot_read_classification_plot", height = "100%"))
#         #          ,
#         #          column(id="realtime_plot", 6,
#         #                 tags$b("Q-SCORE HISTOGRAM"),
#         #                 downloadButton(outputId = "download_qscore_plot",
#         #                                label = "Download Q score Plot (PNG)"),
#         #                 plotOutput("plot_q_score_plot", height = "100%"))),
#         # br(),
#         # fluidRow(id="realtime_read_classification",
#         #          column(id="realtime_readclassification_plot", 6,
#         #                 tags$b("READS CLASSIFICATION"),
#         #                 downloadButton(outputId = "download_classification_plot",
#         #                                label = "Download Reads Classification Plot (PNG)"),
#         #                 plotOutput("plot_classification_plot",height = "100%", width = "100%")),
#         #          column(id="realtime_sankey_plot", 6,
#         #                 tags$b("SANKEY PLOT"),
#         #                 downloadButton(outputId = "download_sankey_plot",
#         #                                label = "Download Sankey Classification Plot (PNG)"),
#         #                 plotOutput("plot_sankey_plot",height = "100%", width = "100%")))
#       ),
#   actionButton("start_analysis", "Start Analysis"),
#   # verbatimTextOutput("sequencing_status", placeholder = T),
#   # textOutput("active_pores"),
#   # selectInput("barcode_select", "Select Barcode", choices = paste0("barcode",sprintf("%02d", 1:9))),
#   # tableOutput("passed_df"),
#   # plotOutput("analysis_plot"),
#   # textOutput("per_barcode_data")
# )
# 
# # Server
# server <- function(input, output, session) {
# 
#   state <- reactiveVal("Initializing")
# 
#   is_running <- reactiveVal(FALSE)
# 
#   timer_30s <- reactiveTimer(1000)  # 30 seconds
# 
#   timer_2m <- reactiveTimer(1000)  # 2 minutes
# 
#   timer_5m <- reactiveTimer(1000)  # 5 minutes
# 
#   pores <- reactiveVal()
# 
#   passed_counts <- reactiveVal()
# 
#   mean_df <- reactiveVal()
# 
#   classified_df <- reactiveVal()
# 
#   species_df <- reactiveVal()
# 
#   new_files <- reactiveVal(character(0))
# 
#   observeEvent(input$start_analysis, {
#     is_running(!is_running())
#     updateActionButton(session, "start_analysis",
#                        label = if (is_running()) "Stop Analysis" else "Start Analysis")
#   })
# 
#   observe({
#     if (is_running()) {
#       timer_30s()
#       check_sequencing_status()
#     }
#   })
# 
#   observe({
#     if (is_running()) {
#       timer_2m()
#       run_analysis_scripts()
#     }
#   })
# 
#   # observe({
#   #   if (is_running()) {
#   #     timer_5m()
#   #     check_and_analyze_files()
#   #   }
#   # })
# 
#   # Function to run Python script and get sequencing status
#   check_sequencing_status <- function() {
#     # Replace with your actual Python script execution
# 
#     path <- paste0(getwd(),"/Scripts")
# 
#     script <- paste0(path,"/state_dummy.py")
# 
#     status <- py_run_file(script)
# 
#     state(status$status)
#   }
# 
#   run_analysis_scripts <- function() {
# 
#     path <- paste0(getwd(),"/Scripts")
# 
#     script_1 <- paste0(path,"/pores_dummy.py")
# 
#     script_2 <- paste0(path,"/barcode_dummy.py")
# 
#     script_3 <- paste0(path,"/n50_dummy.py")
# 
#     script_4 <- paste0(path,"/classified_reads_dummy.py")
# 
#     script_5 <- paste0(path,"/unique_taxa_counts_dummy.py")
# 
#     pores(py_run_file(script_1)$pore_counts)
# 
#     passed_counts(py_run_file(script_2)$df)
# 
#     mean_df(py_run_file(script_3)$df)
# 
#     classified_df(py_run_file(script_4)$df)
# 
#     species_df(py_run_file(script_5)$df)
# 
#   }
# 
#   # Display sequencing status
#   output$state <- renderText({
#     state()
#   })
# 
#   output$pores <- renderText({
#     pores()
#   })
# 
#   # output$passed_df <- renderTable({
#   #   passed_counts()
#   # })
# 
#   output$passed_reads <- renderText({
#     req(input$barcode_select, passed_counts())
#     selected_count <- passed_counts()[passed_counts()$Barcode == input$barcode_select, "Counts"]
#   })
# 
#   output$mean_length <- renderText({
#     req(input$barcode_select, mean_df())
#     mean_length <- mean_df()[mean_df()$Barcode == input$barcode_select, "Mean_Length"]
#   })
# 
#   output$classified_reads <- renderText({
#     req(input$barcode_select, classified_df())
#     classified_reads <- classified_df()[classified_df()$Barcode == input$barcode_select, "Classified"]
#   })
# 
#   output$species <- renderText({
#     req(input$barcode_select, species_df())
#     unique_taxa_counts <- species_df()[species_df()$Barcode == input$barcode_select, "Counts"]
#   })
# 
# 
# 
# 
# }
# 
# # Run the app
# shinyApp(ui = ui, server = server)
# 
# # Set up Python environment
# 
# # conda_path <- system('which conda | sed "s/condabin.*$//g"', intern = TRUE)
# #
# # use_python(paste0(conda_path,"envs/realtime/bin/python"))
# 
# # Define the server logic
# # server <- function(input, output, session) {
# #
# #   # Variables to store reactive values
# #
# #   workdir <- getwd()
# #
# #   script_path <- paste0(workdir,"/Scripts")
# #
# #   status <- reactiveVal(NULL)
# #
# #   active_pores <- reactiveVal(NULL)
# #
# #   passed_counts <- reactiveVal(NULL)
# #
# #   barcodes <- reactiveVal(NULL)
# #
# #   # Start analysis when button is clicked
# #   observeEvent(input$start, {
# #     # Schedule the first script to check status every 30 seconds
# #     invalidateLater(100, session)
# #     run_status_script(script_path)
# #   })
# #
# #   # Function to run the first Python script
# #   run_status_script <- function(path) {
# #
# #     # Run the first Python script
# #
# #     script <- paste0(path, "/state_dummy.py")
# #
# #     status_script <- py_run_file(script)
# #
# #     status(status_script$status)
# #
# #     if (status() == "Active") {
# #       # Schedule the second and third scripts every 2 minutes
# #       invalidateLater(100, session)
# #       run_second_script(script_path)
# #       run_third_script(script_path)
# #     }
# #   }
# #
# #   # Function to run the second Python script (Active Pores)
# #
# #   run_second_script <- function(path) {
# #
# #     script <- paste0(path, "/pores_dummy.py")
# #
# #     pores_script <- py_run_file(script)
# #
# #     active_pores(pores_script$pore_counts)
# #   }
# #
# #   # Function to run the third Python script (Passed Counts)
# #   run_third_script <- function(path) {
# #
# #     script <- paste0(path, "/barcode_dummy.py")
# #
# #     counts_script <- py_run_file(script)
# #
# #     passed_counts(as.data.table(counts_script$df))
# #
# #     print(as.data.table(counts_script$df))
# #   }
# #
# #   # Function to check and analyze new FASTQ files
# #   # analyze_fastq_files <- function() {
# #   #   invalidateLater(300000, session)
# #   #   fastq_script <- py_run_file("analyze_fastq_script.py")
# #   #   # Assuming the script returns plots, you can render them here
# #   #   output$plot1 <- renderPlot({ fastq_script$plot1 })
# #   #   output$plot2 <- renderPlot({ fastq_script$plot2 })
# #   # }
# #
# #   # Reactive display of sequencing status
# #   output$status <- renderText({
# #     paste("Sequencing Status:", status())
# #   })
# #
# #   # Reactive display of active pores
# #   output$active_pores <- renderText({
# #     paste("Active Pores:", active_pores())
# #   })
# #
# #   # Reactive display of passed counts for the selected barcode
# #   # output$passed_reads <- renderText({
# #   #   req(input$barcode)
# #   #   counts <- passed_counts()[Barcode == input$barcode, Passed_Counts]
# #   #   paste("Passed Counts for", input$barcode, ":", counts)
# #   # })
# #
# #   # Start the file analysis if the status is active
# #   # observe({
# #   #   if (status() == "Active") {
# #   #     analyze_fastq_files()
# #   #   }
# #   # })
# # }
# #
# # # Run the application
# # shinyApp(ui = ui, server = server)


library(shiny)
library(tidyverse)

df <- structure(list(site = structure(c(1L, 2L, 4L, 5L, 6L, 7L, 8L, 
                                        9L, 10L, 11L, 3L), .Label = c("allsites", "site1", "site10", 
                                                                      "site2", "site3", "site4", "site5", "site6", "site7", "site8", 
                                                                      "site9"), class = "factor"), lower = c(40.97, 1.92, 7.5, 2.66, 
                                                                                                             1.18, 0.72, 6.92, 6.87, 3.41, 2.17, 4.2), median = c(43.18, 2.56, 
                                                                                                                                                                  8.87, 3.17, 1.84, 1.04, 8.14, 8.1, 4.96, 3.03, 5.87), upper = c(45.54, 
                                                                                                                                                                                                                                  3.64, 10.59, 3.75, 2.63, 1.65, 9.49, 9.45, 6.18, 4.04, 8.15)), class = "data.frame", row.names = c(NA, 
                                                                                                                                                                                                                                                                                                                                     -11L))




ui <- fluidPage(
  fluidRow(
    column(2,
           radioButtons(inputId = "site",
                        label = "Select site",
                        choiceNames = c("allsites", "site1", "site2", "site3", "site4", "site5", "site6", "site7", "site8", "site9", "site10"),
                        choiceValues = c("allsites", "site1", "site2", "site3", "site4", "site5", "site6", "site7", "site8", "site9", "site10"),
                        selected = "allsites")
    ),
    column(10,
           fluidRow( 
             tabsetPanel(type = "tabs",
                         tabPanel("Bar chart", 
                                  tags$br(),
                                  tags$br(),
                                  plotOutput("bar_chart", width = "50%")
                         )
             )
           )
    )
  )
)


server <- function(input, output, session){
  
  site <- reactive({
    input$site
  })
  
  
  output$bar_chart <- renderPlot({
    
    df1 <- df %>%
      dplyr::filter(site %in% c("allsites", site()))
    
    
    ggplot(df1, aes(site, median, fill = site))+
      geom_bar(stat = "identity", width = 0.2)+
      geom_errorbar(aes(ymin = lower, ymax = upper),
                    width = .09, size=0.2)+
      facet_wrap(~site, ncol = 2, scales ="free")
    
  })
  
}

shinyApp(ui, server)