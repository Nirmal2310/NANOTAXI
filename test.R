library(shiny)
library(dplyr)
library(purrr)

server <- function(input, output, session) {
  # Reactive expression for selected taxon
  select_taxon <- reactive({
    req(input$taxa)  # Ensure input$taxa is available
    return(list('taxon' = input$taxa))
  })
  
  # Observe the selected taxon and print it
  observe({
    taxon_value <- select_taxon()$taxon
    print(paste("Selected taxon:", taxon_value))
  })
  
  # Event reactive for analyzing data
  analyze_data_reactive <- eventReactive(input$upload_data, {
    withProgress(message = "Analyzing data, please wait", {
      print("analysisCountDataReactive")
      
      if (input$data_file_type == "examplelist") {
        file_list <- gsub(".txt", "", list.files("Example/")[grep("\\_final_result.txt$", list.files("Example/"))])
        samples_header <- gsub("_final_result", "", file_list)
        
        sample_info <- read.csv("Example/Sample_Information.csv")
        rownames(sample_info) <- sample_info$Sample_Id
        
        sample_metadata <- read.csv("Example/Sample_Information.csv")
        
        sample_data_list <- list()
        abundance_data_list <- list()
        
        lineage <- select_taxon()$taxon
        
        for (i in seq_along(file_list)) {
          sample_data_list[[i]] <- read.delim(file = paste0("Example/", file_list[i], ".txt"), header = FALSE)
          colnames(sample_data_list[[i]]) <- c("READ_ID", "TAX_ID", "Superkingdom", "Phylum",
                                               "Class", "Order", "Family", "Genus", "Species")
          
          abundance_data_list[[i]] <- as.data.frame(table(sample_data_list[[i]][, which(colnames(sample_data_list[[i]]) == lineage)]))
          abundance_data_list[[i]] <- abundance_data_list[[i]] %>% mutate(Freq = (Freq / sum(Freq)) * 100)
          colnames(abundance_data_list[[i]]) <- c(lineage, samples_header[i])
        }
        
        abundance_data <- reduce(abundance_data_list, full_join, by = lineage)
        
        return(list('abundance_data' = abundance_data, 'sample_metadata' = sample_metadata))
      }
    })
  })
  
  # Observe the selected taxon to recreate abundance data
  observe({
    req(select_taxon())
    if (input$data_file_type == "examplelist") {
      file_list <- gsub(".txt", "", list.files("Example/")[grep("\\_final_result.txt$", list.files("Example/"))])
      samples_header <- gsub("_final_result", "", file_list)
      
      sample_info <- read.csv("Example/Sample_Information.csv")
      rownames(sample_info) <- sample_info$Sample_Id
      
      sample_metadata <- read.csv("Example/Sample_Information.csv")
      
      sample_data_list <- list()
      abundance_data_list <- list()
      
      lineage <- select_taxon()$taxon
      
      for (i in seq_along(file_list)) {
        sample_data_list[[i]] <- read.delim(file = paste0("Example/", file_list[i], ".txt"), header = FALSE)
        colnames(sample_data_list[[i]]) <- c("READ_ID", "TAX_ID", "Superkingdom", "Phylum",
                                             "Class", "Order", "Family", "Genus", "Species")
        
        abundance_data_list[[i]] <- as.data.frame(table(sample_data_list[[i]][, which(colnames(sample_data_list[[i]]) == lineage)]))
        abundance_data_list[[i]] <- abundance_data_list[[i]] %>% mutate(Freq = (Freq / sum(Freq)) * 100)
        colnames(abundance_data_list[[i]]) <- c(lineage, samples_header[i])
      }
      
      abundance_data <- reduce(abundance_data_list, full_join, by = lineage)
      
      output$abundance_data <- renderTable({
        abundance_data
      })
    }
  })
  
  # Output the selected taxon
  output$selected_taxon <- renderText({
    select_taxon()$taxon
  })
  
  # Output the abundance data initially
  output$abundance_data <- renderTable({
    analyze_data_reactive()$abundance_data
  })
}

shinyApp(ui = ui, server = server)