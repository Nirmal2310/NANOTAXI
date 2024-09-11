observe({
  if (input$data_file_type == "upload") {
    insertTab(
      inputId = "main_navbar",
      tabPanel
      (
        "Real-Time Analysis",
        fluidRow(id="realtime_analysis",
                 column(id="realtime_one", 6,
                        column(id="realtime_barcode", 3,
                               fluidRow(
                                div(selectInput("barcode_select", "BARCODE", choices = paste0("barcode",sprintf("%02d", 1:24)),
                                selected = "barcode01"),
                                style="padding-left: 2%; text-align: center; width: 95%; font-size: 13px; font-weight: bold;")
                              ),
                              tags$head(tags$style(HTML(".selectize-input {height: 40px}")))
                              ),
                              column(id="realtime_taxa", 3,
                              div(selectInput("taxon_select", "TAXON", choices = list("Species", "Genus", "Family",
                              "Order", "Class", "Phylum", "Kingdom"), selected = "Species"),
                              style="padding-left: 2%; height: 40px; text-align: center; font-size: 13px; font-weight: bold;")
                            ),
                            column(id="realtime_status", 3,
                            div(tags$b("STATE"), style="padding-bottom: 4px; text-align: center"),
                            verbatimTextOutput("state", placeholder = T),
                            tags$head(
                              tags$style(
                                HTML(
                                  "#state {height: 40px; font-size: 13px; text-align: center; font-weight: bold;}")
                                  )
                                  )
                            ),
                            column(id="realtime_pores", 3,
                            div(tags$b("PORES"), style="padding-bottom: 4px; text-align: center"),
                            verbatimTextOutput("pores", placeholder = T),
                            tags$head(
                              tags$style(
                                HTML("#pores {height: 40px; font-size: 13px; text-align: center; font-weight: bold;}"
                                )
                                )
                                )
                              )
                          ),
                          column(id="realtime_two", 6,
                          column(id="realtime_reads", 3,
                          div(tags$b("READS"), style="padding-bottom: 4px; text-align: center"),
                          verbatimTextOutput("passed_reads", placeholder = T),
                          tags$head(
                            tags$style(
                              HTML(
                                "#passed_reads {height: 40px; font-size: 13px; text-align: center; font-weight: bold;}"
                                )
                                )
                                )
                              ),
                              column(id="realtime_mean", 3,
                              div(tags$b("MEAN LENGTH"), style="padding-bottom: 4px; text-align: center"),
                              verbatimTextOutput("mean_length", placeholder = T),
                              tags$head(
                                tags$style(
                                  HTML(
                                    "#mean_length {height: 40px; font-size: 13px; text-align: center; font-weight: bold;}"
                                    )
                                    )
                                    )
                                ),
                                column(id="realtime_classified", 3,
                                div(tags$b("CLASSIFIED"), style="padding-bottom: 4px; text-align: center"),
                                verbatimTextOutput("classified_reads", placeholder = T),
                                tags$head(
                                  tags$style(
                                    HTML(
                                     "#classified_reads {height: 40px; font-size: 13px; text-align: center; font-weight: bold;}"
                                     )
                                     )
                                     )
                                ),
                                column(id="realtime_species", 3,
                                div(tags$b("TAXON COUNT"), style="padding-bottom: 4px; text-align: center"),
                                verbatimTextOutput("taxon_count", placeholder = T),
                                tags$head(
                                  tags$style(
                                    HTML("#taxon_count {height: 40px; font-size: 13px; text-align: center; font-weight: bold;}"
                                    )
                                    )
                                    )
                                )
                        )
        ),
        br(),
        fluidRow(id="realtime_readlength_analysis",
                 column(width = 12,
                        tabsetPanel(id="realtime_plots",
                        tabPanel(title = div("TAXON TABLE", style="font-size: 15px; font-weight: bold; color: #000000; font-family: serif"),
                        div(downloadButton(outputId = "download_taxa_table",
                        label = "Download Table (TSV)"), style = "padding-top: 10px;"),
                        DT::dataTableOutput("taxa_table", width = "50%", height = "100%")
                        ),
                        tabPanel(title = div("READ LENGTH", style="font-size: 15px; font-weight: bold; color: #000000; font-family: serif"),
                        div(downloadButton(outputId = "download_readlength_plot",
                        label = "Download Plot (PNG)"), style = "padding-top: 10px;"),
                        plotOutput("plot_read_length", height = "100%")
                        ),
                        tabPanel(title = div("Q-SCORE", style="font-size: 15px; font-weight: bold; color: #000000; font-family: serif"),
                        div(downloadButton(outputId = "download_qscore_plot",
                        label = "Download Plot (PNG)"), style = "padding-top: 10px;"),
                        plotOutput("plot_q_score_plot", height = "100%")
                        ),
                        tabPanel(title = div("CLASSIFICATION", style="font-size: 15px; font-weight: bold; color: #000000; font-family: serif"),
                        fluidRow(column(width = 2,
                        sliderInput("top_n", "Top N", min = 1, max = 25, value = 10)),
                        column(width=10,
                        div(downloadButton(outputId = "download_classification_plot",
                        label = "Download Plot (PNG)"), style = "float: right; padding-right: 15px; padding-top: 10px;"),
                        plotOutput("plot_classification_plot", height = "100%")))
                        )
                        
                 )
                        )
                        )
    ),
      target = "Cohort Analysis",
      position = "before"
    )
    }
  else {
    removeTab(inputId = "main_navbar", target = "Real-Time Analysis")
  }
})

observe(
  {
    validate(
      need((input$data_file_type == "examplelist") | (!is.null(input$inputfile))
           | (!is.null(input$metafile)), message = "Please select a File")
    )
  })

input_data_reactive <-  reactive({
  
  print("inputting data")
  
  validate(
    need((input$data_file_type == "examplelist") | (!is.null(input$inputfile)) |
           (!is.null(input$metafile)), "Please Select A File")
  )
  
  if (input$data_file_type == "examplelist") {
    
    seqdata <- read.csv("Example/Sample_Information.csv")
    
    print("Uploaded Example List")
    
    return(list('data' = seqdata))
    
  }
  else if(input$data_file_type == "upload"){
    
    if (!is.null(input$inputfile)) {
      
      file <- input$inputfile
      
      seqdata <- read.csv(file$datapath)
      
      print("Uploaded User Sample List")
      
      validate(need(ncol(seqdata) > 1,
                    message = "File appears to be one column. Check that it is a .csv file.")
      )
      return(list('data' = seqdata))
    }
  }
  else{
    
    if (!is.null(input$metafile)){
      file <- input$metafile
      
      seqdata <- read.csv(file$datapath)
      
      print("Uploaded User Sample List")
      
      validate(need(ncol(seqdata) > 1,
                    message = "File appears to be one column. Check that it is a .csv file.")
      )
      return(list('data' = seqdata))
    }
  }
})


output$fileUploaded <- reactive({
  
  return(!is.null(input_data_reactive()))
  
})

outputOptions(output, "fileUploaded", suspendWhenHidden = FALSE)


select_taxon <- reactive({
  
  return(list('taxon'=input$taxa))
  
})

observe({
  
  taxon_value <- select_taxon()$taxon
  print(taxon_value)

})

analyze_data_reactive <-
  eventReactive(c(input$upload_data,input$taxa), ignoreNULL = TRUE,
                {

                  withProgress(message = "Analyzing data, please wait",{

                    print("analysisCountDataReactive")

                    mode("Offline")

                    print(mode())

                    if (input$data_file_type == "examplelist") {
                        
                      file_list <- gsub(".txt", "", list.files("Example/")[grep("\\_final_blast_result.txt$",
                                                                                                list.files("Example/"))])
                      samples_header <- gsub("_final_blast_result","",file_list)
                                      
                      sample_info <- read.csv("Example/Sample_Information.csv")
                                      
                      rownames(sample_info) <- sample_info$Sample_Id
                                      
                      sample_metadata <- read.csv("Example/Sample_Information.csv")

                      sample_metadata$Group <- factor(sample_metadata$Group)
                                      
                      sample_data_list <- list()
                      
                      abundance_data_list <- list()

                      rel_abundance_data_list <- list()
                                      
                      lineage <- select_taxon()$taxon
                      
                      for (i in 1:length(file_list))
                        {
                        
                        sample_data_list[[i]] <- read.delim(file = paste0("Example/",file_list[i],".txt"), header = FALSE)
                        
                        colnames(sample_data_list[[i]]) <- c("TAX_ID", "Counts", "Superkingdom", "Phylum",
                                                                              "Class", "Order", "Family", "Genus", "Species")

                        sample_data_list[[i]][sample_data_list[[i]]==""] <- "Unclassified"

                        col_num <- which(colnames(sample_data_list[[i]])==lineage)

                        abundance_data_list[[i]] <- sample_data_list[[i]] %>% group_by(sample_data_list[[i]][,col_num]) %>% 
                        summarise(Counts = sum(Counts))

                        rel_abundance_data_list[[i]] <- abundance_data_list[[i]] %>% mutate(Freq = (Counts/sum(Counts))*100)

                        rel_abundance_data_list[[i]] <- rel_abundance_data_list[[i]][,-2]

                        colnames(rel_abundance_data_list[[i]]) <- c(lineage, samples_header[i])

                        colnames(abundance_data_list[[i]]) <- c(lineage, samples_header[i])
                    
                      }
                                      
                      rel_abundance_data <- rel_abundance_data_list %>% purrr::reduce(full_join, by=lineage)

                      rel_abundance_data <- as.data.frame(rel_abundance_data)

                      rel_abundance_data[is.na(rel_abundance_data)] <- 0

                      rel_abundance_matrix <- as.matrix(rel_abundance_data[,-1])

                      rownames(rel_abundance_matrix) <- rel_abundance_data[,which(colnames(rel_abundance_data)==lineage)]

                      rel_abundance_filtered_matrix <- rel_abundance_matrix[apply(rel_abundance_matrix, 1, function(row) any(mean(row) > 0.1 )), ]

                      abundance_data <- abundance_data_list %>% purrr::reduce(full_join, by=lineage)

                      abundance_data <- as.data.frame(abundance_data)

                      abundance_data[is.na(abundance_data)] <- 0

                      abundance_matrix <- abundance_data[,-1]

                      rownames(abundance_matrix) <- abundance_data[,which(colnames(abundance_data)==lineage)]

                      rel_abundance_filtered_data <- rel_abundance_data[which(rel_abundance_data[,1] %in% rownames(rel_abundance_matrix)),]

                      return(list('rel_abundance_data'=rel_abundance_data,
                      'rel_abundance_matrix'=rel_abundance_matrix,
                      'abundance_matrix'=abundance_matrix,
                      'sample_metadata'=sample_metadata,
                      'lineage'=lineage,
                      'rel_abundance_filtered_data'=rel_abundance_filtered_data,
                      'rel_abundance_filtered_matrix'=rel_abundance_filtered_matrix))
                    }
                    
                    else if(input$data_file_type == "precomputed")
                    {
                      work_dir <- getwd()
                      
                      pipeline <- input$pipeline
                      
                      if(pipeline=="BLASTn + 16s DB")
                      
                      {
                        if(input$setup)
                        {
                          print("Setting up Blastn for 16S Data Analysis")
                          
                          setwd(paste0(work_dir,"/Installation"))
                          
                          system('bash blast_install.sh')
                          
                          setwd(work_dir)
                          
                          print("Analyzing Data......")
                          
                          system(paste0('bash blastn_run.sh -p ', input$fastqdir, ' -t ',
                                        input$threads, ' -m ', input$min, ' -M ',
                                        input$max,' -i ', input$iden))
                          
                          result_dir <- paste0(input$fastqdir,"/Blast_Results/")
                          
                          
                          file_list <- gsub(".txt", "", list.files(result_dir)[grep("\\_final_blast_result.txt$",
                          list.files(result_dir))])
                          
                          samples_header <- gsub("_final_blast_result","",file_list)
                          
                          sample_info <- input_data_reactive()$data
                          
                          rownames(sample_info) <- sample_info$Sample_Id
                          
                          sample_metadata <- input_data_reactive()$data
                          
                          sample_data_list <- list()
                          
                          abundance_data_list <- list()

                          rel_abundance_data_list <- list()
                          
                          lineage <- select_taxon()$taxon
                          
                          for (i in 1:length(file_list))
                          {
                            
                            sample_data_list[[i]] <- read.delim(file = paste0(result_dir, "/", file_list[i],".txt"), header = FALSE)
                            
                            colnames(sample_data_list[[i]]) <- c("TAX_ID", "Counts", "Superkingdom", "Phylum",
                                                                "Class", "Order", "Family", "Genus", "Species")
                            
                            sample_data_list[[i]][sample_data_list[[i]]==""] <- "Unclassified"
                            
                            col_num <- which(colnames(sample_data_list[[i]])==lineage)

                            abundance_data_list[[i]] <- sample_data_list[[i]] %>% group_by(sample_data_list[[i]][,col_num]) %>% 
                            summarise(Counts = sum(Counts))

                            rel_abundance_data_list[[i]] <- abundance_data_list[[i]] %>% mutate(Freq = (Counts/sum(Counts))*100)

                            rel_abundance_data_list[[i]] <- rel_abundance_data_list[[i]][,-2]
    
                            colnames(rel_abundance_data_list[[i]]) <- c(lineage, samples_header[i])
    
                            colnames(abundance_data_list[[i]]) <- c(lineage, samples_header[i])
                            
                          }
                          
                          rel_abundance_data <- rel_abundance_data_list %>% purrr::reduce(full_join, by=lineage)

                          rel_abundance_data <- as.data.frame(rel_abundance_data)

                          rel_abundance_data[is.na(rel_abundance_data)] <- 0

                          rel_abundance_matrix <- as.matrix(rel_abundance_data[,-1])

                          rownames(rel_abundance_matrix) <- rel_abundance_data[,which(colnames(rel_abundance_data)==lineage)]

                          rel_abundance_filtered_matrix <- rel_abundance_matrix[apply(rel_abundance_matrix, 1, function(row) any(mean(row) > 0.1 )), ]

                          abundance_data <- abundance_data_list %>% purrr::reduce(full_join, by=lineage)

                          abundance_data <- as.data.frame(abundance_data)

                          abundance_data[is.na(abundance_data)] <- 0

                          abundance_matrix <- abundance_data[,-1]

                          rownames(abundance_matrix) <- abundance_data[,which(colnames(abundance_data)==lineage)]
                          
                          rel_abundance_filtered_data <- rel_abundance_data[which(rel_abundance_data[,1] %in% rownames(rel_abundance_matrix)),]

                          return(list('rel_abundance_data'=rel_abundance_data,
                          'rel_abundance_matrix'=rel_abundance_matrix,
                          'abundance_matrix'=abundance_matrix,
                          'sample_metadata'=sample_metadata,
                          'lineage'=lineage,
                          'rel_abundance_filtered_data'=rel_abundance_filtered_data,
                          'rel_abundance_filtered_matrix'=rel_abundance_filtered_matrix))
                        }
                        
                        else
                        {
                          setwd(work_dir)
                        
                          print("Analyzing Data......")
                          
                          system(paste0('bash blastn_run.sh -p ', input$fastqdir, ' -t ',
                                        input$threads, ' -m ', input$min, ' -M ',
                                        input$max,' -i ', input$iden))
                          
                          result_dir <- paste0(input$fastqdir,"/Blast_Results/")
                            
                          file_list <- gsub(".txt", "", list.files(result_dir)[grep("\\_final_blast_result.txt$",
                          list.files(result_dir))])
                          
                          samples_header <- gsub("_final_blast_result","",file_list)
                          
                          sample_info <- input_data_reactive()$data
                          
                          rownames(sample_info) <- sample_info$Sample_Id
                          
                          sample_metadata <- input_data_reactive()$data
                          
                          sample_data_list <- list()
                          
                          abundance_data_list <- list()

                          rel_abundance_data_list <- list()
                          
                          lineage <- select_taxon()$taxon
                          
                          for (i in 1:length(file_list))
                          {
                            
                            sample_data_list[[i]] <- read.delim(file = paste0(result_dir, "/", file_list[i],".txt"), header = FALSE)
                            
                            colnames(sample_data_list[[i]]) <- c("TAX_ID", "Counts", "Superkingdom", "Phylum",
                                                                "Class", "Order", "Family", "Genus", "Species")
                            
                            sample_data_list[[i]][sample_data_list[[i]]==""] <- "Unclassified"
                            
                            col_num <- which(colnames(sample_data_list[[i]])==lineage)

                            abundance_data_list[[i]] <- sample_data_list[[i]] %>% group_by(sample_data_list[[i]][,col_num]) %>% 
                            summarise(Counts = sum(Counts))

                            rel_abundance_data_list[[i]] <- abundance_data_list[[i]] %>% mutate(Freq = (Counts/sum(Counts))*100)

                            rel_abundance_data_list[[i]] <- rel_abundance_data_list[[i]][,-2]
    
                            colnames(rel_abundance_data_list[[i]]) <- c(lineage, samples_header[i])
    
                            colnames(abundance_data_list[[i]]) <- c(lineage, samples_header[i])
                            
                          }
                          
                          rel_abundance_data <- rel_abundance_data_list %>% purrr::reduce(full_join, by=lineage)

                          rel_abundance_data <- as.data.frame(rel_abundance_data)

                          rel_abundance_data[is.na(rel_abundance_data)] <- 0

                          rel_abundance_matrix <- as.matrix(rel_abundance_data[,-1])

                          rownames(rel_abundance_matrix) <- rel_abundance_data[,which(colnames(rel_abundance_data)==lineage)]

                          rel_abundance_filtered_matrix <- rel_abundance_matrix[apply(rel_abundance_matrix, 1, function(row) any(mean(row) > 0.1 )), ]

                          abundance_data <- abundance_data_list %>% purrr::reduce(full_join, by=lineage)

                          abundance_data <- as.data.frame(abundance_data)

                          abundance_data[is.na(abundance_data)] <- 0

                          abundance_matrix <- abundance_data[,-1]

                          rownames(abundance_matrix) <- abundance_data[,which(colnames(abundance_data)==lineage)]

                          rel_abundance_filtered_data <- rel_abundance_data[which(rel_abundance_data[,1] %in% rownames(rel_abundance_matrix)),]

                          return(list('rel_abundance_data'=rel_abundance_data,
                          'rel_abundance_matrix'=rel_abundance_matrix,
                          'abundance_matrix'=abundance_matrix,
                          'sample_metadata'=sample_metadata,
                          'lineage'=lineage,
                          'rel_abundance_filtered_data'=rel_abundance_filtered_data,
                          'rel_abundance_filtered_matrix'=rel_abundance_filtered_matrix))
                        }
                      }
                      
                      else if(pipeline == "Kraken2 + Greengenes")
                      {
                        if(input$setup)
                        {
                          print("Setting up Kraken2 for 16S Data Analysis")
                          
                          setwd(paste0(work_dir,"/Installation"))
                          
                          system('bash kraken_install.sh')
                          
                          setwd(work_dir)
                          
                          print("Analyzing Data......")
                          
                          system(paste0('bash kraken_run.sh -p ', input$fastqdir, ' -t ',
                                        input$threads, ' -m ', input$min, ' -M ',
                                        input$max, '-r ', input$tax))

                          result_dir <- paste0(input$fastqdir,"/Kraken2_Results/")

                          
                          file_list <- gsub(".txt", "", list.files(result_dir)[grep("\\_final_kraken2_result.txt$",
                          list.files(result_dir))])
                          
                          samples_header <- gsub("_final_kraken2_result","",file_list)
                          
                          sample_metadata <- input_data_reactive()$data
                          
                          sample_data_list <- list()
                          
                          abundance_data_list <- list()

                          rel_abundance_data_list <- list()

                          lineage <- select_taxon()$taxon
                          
                          for (i in 1:length(file_list))
                          {
                            
                            sample_data_list[[i]] <- read.delim(file = paste0(result_dir, "/", file_list[i],".txt"), header = FALSE)
                            
                            colnames(sample_data_list[[i]]) <- c("TAX_ID", "Counts", "Kingdom", "Phylum",
                                                                "Class", "Order", "Family", "Genus", "Species")
                            
                            sample_data_list[[i]][sample_data_list[[i]]==""] <- "Unclassified"

                            col_num <- which(colnames(sample_data_list[[i]])==lineage)
                            
                            abundance_data_list[[i]] <- sample_data_list[[i]] %>% group_by(sample_data_list[[i]][,col_num]) %>% 
                            summarise(Counts = sum(Counts))
                            
                            rel_abundance_data_list[[i]] <- abundance_data_list[[i]] %>% mutate(Freq = (Counts/sum(Counts))*100)
                            
                            rel_abundance_data_list[[i]] <- rel_abundance_data_list[[i]][,-2]

                            colnames(rel_abundance_data_list[[i]]) <- c(lineage, samples_header[i])

                            colnames(abundance_data_list[[i]]) <- c(lineage, samples_header[i])
                            
                          }
                          
                          rel_abundance_data <- rel_abundance_data_list %>% purrr::reduce(full_join, by=lineage)

                          rel_abundance_data <- as.data.frame(rel_abundance_data)

                          rel_abundance_data[is.na(rel_abundance_data)] <- 0

                          rel_abundance_matrix <- as.matrix(rel_abundance_data[,-1])

                          rownames(rel_abundance_matrix) <- rel_abundance_data[,which(colnames(rel_abundance_data)==lineage)]

                          rel_abundance_filtered_matrix <- rel_abundance_matrix[apply(rel_abundance_matrix, 1, function(row) any(mean(row) > 0.1 )), ]

                          abundance_data <- abundance_data_list %>% purrr::reduce(full_join, by=lineage)

                          abundance_data <- as.data.frame(abundance_data)

                          abundance_data[is.na(abundance_data)] <- 0

                          abundance_matrix <- abundance_data[,-1]

                          rownames(abundance_matrix) <- abundance_data[,which(colnames(abundance_data)==lineage)]

                          rel_abundance_filtered_data <- rel_abundance_data[which(rel_abundance_data[,1] %in% rownames(rel_abundance_matrix)),]

                          return(list('rel_abundance_data'=rel_abundance_data,
                          'rel_abundance_matrix'=rel_abundance_matrix,
                          'abundance_matrix'=abundance_matrix,
                          'sample_metadata'=sample_metadata,
                          'lineage'=lineage,
                          'rel_abundance_filtered_data'=rel_abundance_filtered_data,
                          'rel_abundance_filtered_matrix'=rel_abundance_filtered_matrix))  
                        }
                        
                        else {
                           
                          setwd(work_dir)
                          
                          print("Analyzing Data......")
                          
                          system(paste0('bash kraken_run.sh -p ', input$fastqdir, ' -t ',
                                        input$threads, ' -m ', input$min, ' -M ',
                                        input$max, '-r ', input$tax))

                          result_dir <- paste0(input$fastqdir,"/Kraken2_Results/")
                          
                          file_list <- gsub(".txt", "", list.files(result_dir)[grep("\\_final_kraken2_result.txt$",
                          list.files(result_dir))])
                          
                          samples_header <- gsub("_final_kraken2_result","",file_list)
                          
                          sample_metadata <- input_data_reactive()$data
                          
                          sample_data_list <- list()
                          
                          abundance_data_list <- list()

                          rel_abundance_data_list <- list()

                          lineage <- select_taxon()$taxon
                          
                          for (i in 1:length(file_list))
                          {
                            
                            sample_data_list[[i]] <- read.delim(file = paste0(result_dir, "/", file_list[i],".txt"), header = FALSE)
                            
                            colnames(sample_data_list[[i]]) <- c("TAX_ID", "Counts", "Kingdom", "Phylum",
                                                                "Class", "Order", "Family", "Genus", "Species")
                            
                            sample_data_list[[i]][sample_data_list[[i]]==""] <- "Unclassified"

                            col_num <- which(colnames(sample_data_list[[i]])==lineage)
                            
                            abundance_data_list[[i]] <- sample_data_list[[i]] %>% group_by(sample_data_list[[i]][,col_num]) %>% 
                            summarise(Counts = sum(Counts))
                            
                            rel_abundance_data_list[[i]] <- abundance_data_list[[i]] %>% mutate(Freq = (Counts/sum(Counts))*100)
                            
                            rel_abundance_data_list[[i]] <- rel_abundance_data_list[[i]][,-2]

                            colnames(rel_abundance_data_list[[i]]) <- c(lineage, samples_header[i])

                            colnames(abundance_data_list[[i]]) <- c(lineage, samples_header[i])
                            
                          }
                          
                          rel_abundance_data <- rel_abundance_data_list %>% purrr::reduce(full_join, by=lineage)

                          rel_abundance_data <- as.data.frame(rel_abundance_data)

                          rel_abundance_data[is.na(rel_abundance_data)] <- 0

                          rel_abundance_matrix <- as.matrix(rel_abundance_data[,-1])

                          rownames(rel_abundance_matrix) <- rel_abundance_data[,which(colnames(rel_abundance_data)==lineage)]

                          rel_abundance_filtered_matrix <- rel_abundance_matrix[apply(rel_abundance_matrix, 1, function(row) any(mean(row) > 0.1 )), ]

                          abundance_data <- abundance_data_list %>% purrr::reduce(full_join, by=lineage)

                          abundance_data <- as.data.frame(abundance_data)

                          abundance_data[is.na(abundance_data)] <- 0

                          abundance_matrix <- abundance_data[,-1]

                          rownames(abundance_matrix) <- abundance_data[,which(colnames(abundance_data)==lineage)]

                          rel_abundance_filtered_data <- rel_abundance_data[which(rel_abundance_data[,1] %in% rownames(rel_abundance_matrix)),]

                          return(list('rel_abundance_data'=rel_abundance_data,
                          'rel_abundance_matrix'=rel_abundance_matrix,
                          'abundance_matrix'=abundance_matrix,
                          'sample_metadata'=sample_metadata,
                          'lineage'=lineage,
                          'rel_abundance_filtered_data'=rel_abundance_filtered_data,
                          'rel_abundance_filtered_matrix'=rel_abundance_filtered_matrix))
                          
                          
                        }
                      }

                      else if(pipeline == "EMU + Standard DB")
                      {
                        if(input$setup)
                        {
                          print("Setting up EMU pipeline for 16S Data Analysis")
                          
                          setwd(paste0(work_dir,"/Installation"))
                          
                          system('bash emu_install.sh')
                          
                          setwd(work_dir)
                          
                          print("Analyzing Data......")
                          
                          system(paste0('bash emu_run.sh -p ', input$fastqdir, ' -t ',
                                        input$threads, ' -m ', input$min, ' -M ',
                                        input$max))

                          result_dir <- paste0(input$fastqdir,"/EMU_Results/")
                            
                          file_list <- gsub(".txt", "", list.files(result_dir)[grep("\\_final_emu_result.txt$",
                          list.files(result_dir))])
                          
                          samples_header <- gsub("_final_emu_result","",file_list)
                                                                
                          sample_info <- input_data_reactive()$data
                                                                
                          rownames(sample_info) <- sample_info$Sample_Id
                                                                
                          sample_metadata <- input_data_reactive()$data
                                                                
                          sample_data_list <- list()
                                                                
                          abundance_data_list <- list()

                          lineage <- select_taxon()$taxon

                          rel_abundance_data_list <- list()

                          for (i in 1:length(file_list))
                          {
                            sample_data_list[[i]] <- read.delim(file = paste0(result_dir, "/", file_list[i],".txt"), header = FALSE)
                            
                            colnames(sample_data_list[[i]]) <- c("TAX_ID", "Counts", "Kingdom", "Phylum",
                                                                "Class", "Order", "Family", "Genus", "Species")

                            sample_data_list[[i]][sample_data_list[[i]]==""] <- "Unclassified"

                            col_num <- which(colnames(sample_data_list[[i]]==lineage))

                            abundance_data_list[[i]] <- sample_data_list[[i]] %>% group_by(sample_data_list[[i]][,col_num]) %>% 
                            summarise(Counts = sum(Counts))

                            rel_abundance_data_list[[i]] <- abundance_data_list[[i]] %>% mutate(Freq = (Counts/sum(Counts))*100)

                            rel_abundance_data_list[[i]] <- rel_abundance_data_list[[i]][,-2]

                            colnames(rel_abundance_data_list[[i]]) <- c(lineage, samples_header[i])

                            colnames(abundance_data_list[[i]]) <- c(lineage, samples_header[i])

                          }

                          rel_abundance_data <- rel_abundance_data_list %>% purrr::reduce(full_join, by=lineage)

                          rel_abundance_data <- as.data.frame(rel_abundance_data)

                          rel_abundance_data[is.na(rel_abundance_data)] <- 0

                          rel_abundance_matrix <- as.matrix(rel_abundance_data[,-1])

                          rownames(rel_abundance_matrix) <- rel_abundance_data[,which(colnames(rel_abundance_data)==lineage)]

                          rel_abundance_filtered_matrix <- rel_abundance_matrix[apply(rel_abundance_matrix, 1, function(row) any(mean(row) > 0.1 )), ]

                          abundance_data <- abundance_data_list %>% purrr::reduce(full_join, by=lineage)

                          abundance_data <- as.data.frame(abundance_data)

                          abundance_data[is.na(abundance_data)] <- 0

                          abundance_matrix <- abundance_data[,-1]

                          rownames(abundance_matrix) <- abundance_data[,which(colnames(abundance_data)==lineage)]

                          rel_abundance_filtered_data <- rel_abundance_data[which(rel_abundance_data[,1] %in% rownames(rel_abundance_matrix)),]

                          return(list('rel_abundance_data'=rel_abundance_data,
                          'rel_abundance_matrix'=rel_abundance_matrix,
                          'abundance_matrix'=abundance_matrix,
                          'sample_metadata'=sample_metadata,
                          'lineage'=lineage,
                          'rel_abundance_filtered_data'=rel_abundance_filtered_data,
                          'rel_abundance_filtered_matrix'=rel_abundance_filtered_matrix))
                        }
                        
                        else
                        {
                          setwd(work_dir)
                          
                          print("Analyzing Data......")
                          
                          system(paste0('bash emu_run.sh -p ', input$fastqdir, ' -t ',
                                        input$threads, ' -m ', input$min, ' -M ',
                                        input$max))

                          result_dir <- paste0(input$fastqdir,"/EMU_Results/")
                          
                          file_list <- gsub(".txt", "", list.files(result_dir)[grep("\\_final_emu_result.txt$",
                          list.files(result_dir))])
                          
                          samples_header <- gsub("_final_emu_result","",file_list)
                                                                
                          sample_info <- input_data_reactive()$data
                                                                
                          rownames(sample_info) <- sample_info$Sample_Id
                                                                
                          sample_metadata <- input_data_reactive()$data
                                                                
                          sample_data_list <- list()
                                                                
                          abundance_data_list <- list()

                          lineage <- select_taxon()$taxon

                          rel_abundance_data_list <- list()

                          for (i in 1:length(file_list))
                          {
                            sample_data_list[[i]] <- read.delim(file = paste0(result_dir, "/", file_list[i],".txt"), header = FALSE)
                            
                            colnames(sample_data_list[[i]]) <- c("TAX_ID", "Counts", "Kingdom", "Phylum",
                                                                "Class", "Order", "Family", "Genus", "Species")

                            sample_data_list[[i]][sample_data_list[[i]]==""] <- "Unclassified"

                            col_num <- which(colnames(sample_data_list[[i]]==lineage))

                            abundance_data_list[[i]] <- sample_data_list[[i]] %>% group_by(sample_data_list[[i]][,col_num]) %>% 
                            summarise(Counts = sum(Counts))

                            rel_abundance_data_list[[i]] <- abundance_data_list[[i]] %>% mutate(Freq = (Counts/sum(Counts))*100)

                            rel_abundance_data_list[[i]] <- rel_abundance_data_list[[i]][,-2]

                            colnames(rel_abundance_data_list[[i]]) <- c(lineage, samples_header[i])

                            colnames(abundance_data_list[[i]]) <- c(lineage, samples_header[i])

                          }

                          rel_abundance_data <- rel_abundance_data_list %>% purrr::reduce(full_join, by=lineage)

                          rel_abundance_data <- as.data.frame(rel_abundance_data)

                          rel_abundance_data[is.na(rel_abundance_data)] <- 0

                          rel_abundance_matrix <- as.matrix(rel_abundance_data[,-1])

                          rownames(rel_abundance_matrix) <- rel_abundance_data[,which(colnames(rel_abundance_data)==lineage)]

                          rel_abundance_filtered_matrix <- rel_abundance_matrix[apply(rel_abundance_matrix, 1, function(row) any(mean(row) > 0.1 )), ]

                          abundance_data <- abundance_data_list %>% purrr::reduce(full_join, by=lineage)

                          abundance_data <- as.data.frame(abundance_data)

                          abundance_data[is.na(abundance_data)] <- 0

                          abundance_matrix <- abundance_data[,-1]

                          rownames(abundance_matrix) <- abundance_data[,which(colnames(abundance_data)==lineage)]

                          rel_abundance_filtered_data <- rel_abundance_data[which(rel_abundance_data[,1] %in% rownames(rel_abundance_matrix)),]

                          return(list('rel_abundance_data'=rel_abundance_data,
                          'rel_abundance_matrix'=rel_abundance_matrix,
                          'abundance_matrix'=abundance_matrix,
                          'sample_metadata'=sample_metadata,
                          'lineage'=lineage,
                          'rel_abundance_filtered_data'=rel_abundance_filtered_data,
                          'rel_abundance_filtered_matrix'=rel_abundance_filtered_matrix))
                        }
                      }
                    }
                  })
                })

  observeEvent(input$real_time, {
    if(input$setup)
    {
      print("Setting up EMU pipeline for 16S Data Analysis")

      work_dir <- getwd()
                          
      system(paste0('bash ',paste0(work_dir,"/Installation/kraken_ncbi_install.sh")))

      is_running(!is_running())
      updateActionButton(session, "start_analysis",
                       label = if (is_running()) "Stop Analysis" else "Start Analysis")
      if (is_running()) {
      # Reset initial delay status
        real_time_delay_done(FALSE)
      
        status_checked(FALSE)

        cohort_analysis_delay(FALSE)
      
        print(paste0("Initial Delay Started at ", Sys.time()))

        mode("Realtime")
      
        # Use shinyjs::delay to handle the initial 10-minute delay
        delay(600000, {  # 10-minute delay
          real_time_delay_done(TRUE)  # Mark initial delay as done
        })
        # Use shinyjs::delay to handle the initial 30-minute delay  
        delay(1800000, { # 30-minute delay
          cohort_analysis_delay(TRUE) # Mark initial delay as done
        })
      }
    
    }

    else {
       is_running(!is_running())
      updateActionButton(session, "start_analysis",
                       label = if (is_running()) "Stop Analysis" else "Start Analysis")
      if (is_running()) {
      # Reset initial delay status
        real_time_delay_done(FALSE)
      
        status_checked(FALSE)

        cohort_analysis_delay(FALSE)
      
        print(paste0("Initial Delay Started at ", Sys.time()))

        mode("Realtime")

        delay(600000, {
          real_time_delay_done(TRUE)
        })
        
        delay(1800000, {
          cohort_analysis_delay(TRUE)
        })
      }
    }
    
  })

  observe({
    if (is_running()) {
      
      check_sequencing_status()
      
      timer_10s()
      
      print("Sequencing Status")
      
      status_checked(TRUE)
    }
  })
  
  observe({
    if (is_running() && status_checked() && state() == "Sequencing") {
      
      run_analysis_scripts()
      
      timer_30s()
      
      print("Pore Check")
    
    }
  })
  
  observe({

    req(is_running(), real_time_delay_done(), status_checked(), state() == "Sequencing")
    
    run_realtime_scripts()

    invalidateLater(120000, session)
    
  })

  observe({

    req(is_running(), cohort_analysis_delay(), status_checked(), state() == "Sequencing")
    
    run_cohort_analysis()

    invalidateLater(600000, session)
    
  })


  
  # Function to run Python script and get sequencing status
  
  check_sequencing_status <- function() {
    
    path <- paste0(getwd(),"/Scripts")
      
    sequencing_state_check <- paste0(path,"/get_sequencing_state.py")
    
    state(py_run_file(sequencing_state_check)$status)
  }
  
  ########## Function to run python scripts to get active pores and passed reads ##########
  
  run_analysis_scripts <- function() {
    
    script_path <- paste0(getwd(),"/Scripts")
    
    pore_check <- paste0(script_path,"/get_active_pores.py")
    
    passed_reads_check <- paste0(script_path,"/get_per_barcode_reads.py")
    
    pores(py_run_file(pore_check)$pore_counts)
    
    passed_counts(py_run_file(passed_reads_check)$df)

  }

########## Function to run real-time analysis scripts ##########
  
  run_realtime_scripts <- function() {

    script_path <- paste0(getwd(), "/Scripts")

    pipeline_path <- paste0(getwd(), "/Pipelines")
    
    output_path_check <- paste0(script_path, "/get_minknow_output_dir.py")

    reads_path <- paste0(py_run_file(output_path_check)$data_directory, "/fastq_pass")
    
    system(paste0("bash ",pipeline_path,"/main_run.sh -d ",reads_path," -s ", pipeline_path, " -m 1400 -M 1800 -t Species"))

    length_list <- gsub("/.*$", "", list.files(reads_path, pattern = "average_length.txt", recursive = TRUE))

    sample_list(data.frame(Barcode = length_list))
    
    mean_read_length_df <- data.frame()

    quality_data_list <- list()

    hist_data_list <- list()

    for(i in 1:length(length_list))
    {
      
      tmp <- read.table(file = paste0(reads_path, "/", length_list[i], "/", length_list[i], "_average_length.txt"), sep = "\t", header = FALSE)

      mean_read_length_df <- rbind(mean_read_length_df, tmp)

      quality_data_list[[i]] <- read.table(file = paste0(reads_path, "/", length_list[i], "/", length_list[i], "_quality.txt"), sep = "\t", header = FALSE)

      colnames(quality_data_list[[i]]) <- "Phred"

      hist_data_list[[i]] <- read.table(file = paste0(reads_path, "/", length_list[i], "/", length_list[i], "_hist.txt"), sep = "\t", header = FALSE)

      colnames(hist_data_list[[i]]) <- c("Bin", "Counts")
    
    }
    
    colnames(mean_read_length_df) <- c("Barcode", "Mean_Length")

    mean_df(mean_read_length_df)

    hist_list(hist_data_list)

    qual_list(quality_data_list)

    classification_data_list <- list()

    classified_samples <- gsub("/.*$", "", list.files(reads_path, pattern = "kraken2_report.txt", recursive = TRUE))

    classified_samples_list(data.frame(Barcode = classified_samples))

    for (i in 1:length(classified_samples))
    {
        classification_data_list[[i]] <- read.delim(file = paste0(reads_path, "/", classified_samples[i], "/", classified_samples[i],
        "_final_kraken2_result.txt"), header = FALSE, sep = "\t")
        
        colnames(classification_data_list[[i]]) <- c("TAX_ID", "Counts", "Kingdom", "Phylum",
                                            "Class", "Order", "Family", "Genus", "Species")
        
        classification_data_list[[i]][classification_data_list[[i]]==""] <- "Unclassified"
    }

    classified_list(classification_data_list)

  }

  ############### Function to do cohort level analysis on Real-time Data after the lag of 30 minutes ###############

  run_cohort_analysis <- function()
  {
    req(classified_samples_list(), classified_list(), input$taxon_select)

    abundance_data_list <- list()

    rel_abundance_data_list <- list()

    for(i in 1:length(classified_samples_list()))
    {
      
      lineage <- input$taxon_select
      
      abundance_data_list[[i]] <- classified_list()[[i]] %>% group_by(!!sym(lineage)) %>% 
                            summarise(Counts = sum(Counts))

      rel_abundance_data_list[[i]] <- abundance_data_list[[i]] %>% mutate(Freq = (Counts/sum(Counts))*100)

      rel_abundance_data_list[[i]] <- rel_abundance_data_list[[i]][,-2]

      colnames(rel_abundance_data_list[[i]]) <- c(!!sym(lineage), classified_samples_list()$Barcode[i])

      colnames(abundance_data_list[[i]]) <- c(!!sym(lineage), classified_samples_list()$Barcode[i])
    
    }

    rel_abundance_data <- rel_abundance_data_list %>% purrr::reduce(full_join, by=!!sym(lineage))

    rel_abundance_data <- as.data.frame(rel_abundance_data)

    rel_abundance_data[is.na(rel_abundance_data)] <- 0

    rel_abundance_matrix <- as.matrix(rel_abundance_data[,-1])

    rownames(rel_abundance_matrix) <- rel_abundance_data[,which(colnames(rel_abundance_data)==!!sym(lineage))]

    rel_abundance_filtered_matrix <- rel_abundance_matrix[apply(rel_abundance_matrix, 1, function(row) any(mean(row) > 0.1 )), ]

    abundance_data <- abundance_data_list %>% purrr::reduce(full_join, by=!!sym(lineage))

    abundance_data <- as.data.frame(abundance_data)

    abundance_data[is.na(abundance_data)] <- 0

    abundance_matrix <- abundance_data[,-1]

    rownames(abundance_matrix) <- abundance_data[,which(colnames(abundance_data)==!!sym(lineage))]

    rel_abundance_filtered_data <- rel_abundance_data[which(rel_abundance_data[,1] %in% rownames(rel_abundance_matrix)),]

    real_rel_abundance_val(rel_abundance_data)

    real_abundance_matrix_val(abundance_matrix)

    real_rel_abundance_filtered_val(rel_abundance_filtered_data)

    real_rel_abundance_filtered_matrix_val(rel_abundance_filtered_matrix)
  
  }

  # Display sequencing status
  output$state <- renderText({
    state()
  })
  
  output$pores <- renderText({
    pores()
  })
  
  
  output$passed_reads <- renderText({
    
    req(input$barcode_select, passed_counts())
    
    selected_count <- passed_counts()[passed_counts()$Barcode == input$barcode_select, "Counts"]
    
    selected_count
  })

  output$mean_length <- renderText({
    
    req(input$barcode_select, mean_df(), sample_list())

    if(!input$barcode_select %in% sample_list()$Barcode)
    {
      
      mean_length <- 0
    
    }
    else {
      
      mean_length <- mean_df()[mean_df()$Barcode == input$barcode_select, "Mean_Length"]
      
      round(mean_length, digits = 0)
    
    }
  })

  output$classified_reads <- renderText({
    
    req(input$barcode_select, classified_list(), input$taxon_select, classified_samples_list())
    
    if(!input$barcode_select %in% classified_samples_list()$Barcode) {
      
      classified_reads <- 0
    
    }
    
    else {
      
      df <- classified_list()[[which(classified_samples_list()$Barcode == input$barcode_select)]]

      df <- df %>% dplyr::select(c(!!sym(input$taxon_select), Counts)) %>% 
      dplyr::filter(!!sym(input$taxon_select) !="Unclassified") %>% dplyr::summarise(Counts = sum(Counts))
    
      classified_reads <- df$Counts

      classified_reads
    }
  
  })

  output$taxon_count <- renderText({
    
    req(input$barcode_select, classified_list(), input$taxon_select, classified_samples_list())
    
    if(!input$barcode_select %in% classified_samples_list()$Barcode) {

      unique_taxon <- 0

    }

    else {
      
      df <- classified_list()[[which(classified_samples_list()$Barcode == input$barcode_select)]]

      taxa <- input$taxon_select
      
      unique_taxon <- df %>% dplyr::select(c(!!sym(taxa), Counts)) %>% 
      dplyr::filter(!!sym(taxa) != "Unclassified") %>% dplyr::group_by(!!sym(taxa)) %>% dplyr::summarise(Counts = sum(Counts)) %>% nrow()

      unique_taxon

    }
  })

  output$plot_classification_plot <- renderPlot({

    req(input$barcode_select, classified_list(), input$taxon_select, classified_samples_list(), input$top_n)

    if(!input$barcode_select %in% classified_samples_list()$Barcode) {

      return(NULL)

    }

    else {
      
      df <- classified_list()[[which(classified_samples_list()$Barcode == input$barcode_select)]]

      taxa <- input$taxon_select

      num <- input$top_n

      df <- df %>% dplyr::select(c(!!sym(taxa), Counts)) %>%
      filter(!!sym(taxa) != "Unclassified") %>% 
      group_by(!!sym(taxa)) %>% summarise(Counts = sum(Counts)) %>%
      mutate(Abundance = (Counts/sum(Counts))*100) %>% arrange(desc(Counts)) %>% head(n=num)

      classification_plot <- ggplot(df, aes(!!(sym(taxa)), Counts)) +
      geom_bar(stat = "identity", fill = "#00789A", color = "black") +
      theme_minimal() +
      ylab("Read Counts") +
      theme(
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 15, face = "bold", colour = "black"),
        axis.line = element_line(colour = "black", linewidth = 0.5, linetype = "solid" ),
        strip.text.x = element_text(size = 15, face = "bold", colour = "black"),
        axis.text.y= element_text(size=15, face = "bold", colour = "black"),
        axis.text.x= element_text(size=15, face = "bold", colour = "black", angle = 45, vjust = 1, hjust=1),
        legend.position = "none",
        title = element_text(size = 15, face = "bold")
      )

      plot_taxa_bar(classification_plot)

      classification_plot

    }

  }, height = 500)

  output$plot_read_length <- renderPlot({

    req(input$barcode_select, sample_list(), hist_list())

    if(!input$barcode_select %in% sample_list()$Barcode) {

      return(NULL)

    }

    else {
      
      df <- hist_list()[[which(sample_list()$Barcode == input$barcode_select)]]

      df <- df %>% group_by(Bin) %>% summarise(Counts = sum(Counts)) %>% filter(Bin>=1000 & Bin<=2000) %>% 
      mutate(Selected = ifelse(Bin>=1400 & Bin<=1800, "Yes", "No"))

      cols <- c("Yes" = "#631879FF", "No" = "#008280FF")

      histogram_plot <- ggplot(df, aes(Bin, Counts, fill = Selected)) +
      geom_bar(stat = "identity", color = "black") +
      geom_text(aes(label=Counts), fontface = "bold", position=position_dodge(width=0.9), vjust=-0.25, color = "black", size = 5) +
      scale_fill_manual(values = cols) +
      theme_minimal() +
      theme(
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        axis.text.y=element_text(size=14, face = "bold"),
        axis.text.x=element_text(size=14, face = "bold"),
        axis.line = element_line(colour = "black", linewidth = 0.5, linetype = "solid" ),
        legend.title=element_blank(),
        legend.text=element_blank(),
        plot.title = element_text(hjust = 0.5, size = rel(2)),
        legend.position = "none"
      ) +
      xlab("Read Length") +
      ylab("Read Counts")
      
      plot_read_histogram(histogram_plot)
      
      histogram_plot
      
    }
  
  }, height = 500)

  output$plot_q_score_plot <- renderPlot({

    req(input$barcode_select, sample_list(), qual_list())

    if(!input$barcode_select %in% sample_list()$Barcode) {

      return(NULL)

    }

    else {
       
      df <- qual_list()[[which(sample_list()$Barcode == input$barcode_select)]]

      df <- df %>%
      mutate(Bin = cut(Phred, breaks = c(seq(1, max(Phred), 1), max(Phred)), include.lowest = TRUE, right = TRUE, labels = as.character(c(seq(1, max(Phred), 1))))) %>%
      select(-Phred) %>%
      group_by(Bin) %>% summarise(Counts = n())

      colnames(df) <- c("Q_Score", "Frequency")

      df$Q_Score <- as.numeric(df$Q_Score)

      df <- df %>% mutate(Selected = ifelse(Q_Score>=10, "Yes", "No"))

      cols <- c("Yes" = "#631879FF", "No" = "#008280FF")

      q_score_plot <- df %>% ggplot(aes(Q_Score, Frequency, fill = Selected)) +
        geom_bar(stat = "identity", color = "black") +
        theme_minimal() +
        xlab("Phred Score") +
        ylab("Read Counts") +
        theme(
          axis.title.x = element_text(size = 14, face = "bold", colour = "black"),
          axis.title.y = element_text(size = 14, face = "bold", colour = "black"),
          strip.text.x = element_text(size = 14, face = "bold", colour = "black"),
          axis.text.y= element_text(size=14, face = "bold", colour = "black"),
          axis.text.x= element_text(size=14, face = "bold", colour = "black"),
          axis.line = element_line(colour = "black", linewidth = 0.5, linetype = "solid" ),
          strip.background = element_blank(),
          legend.position = "none"
        ) +
        scale_fill_manual(values = cols)
      
      plot_q_histogram(q_score_plot)

      q_score_plot

    }

  }, height = 500)

  output$taxa_table <- renderDataTable({

    req(input$barcode_select, classified_list(), input$taxon_select, classified_samples_list())

    if(!input$barcode_select %in% classified_samples_list()$Barcode) {

      return(NULL)

    }

    else {
      
      df <- classified_list()[[which(classified_samples_list()$Barcode == input$barcode_select)]]

      taxa <- input$taxon_select

      df <- df  %>% dplyr::select(c(!!sym(taxa), Counts)) %>% 
      filter(!!sym(taxa) != "Unclassified") %>% 
      group_by(!!sym(taxa)) %>% summarise(Counts = sum(Counts)) %>% 
      mutate(Abundance = (Counts/sum(Counts))*100)

      df <- df %>% arrange(desc(Counts))

      df <- df %>% select(!!sym(taxa), Counts)

      taxa_count_table(df)

      formattable(df, list(area(col = Counts) ~ normalize_bar("pink", 0.2))) %>% 
      as.datatable(escape = FALSE, options = list(scrollX = TRUE), rownames = FALSE)
    }

  }, height = 500)

  output$download_readlength_plot <- downloadHandler(
    req(input$barcode_select),
    filename = function() {
      paste0("Read_Length_Histogram_", input$barcode_select, "_", Sys.Date(), ".png")
    },
    content = function(file) {
      ggsave(file, plot_read_histogram(), width = 13.69, height = 8.27,
      units = "in", dpi = "retina", bg = "white")
    }
  )

  output$download_qscore_plot <- downloadHandler(
    req(input$barcode_select),
    filename = function() {
      paste0("Phred_Score_Histogram_", input$barcode_select, "_", Sys.Date(), ".png")
    },
    content = function(file) {
      ggsave(file, plot_q_histogram(), width = 13.69, height = 8.27,
      units = "in", dpi = "retina", bg = "white")
    }
  )

  output$download_classification_plot <- downloadHandler(
    req(input$barcode_select, input$taxon_select),
    filename = function() {
      paste0(input$taxon_select,"_Classification_", input$barcode_select, "_", Sys.Date(), ".png")
    },
    content = function(file) {
      ggsave(file, plot_taxa_bar(), width = 13.69, height = 8.27,
      units = "in", dpi = "retina", bg = "white")
    }
  )

  output$download_taxa_table <- downloadHandler(
    req(input$barcode_select, input$taxon_select),
    filename = function() {
      paste0(input$taxon_select, "_Counts_", input$barcode_select, "_", Sys.Date(), ".csv")
    },
    content = function(file) {
      write.table(taxa_count_table(), file, row.names = FALSE, sep = "\t", quote = FALSE, col.names = TRUE)
    }
  )


output$sampleinfo <- renderDataTable({
  
  temp <- input_data_reactive()
  
  if(!is.null(temp)) temp$data
  
})

observeEvent(input$upload_data,
             ({
               
               updateCollapse(session,id = "input_collapse_panel", open = "analysis_panel",
                              style = list("analysis_panel" = "success", "data_panel" = "primary"))  
             }))

output$analysisoutput <- renderDataTable({

  if(mode() == "Offline") {
    print("Data Engineering Output")
  
    analyze_data_reactive()$rel_abundance_data
  }

  else if(mode() == "Realtime") {
    
    print("Real-time Analysis Output")

    req(real_rel_abundance_val())

    real_rel_abundance_val()
  }

})

output$download_results_CSV <- downloadHandler(
  filename = paste0("abundance_data_", Sys.Date(), ".tsv"),
  
  content = function(file) {

    if(mode() == "Offline")
    {
      write.tsv(analyze_data_reactive()$rel_abundance_data,
              
              file, row.names = FALSE, quote = FALSE)
    }
    else if(mode() == "Realtime") {
        
        req(real_rel_abundance_val())

        write.tsv(real_rel_abundance_val(),
                
                file, row.names = FALSE, quote = FALSE)
    }
  }
)