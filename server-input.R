observe({
  if (input$data_file_type == "upload") {
    insertTab(
      inputId = "main_navbar",
      tabPanel(
        "Real-Time Analysis",
        fluidRow(id="realtime_analysis",
                 column(id="realtime_one", 6,
                        column(id="realtime_barcode", 3,
                          fluidRow(
                            div(selectInput("select", "BARCODE", choices = list("barcode1"=1,   "barcode2"=2,   "barcode3"=3,   "barcode4"=4,   "barcode5"=5,   "barcode6"=6,
                                                                                "barcode7"=7,   "barcode8"=8,   "barcode9"=9,   "barcode10"=10, "barcode11"=11, "barcode12"=12,
                                                                                "barcode13"=13, "barcode14"=14, "barcode15"=15, "barcode16"=16, "barcode17"=17, "barcode18"=18,
                                                                                "barcode19"=19, "barcode20"=20, "barcode21"=21, "barcode22"=22, "barcode23"=23, "barcode24"=24),
                                            selected = 1),
                                style="padding-left: 2%; text-align: center; width: 95%")
                          ),
                          tags$head(tags$style(HTML(".selectize-input {height: 40px}")))
                        ),
                        column(id="realtime_taxa", 3,
                               div(selectInput("select", "TAXON", choices = list("Species"=1, "Genus"=2, "Family"=3,
                                                                                 "Order"=4), selected = 1),
                                   style="padding-left: 2%; height: 40px; text-align: center")
                        ),
                        column(id="realtime_status", 3,
                               div(tags$b("STATE"), style="padding-bottom: 4px; text-align: center"),
                               verbatimTextOutput("state", placeholder = T),
                               tags$head(
                                 tags$style(
                                   HTML(
                                     "#state {
        height: 40px;
        }"
                                   )))),
                        column(id="realtime_pores", 3,
                               div(tags$b("PORES"), style="padding-bottom: 4px; text-align: center"),
                               verbatimTextOutput("pores", placeholder = T),
                               tags$head(
                                 tags$style(
                                   HTML(
                                     "#pores {
        height: 40px;
        }"))))
                        ),
                 column(id="realtime_two", 6,
                        column(id="realtime_reads", 3,
                               div(tags$b("PASSED READS"), style="padding-bottom: 4px; text-align: center"),
                               verbatimTextOutput("passed_reads", placeholder = T),
                               tags$head(
                                 tags$style(
                                   HTML(
                                     "#passed_reads {
        height: 40px;
        }")))),
                        column(id="realtime_n50", 3,
                               div(tags$b("N50 LENGTH"), style="padding-bottom: 4px; text-align: center"),
                               verbatimTextOutput("n50_length", placeholder = T),
                               tags$head(
                                 tags$style(
                                   HTML(
                                     "#n50_length {
        height: 40px;
        }")))),
                        column(id="realtime_classified", 3,
                               div(tags$b("CLASSIFIED READS"), style="padding-bottom: 4px; text-align: center"),
                               verbatimTextOutput("classified_reads", placeholder = T),
                               tags$head(
                                 tags$style(
                                   HTML(
                                     "#classified_reads {
        height: 40px;
        }")))),
                        column(id="realtime_species", 3,
                               div(tags$b("# SPECIES"), style="padding-bottom: 4px; text-align: center"),
                               verbatimTextOutput("species", placeholder = T),
                               tags$head(
                                 tags$style(
                                   HTML(
                                     "#species {
        height: 40px;
        }"))))
                        )
        ),
        fluidRow(id="realtime_readlength_analysis",
                 br(),
                 column(id="realtime_readlength_plot", 6,
                        tags$b("READ LENGTH DISTRIBUTION"),
                        downloadButton(outputId = "download_readlength_plot",
                                       label = "Download Read Length Plot (PNG)"),
                        plotOutput("plot_read_classification_plot", height = "100%"))
                 ,
                 column(id="realtime_plot", 6,
                        tags$b("Q-SCORE HISTOGRAM"),
                        downloadButton(outputId = "download_qscore_plot",
                                       label = "Download Q score Plot (PNG)"),
                        plotOutput("plot_q_score_plot", height = "100%"))),
        br(),
        fluidRow(id="realtime_read_classification",
                 column(id="realtime_readclassification_plot", 6,
                        tags$b("READS CLASSIFICATION"),
                        downloadButton(outputId = "download_classification_plot",
                                       label = "Download Reads Classification Plot (PNG)"),
                        plotOutput("plot_classification_plot",height = "100%", width = "100%")),
                 column(id="realtime_sankey_plot", 6,
                        tags$b("SANKEY PLOT"),
                        downloadButton(outputId = "download_sankey_plot",
                                       label = "Download Sankey Classification Plot (PNG)"),
                        plotOutput("plot_sankey_plot",height = "100%", width = "100%")))
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

                          abundance_data <- abundance_data_list %>% purrr::reduce(full_join, by=lineage)

                          abundance_data <- as.data.frame(abundance_data)

                          abundance_data[is.na(abundance_data)] <- 0

                          rel_abundance_matrix <- as.matrix(rel_abundance_data[,-1])

                          rownames(rel_abundance_matrix) <- rel_abundance_data[,which(colnames(rel_abundance_data)==lineage)]

                          abundance_matrix <- abundance_data[,-1]

                          rownames(abundance_matrix) <- abundance_data[,which(colnames(abundance_data)==lineage)]

                          return(list('abundance_data' = abundance_data,
                          'sample_metadata' = sample_metadata,
                          'lineage' = lineage,
                          'rel_abundance_data' = rel_abundance_data,
                          'abundance_matrix' = abundance_matrix,
                          'rel_abundance_matrix' = rel_abundance_matrix))
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

                          abundance_data <- abundance_data_list %>% purrr::reduce(full_join, by=lineage)

                          abundance_data <- as.data.frame(abundance_data)

                          abundance_data[is.na(abundance_data)] <- 0

                          rel_abundance_matrix <- as.matrix(rel_abundance_data[,-1])

                          rownames(rel_abundance_matrix) <- rel_abundance_data[,which(colnames(rel_abundance_data)==lineage)]

                          abundance_matrix <- abundance_data[,-1]

                          rownames(abundance_matrix) <- abundance_data[,which(colnames(abundance_data)==lineage)]
                          
                          return(list('sample_metadata' = sample_metadata,
                          'lineage' = lineage,
                          'rel_abundance_data' = rel_abundance_data,
                          'abundance_matrix' = abundance_matrix,
                          'rel_abundance_matrix' = rel_abundance_matrix)
                          )}
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

                          abundance_data <- abundance_data_list %>% purrr::reduce(full_join, by=lineage)

                          abundance_data <- as.data.frame(abundance_data)

                          abundance_data[is.na(abundance_data)] <- 0

                          rel_abundance_matrix <- as.matrix(rel_abundance_data[,-1])

                          rownames(rel_abundance_matrix) <- rel_abundance_data[,which(colnames(rel_abundance_data)==lineage)]

                          abundance_matrix <- abundance_data[,-1]

                          rownames(abundance_matrix) <- abundance_data[,which(colnames(abundance_data)==lineage)]

                          return(list('abundance_data' = abundance_data,
                          'sample_metadata' = sample_metadata,
                          'lineage' = lineage,
                          'rel_abundance_data' = rel_abundance_data,
                          'abundance_matrix' = abundance_matrix,
                          'rel_abundance_matrix' = rel_abundance_matrix))  
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

                          abundance_data <- abundance_data_list %>% purrr::reduce(full_join, by=lineage)

                          abundance_data <- as.data.frame(abundance_data)

                          abundance_data[is.na(abundance_data)] <- 0

                          rel_abundance_matrix <- as.matrix(rel_abundance_data[,-1])

                          rownames(rel_abundance_matrix) <- rel_abundance_data[,which(colnames(rel_abundance_data)==lineage)]

                          abundance_matrix <- abundance_data[,-1]

                          rownames(abundance_matrix) <- abundance_data[,which(colnames(abundance_data)==lineage)]

                          return(list('abundance_data' = abundance_data,
                          'sample_metadata' = sample_metadata,
                          'lineage' = lineage,
                          'rel_abundance_data' = rel_abundance_data,
                          'abundance_matrix' = abundance_matrix,
                          'rel_abundance_matrix' = rel_abundance_matrix))
                          
                          
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

                          abundance_data <- abundance_data_list %>% purrr::reduce(full_join, by=lineage)

                          abundance_data <- as.data.frame(abundance_data)

                          abundance_data[is.na(abundance_data)] <- 0

                          rel_abundance_matrix <- as.matrix(rel_abundance_data[,-1])

                          rownames(rel_abundance_matrix) <- rel_abundance_data[,which(colnames(rel_abundance_data)==lineage)]

                          return(list('abundance_data' = abundance_data,
                          'sample_metadata' = sample_metadata,
                          'lineage' = lineage,
                          'rel_abundance_data' = rel_abundance_data))

                          
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

                          rel_abundance_filtered_matrix <- rel_abundance_matrix[apply(rel_abundance_matrix, 1, function(row) any(mean(row) >0.1 )), ]

                          abundance_data <- abundance_data_list %>% purrr::reduce(full_join, by=lineage)

                          abundance_data <- as.data.frame(abundance_data)

                          abundance_data[is.na(abundance_data)] <- 0

                          abundance_matrix <- abundance_data[,-1]

                          rownames(abundance_matrix) <- abundance_data[,which(colnames(abundance_data)==lineage)]

                          rel_abundance_filtered_data <- rel_abundance_data[which(rel_abundance_data[,1] %in% rownames(rel_abundance_filtered_matrix)),]

                          return(list('abundance_data' = abundance_data,
                          'sample_metadata' = sample_metadata,
                          'lineage' = lineage,
                          'rel_abundance_data' = rel_abundance_data))
                        }
                      }

                    }
                  })
                })


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
  
  print("Data Engineering Output")
  
  analyze_data_reactive()$rel_abundance_data
  
})

output$download_results_CSV <- downloadHandler(
  filename = paste0("abundance_data_", Sys.Date(), ".csv"),
  
  content = function(file) {
    
    write.csv(analyze_data_reactive()$rel_abundance_data,
              
              file, row.names = FALSE)
  }
)