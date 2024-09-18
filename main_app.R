options(shiny.maxRequestSize = 100*1024^2)

source("packages.R")

# UI
ui <- navbarPage(title = div(class="titleimg",img(src="Nanotaxi.png", height=100, width=100), "",
              style="position: absolute; top: 3px; left: 40px"),
              id = "main_navbar",
  tags$head(
    tags$style(HTML('.navbar-nav {
                            height: 110px;
                            padding-left: 25px;
                            padding-top: 30px
                            }'),
               HTML('.navbar-brand {width: 110px; font-size:35px; text-align:center}'))
  ),
  useShinyjs(),
  source("ui-tab-intro.R", local = TRUE)$value,
  source("ui-tab-input.R", local = TRUE)$value,
  source("ui-tab-cohort-analysis.R", local = TRUE)$value,
  source("ui-tab-help.R", local = TRUE)$value,
  source("ui-tab-conditions.R", local = TRUE)$value
)

# Set up Python environment

conda_path <- system('if [ $(which conda | grep "condabin") ]; then conda_path=$(which conda | sed "s/\\/condabin.*$//g"); else conda_path=$(which conda | sed "s/\\/bin.*$//g"); fi && echo $conda_path', intern = TRUE)

Sys.setenv(RETICULATE_PYTHON = paste0(conda_path, "/envs/minknow_api/bin/python"))

# Server

server <- function(input, output, session) {
  observe({
    if (input$data_file_type == "upload") {
      insertTab(
        inputId = "main_navbar",
        source("ui-tab-real-time-analysis.R", local = TRUE)$value,
        target = "Cohort Analysis",
        position = "before")}
        else {
          removeTab(inputId = "main_navbar", target = "Real-Time Analysis")}
          })
  
  state <- reactiveVal()

  route <- reactiveVal()
  
  is_running <- reactiveVal(FALSE)

  cohort_run <- reactiveVal(FALSE)
  
  timer_10s <- reactiveTimer(10000)

  initial_delay_done <- reactiveVal(FALSE)

  cohort_delay_done <- reactiveVal(FALSE)

  status_checked <- reactiveVal(FALSE)
  
  pores <- reactiveVal()
  
  passed_counts <- reactiveVal()

  sample_list <- reactiveVal()

  classified_samples_list <- reactiveVal()

  classified_list <- reactiveVal()

  cohort_analysis_list <- reactiveVal()

  hist_list <- reactiveVal()

  qual_list <- reactiveVal()

  mean_df <- reactiveVal()

  taxa_count_table <- reactiveVal()

  plot_read_histogram <- reactiveVal()

  plot_q_histogram <- reactiveVal()

  plot_taxa_bar <- reactiveVal()

  result_dir_val <- reactiveVal()

  rel_abundance_filtered_val <- reactiveVal()
    
  plot_taxa_stacked <- reactiveVal()

  plot_diversity_box <- reactiveVal()

  plot_pcoa_dot <- reactiveVal()

  plot_nmds_dot <- reactiveVal()

  plot_pca_dot <- reactiveVal()

  table_permanova <- reactiveVal()

  plot_real_heatmap <- reactiveVal()

  observeEvent(input$example_run, {
    
    route("Example")
  
  })

  observeEvent(input$upload_data, {
    
    cohort_run(!cohort_run())
    
    route("Offline")

    updateActionButton(session, "upload_data",
                       label = if(cohort_run()) "Stop Analysis" else "Start Analysis")
  
  })
  
  observeEvent(input$start_analysis, {
    
    is_running(!is_running())
    
    route("Realtime")
    
    updateActionButton(session, "start_analysis", 
                       label = if (is_running()) "Stop Analysis" else "Start Analysis")
    
    if (is_running()) {
      # Reset initial delay status
      initial_delay_done(FALSE)
      
      status_checked(FALSE)
      
      print(paste0("Initial Delay Started at ", Sys.time()))
      
      delay(600000, {
        initial_delay_done(TRUE)
      })

      print(paste0("Cohort Delay Started at ", Sys.time()))

      delay(1800000, {
        cohort_delay_done(TRUE)
      })
    }
  })

  observe({
    if (is_running()) {
      
      check_sequencing_status()
      
      timer_10s()
      
      status_checked(TRUE)
    }
  })
  
  observe({
    if (is_running() && status_checked() && state() == "Sequencing") {
      
      run_analysis_scripts()
      
      timer_10s()
    
    }
  })
  
  observe({

    req(is_running(), initial_delay_done(), status_checked(), state() == "Sequencing")
    
    run_realtime_scripts()

    invalidateLater(120000, session)
    
  })

  observe({
    
    req(is_running(), cohort_delay_done(), status_checked(), state()=="Sequencing", classified_list())

    cohort_analysis_list(classified_list())

    invalidateLater(600000, session)
  
  })

  observe({
    
    req(cohort_run(), route()=="Offline")

    work_dir <- getwd()

    install_dir <- paste0(work_dir, "/Installation")

    pipeline_dir <- paste0(work_dir, "/Pipelines")

    if(input$pipeline=="BLASTn + 16s DB")
    {
      if(input$setup)
      {
        print("Setting up Blastn for 16S Data Analysis....")

        system(paste0('bash ', install_dir,'/blast_install.sh'))

        print("Analyzing Data....")

        system(paste0("bash ", pipeline_dir, "/blast_run.sh -p ", input$fastqdir,
        " -t ", input$threads, " -m ", input$min, " -M ", input$max, " -i ", input$iden))

        result_dir <- paste0(input$fastqdir, "/Blast_Results/")

        result_dir_val(result_dir)

      }
      else
      {
        print("Analyzing Data....")

        system(paste0("bash ", pipeline_dir, "/blast_run.sh -p ", input$fastqdir,
        " -t ", input$threads, " -m ", input$min, " -M ", input$max, " -i ", input$iden))

        result_dir <- paste0(input$fastqdir, "/Blast_Results/")

        result_dir_val(result_dir)
      
      }
    }
    
    else if(input$pipeline == "Kraken2 + Greengenes")
    {
      if(input$setup)
      {
        print("Setting up Kraken2 for 16S Data Analysis")

        system(paste0("bash ", install_dir, "/kraken_install.sh"))

        print("Analyzing Data....")

        system(paste0("bash ", pipeline_dir, "/kraken_run.sh -p ", input$fastqdir, ' -t ',
        input$threads, ' -m ', input$min, ' -M ',
        input$max, ' -r ', input$tax))

        result_dir <- paste0(input$fastqdir, "/Kraken2_Results/")

        result_dir_val(result_dir)

      }
      else
      {
        print("Analyzing Data....")

        system(paste0("bash ", pipeline_dir, "/kraken_run.sh -p ", input$fastqdir, ' -t ',
        input$threads, ' -m ', input$min, ' -M ',
        input$max, ' -r ', input$tax))

        result_dir <- paste0(input$fastqdir, "/Kraken2_Results/")

        result_dir_val(result_dir)
      
      }
    }

    else if(input$pipeline == "EMU + Standard DB")
    {
      if(input$setup)
      {
        print("Setting up EMU Pipeline for 16S Data Analysis....")

        system(paste0("bash ", install_dir, "/emu_install.sh"))

        print("Analyzing Data....")

        system(paste0("bash ", pipeline_dir, "/emu_run.sh -p ", input$fastqdir, ' -t ',
        input$threads, ' -m ', input$min, ' -M ',
        input$max))

        result_dir <- paste0(input$fastqdir,"/EMU_Results/")  

        result_dir_val(result_dir)

      }
      else
      {
        print("Analyzing Data....")

        system(paste0("bash ", pipeline_dir, "/emu_run.sh -p ", input$fastqdir, ' -t ',
        input$threads, ' -m ', input$min, ' -M ',
        input$max))
      
        result_dir <- paste0(input$fastqdir,"/EMU_Results/")

        result_dir_val(result_dir)

      }
    }
  
  })

  check_sequencing_status <- function() {
    
    path <- paste0(getwd(),"/Scripts")
      
    sequencing_state_check <- paste0(path,"/get_sequencing_state.py")
    
    state(py_run_file(sequencing_state_check)$status)
  }
  
  run_analysis_scripts <- function() {
    
    script_path <- paste0(getwd(),"/Scripts")
    
    pore_check <- paste0(script_path,"/get_active_pores.py")
    
    passed_reads_check <- paste0(script_path,"/get_per_barcode_reads.py")
    
    pores(py_run_file(pore_check)$pore_counts)
    
    passed_counts(py_run_file(passed_reads_check)$df)

  }

  run_realtime_scripts <- function() {

    script_path <- paste0(getwd(), "/Scripts")

    pipeline_path <- paste0(getwd(), "/Pipelines")
    
    output_path_check <- paste0(script_path, "/get_minknow_output_dir.py")

    reads_path <- paste0(py_run_file(output_path_check)$data_directory, "/fastq_pass")
    
    system(paste0("bash ",pipeline_path,"/main_run.sh -d ",reads_path," -s ", pipeline_path, " -m 1400 -M 1800 -t Species 2>&1 | tee -a ", reads_path, "/realtime_run.log"))

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

  example_analysis <- function(path, file_list, lineage) 
  {
    samples_header <- gsub("_final_blast_result","",file_list)
    
    sample_data_list <- list()
    
    abundance_data_list <- list()
    
    rel_abundance_data_list <- list()
    
    for(i in 1:length(file_list))
    {
      sample_data_list[[i]] <- read.delim(paste0(path,"/",file_list[i],".txt"), header = FALSE)

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

    rel_abundance_filtered_data <- rel_abundance_data[which(rel_abundance_data[,1] %in% rownames(rel_abundance_filtered_matrix)),] %>% as.data.frame()

    abundance_data <- abundance_data_list %>% purrr::reduce(full_join, by=lineage)

    abundance_data <- as.data.frame(abundance_data)

    abundance_data[is.na(abundance_data)] <- 0

    abundance_matrix <- abundance_data[,-1]

    rownames(abundance_matrix) <- abundance_data[,which(colnames(abundance_data)==lineage)]

    return(list('rel_abundance_data' = rel_abundance_data,
    'rel_abundance_filtered_data' = rel_abundance_filtered_data,
    'rel_abundance_matrix' = rel_abundance_matrix,
    'rel_abundance_filtered_matrix' = rel_abundance_filtered_matrix,
    'abundance_matrix' = abundance_matrix))
  }
  
  cohort_realtime_analysis <- function(classification_list, sample_list, lineage)
  {
    data <- classification_list

    abundance_data_list <- list()

    rel_abundance_data_list <- list()

    for(i in 1:length(sample_list))
      {

        col_num <- which(colnames(data[[i]])==lineage)
        
        abundance_data_list[[i]] <- data[[i]] %>% group_by(data[[i]][, col_num]) %>% 
        summarise(Counts = sum(Counts))

        rel_abundance_data_list[[i]] <- abundance_data_list[[i]] %>% mutate(Freq = (Counts/sum(Counts))*100)

        rel_abundance_data_list[[i]] <- rel_abundance_data_list[[i]][,-2]

        colnames(rel_abundance_data_list[[i]]) <- c(lineage, sample_list[i])

        colnames(abundance_data_list[[i]]) <- c(lineage, sample_list[i])
      
      }

    rel_abundance_data <- rel_abundance_data_list %>% purrr::reduce(full_join, by=lineage)

    rel_abundance_data <- as.data.frame(rel_abundance_data)

    rel_abundance_data[is.na(rel_abundance_data)] <- 0

    rel_abundance_matrix <- as.matrix(rel_abundance_data[,-1])

    rownames(rel_abundance_matrix) <- rel_abundance_data[,which(colnames(rel_abundance_data)==lineage)]

    rel_abundance_filtered_matrix <- rel_abundance_matrix[apply(rel_abundance_matrix, 1, function(row) any(mean(row) > 0.1 )), ]

    rel_abundance_filtered_data <- rel_abundance_data[which(rel_abundance_data[,1] %in% rownames(rel_abundance_filtered_matrix)),] %>% as.data.frame()

    abundance_data <- abundance_data_list %>% purrr::reduce(full_join, by=lineage)

    abundance_data <- as.data.frame(abundance_data)

    abundance_data[is.na(abundance_data)] <- 0

    abundance_matrix <- abundance_data[,-1]

    rownames(abundance_matrix) <- abundance_data[,which(colnames(abundance_data)==lineage)]

    return(list('rel_abundance_data' = rel_abundance_data,
    'rel_abundance_filtered_data' = rel_abundance_filtered_data,
    'rel_abundance_matrix' = rel_abundance_matrix,
    'rel_abundance_filtered_matrix' = rel_abundance_filtered_matrix,
    'abundance_matrix' = abundance_matrix))
  }

  cohort_offline_analysis <- function(result_dir, lineage) {
    
    file_list <- gsub(".txt", "", list.files(result_dir))[grep("\\_final.*result.txt", list.files(result_dir))]

    samples_header <- gsub("_final.*result", "", file_list)

    sample_data_list <- list()

    abundance_data_list <- list()

    rel_abundance_data_list <- list()

    for (i in 1:length(file_list))
    {
      sample_data_list[[i]] <- read.delim(file = paste0(result_dir, "/", file_list[i], ".txt"), header = FALSE)

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

    rel_abundance_filtered_data <- rel_abundance_data[which(rel_abundance_data[,1] %in% rownames(rel_abundance_matrix)),] %>% as.data.frame()

    abundance_data <- abundance_data_list %>% purrr::reduce(full_join, by=lineage)

    abundance_data <- as.data.frame(abundance_data)

    abundance_data[is.na(abundance_data)] <- 0

    abundance_matrix <- abundance_data[,-1]

    rownames(abundance_matrix) <- abundance_data[,which(colnames(abundance_data)==lineage)]

    return(list('rel_abundance_data'=rel_abundance_data,
    'rel_abundance_matrix'=rel_abundance_matrix,
    'rel_abundance_filtered_data'=rel_abundance_filtered_data,
    'rel_abundance_filtered_matrix'=rel_abundance_filtered_matrix,
    'abundance_matrix' = abundance_matrix))
  }
  
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

  output$sampleinfo <- renderDataTable({
    
    temp <- input_data_reactive()
    
    if(!is.null(temp)) temp$data
    
  })

  output$analysisoutput <- DT::renderDataTable({
  
    print("Data Engineering Output")

    if(route()=="Realtime")
    {
      req(cohort_analysis_list(), input$taxa, classified_samples_list(), input_data_reactive())
      
      sample_list <- classified_samples_list()$Barcode

      lineage <- input$taxa

      rel_abundance_filtered_data <- cohort_realtime_analysis(cohort_analysis_list(), sample_list, lineage)$rel_abundance_filtered_data

      rel_abundance_filtered_val(rel_abundance_filtered_data)

      rel_abundance_filtered_data
    
    }
    
    else if(route()=="Example")
    {
      req(input$taxa)
        
      file_list <- gsub(".txt", "", list.files("Example/")[grep("\\_final_blast_result.txt$",
                                                                                                list.files("Example/"))])
                      
      lineage <- input$taxa

      rel_abundance_filtered_data <- example_analysis("Example", file_list, lineage)$rel_abundance_filtered_data

      rel_abundance_filtered_val(rel_abundance_filtered_data)

      rel_abundance_filtered_data
    
    }

    else if(route()=="Offline")
    {
      req(result_dir_val(), input$taxa)

      lineage <- input$taxa

      dir <- result_dir_val()

      rel_abundance_filtered_data <- cohort_offline_analysis(dir, lineage)$rel_abundance_filtered_data

      rel_abundance_filtered_val(rel_abundance_filtered_data)

      rel_abundance_filtered_data
    
    }
  
  })
    
  observeEvent(input$example_run,
  ({
    updateCollapse(session,id = "input_collapse_panel", open = "analysis_panel",
    style = list("analysis_panel" = "success", "data_panel" = "primary"))  
  }))

  observeEvent(input$upload_data,
  ({
    updateCollapse(session,id = "input_collapse_panel", open = "analysis_panel",
    style = list("analysis_panel" = "success", "data_panel" = "primary"))  
  }))

  output$download_results_csv <- downloadHandler(
      filename = function() {
        req(input$taxa)
        lineage <- input$taxa
        paste0(lineage,"_Relative_Abundance_Data_", Sys.Date(), ".csv")
        },
      content = function(file) {
        write.csv(rel_abundance_filtered_val(), file, row.names = FALSE, quote = FALSE)
      }
  )

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
      geom_bar(stat = "identity", fill = "#00789A", color = "black", width = 0.35) +
      theme_minimal() +
      ylab("Read Counts") +
      theme(
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 15, face = "bold", colour = "black"),
        axis.line = element_line(colour = "black", linewidth = 0.5, linetype = "solid" ),
        strip.text.x = element_text(size = 15, face = "bold", colour = "black"),
        axis.text.y= element_text(size=15, face = "bold", colour = "black"),
        axis.text.x= element_text(size=15, face = "bold", colour = "black", angle = 90, vjust = 1, hjust=1),
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
  
  }, height = 350)

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
        geom_bar(stat = "identity", color = "black", width = 0.35) +
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

  }, height = 350)

  output$taxa_table <- DT::renderDataTable({

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

      df <- df %>% dplyr::select(!!sym(taxa), Counts)

      taxa_count_table(df)

      formattable(df, list()) %>% 
      as.datatable(escape = FALSE, options = list(scrollX = TRUE), rownames = FALSE)
    }

  }, height = 350)

  output$plot_stacked_barplot <- renderPlotly({

    if(route()=="Realtime")
    {
      req(cohort_analysis_list(), input$taxa, input$top_taxa, classified_samples_list())

      sample_list <- classified_samples_list()$Barcode

      top_n <- input$top_taxa

      lineage <- input$taxa
      
      stacked_df <- cohort_realtime_analysis(cohort_analysis_list(), sample_list, lineage)$rel_abundance_data

      required_col <- which(colnames(stacked_df)==lineage)

      stacked_df[,required_col] <- as.character(stacked_df[,required_col])

      stacked_df <- stacked_df %>% 
        pivot_longer(cols = -all_of(required_col), names_to = "Sample_Id", values_to = "Abundance") %>% 
        filter(Abundance > 0)

      stacked_df <- stacked_df %>% group_by(Sample_Id) %>%
        arrange(desc(Abundance)) %>% dplyr::slice(1:ifelse(n()<top_n,n(),top_n))

      stacked_df$Sample_Id <- gsub("barcode","", stacked_df$Sample_Id)

      temp <- stacked_df %>% group_by(Sample_Id) %>% summarise(Abundance = 100-sum(Abundance)) %>% mutate(Species="Others") %>% 
        dplyr::select(Species, Sample_Id, Abundance)

      colnames(temp) <- colnames(stacked_df)

      stacked_df <- rbind(stacked_df, temp) %>% arrange(Sample_Id)

      total_abundance <- stacked_df %>%
        group_by(!!sym(lineage)) %>%
        summarise(Total_Abundance = sum(Abundance)) %>%
        arrange(Total_Abundance)

      stacked_df[[required_col]] <- factor(stacked_df[[required_col]], 
                                  levels = c("Others", total_abundance[[required_col]][total_abundance[[required_col]] != "Others"]))

      stacked_df <- stacked_df %>% filter(!!sym(lineage)!="Unclassified")

      fill_colors <- c("Others" = "#D3D3D3",setNames(viridis_pal(option = "D")(
        length(levels(stacked_df[[required_col]]))-1), 
        levels(stacked_df[[required_col]])[levels(stacked_df[[required_col]]) != "Others"]))

      stacked_barplot <- ggplot(stacked_df, aes(x = Sample_Id, y = Abundance, fill = !!sym(lineage))) +
      geom_bar(stat = "identity", position = "stack", color = "black", linewidth=0.2) +
      theme_linedraw() +
      theme(
        axis.text.x = element_text(size = 10, face = "bold"),
        axis.text.y = element_text(size = 10, face = "bold"),
        axis.title.x = element_text(size = 10, face = "bold"),
        axis.title.y = element_text(size = 10, face = "bold"),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(face = "bold"),
        legend.position = "right",
        title = element_text(size = 10, face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      ) +
      guides(fill = guide_legend(ncol = 4)) +
      labs(x = "BARCODE", y="% RELATIVE ABUNDANCE", fill = lineage) + 
      scale_fill_manual(values = fill_colors) +
      scale_y_continuous(expand = c(0, 0), limits = c(0,102))

      plot_taxa_stacked(stacked_barplot)

      ggplotly(stacked_barplot)
    
    }

    else if(route()=="Example")
    {
      req(input$taxa, input$top_taxa)

      file_list <- gsub(".txt", "", list.files("Example/")[grep("\\_final_blast_result.txt$",
                                                                                                list.files("Example/"))])
      lineage <- input$taxa

      top_n <- input$top_taxa
      
      stacked_df <- example_analysis("Example", file_list, lineage)$rel_abundance_data

      required_col <- which(colnames(stacked_df)==lineage)

      stacked_df[,required_col] <- as.character(stacked_df[,required_col])

      stacked_df <- stacked_df %>% 
        pivot_longer(cols = -all_of(required_col), names_to = "Sample_Id", values_to = "Abundance") %>% 
        filter(Abundance > 0)

      stacked_df <- stacked_df %>% group_by(Sample_Id) %>%
        arrange(desc(Abundance)) %>% dplyr::slice(1:ifelse(n()<top_n,n(),top_n))

      stacked_df$Sample_Id <- gsub("barcode","", stacked_df$Sample_Id)

      temp <- stacked_df %>% group_by(Sample_Id) %>% summarise(Abundance = 100-sum(Abundance)) %>% mutate(Species="Others") %>% 
        dplyr::select(Species, Sample_Id, Abundance)

      colnames(temp) <- colnames(stacked_df)

      stacked_df <- rbind(stacked_df, temp) %>% arrange(Sample_Id)

      total_abundance <- stacked_df %>%
        group_by(!!sym(lineage)) %>%
        summarise(Total_Abundance = sum(Abundance)) %>%
        arrange(Total_Abundance)

      stacked_df[[required_col]] <- factor(stacked_df[[required_col]], 
                                  levels = c("Others", total_abundance[[required_col]][total_abundance[[required_col]] != "Others"]))

      stacked_df <- stacked_df %>% filter(!!sym(lineage)!="Unclassified")
      
      fill_colors <- c("Others" = "#D3D3D3",setNames(viridis_pal(option = "D")(
        length(levels(stacked_df[[required_col]]))-1), 
        levels(stacked_df[[required_col]])[levels(stacked_df[[required_col]]) != "Others"]))

      stacked_barplot <- ggplot(stacked_df, aes(x = Sample_Id, y = Abundance, fill = !!sym(lineage))) +
      geom_bar(stat = "identity", position = "stack", color = "black", linewidth=0.2) +
      theme_linedraw() +
      theme(
        axis.text.x = element_text(size = 10, face = "bold"),
        axis.text.y = element_text(size = 10, face = "bold"),
        axis.title.x = element_text(size = 10, face = "bold"),
        axis.title.y = element_text(size = 10, face = "bold"),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(face = "bold"),
        legend.position = "right",
        title = element_text(size = 10, face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      ) +
      guides(fill = guide_legend(ncol = 4)) +
      labs(x = "BARCODE", y="% RELATIVE ABUNDANCE", fill = lineage) + 
      scale_fill_manual(values = fill_colors) +
      scale_y_continuous(expand = c(0, 0), limits = c(0,102))

      plot_taxa_stacked(stacked_barplot)

      ggplotly(stacked_barplot)    
    
    }

    else if(route()=="Offline")
    {
      req(result_dir_val(), input$taxa, input$top_taxa)

      lineage <- input$taxa

      top_n <- input$top_taxa

      dir <- result_dir_val()

      stacked_df <- cohort_offline_analysis(dir, lineage)$rel_abundance_data

      required_col <- which(colnames(stacked_df)==lineage)

      stacked_df[,required_col] <- as.character(stacked_df[,required_col])

      stacked_df <- stacked_df %>% 
        pivot_longer(cols = -all_of(required_col), names_to = "Sample_Id", values_to = "Abundance") %>% 
        filter(Abundance > 0)

      stacked_df <- stacked_df %>% group_by(Sample_Id) %>%
        arrange(desc(Abundance)) %>% dplyr::slice(1:ifelse(n()<top_n,n(),top_n))

      stacked_df$Sample_Id <- gsub("barcode","", stacked_df$Sample_Id)

      temp <- stacked_df %>% group_by(Sample_Id) %>% summarise(Abundance = 100-sum(Abundance)) %>% mutate(Species="Others") %>% 
        dplyr::select(Species, Sample_Id, Abundance)

      colnames(temp) <- colnames(stacked_df)

      stacked_df <- rbind(stacked_df, temp) %>% arrange(Sample_Id)

      total_abundance <- stacked_df %>%
        group_by(!!sym(lineage)) %>%
        summarise(Total_Abundance = sum(Abundance)) %>%
        arrange(Total_Abundance)

      stacked_df[[required_col]] <- factor(stacked_df[[required_col]], 
                                  levels = c("Others", total_abundance[[required_col]][total_abundance[[required_col]] != "Others"]))

      stacked_df <- stacked_df %>% filter(!!sym(lineage)!="Unclassified")

      fill_colors <- c("Others" = "#D3D3D3",setNames(viridis_pal(option = "D")(
        length(levels(stacked_df[[required_col]]))-1), 
        levels(stacked_df[[required_col]])[levels(stacked_df[[required_col]]) != "Others"]))

      stacked_barplot <- ggplot(stacked_df, aes(x = Sample_Id, y = Abundance, fill = !!sym(lineage))) +
      geom_bar(stat = "identity", position = "stack", color = "black", linewidth=0.2) +
      theme_linedraw() +
      theme(
        axis.text.x = element_text(size = 10, face = "bold"),
        axis.text.y = element_text(size = 10, face = "bold"),
        axis.title.x = element_text(size = 10, face = "bold"),
        axis.title.y = element_text(size = 10, face = "bold"),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(face = "bold"),
        legend.position = "right",
        title = element_text(size = 10, face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      ) +
      guides(fill = guide_legend(ncol = 4)) +
      labs(x = "BARCODE", y="% RELATIVE ABUNDANCE", fill = lineage) + 
      scale_fill_manual(values = fill_colors) +
      scale_y_continuous(expand = c(0, 0), limits = c(0,102))

      plot_taxa_stacked(stacked_barplot)

      ggplotly(stacked_barplot)
    }

  })

  output$plot_boxplot <- renderPlot({
    
    if(route()=="Realtime")
    {
      req(cohort_analysis_list(), input$taxa, classified_samples_list(), input_data_reactive())

      sample_list <- classified_samples_list()$Barcode

      lineage <- input$taxa
      
      sample_metadata <- input_data_reactive()$data
      
      alpha_diversity_data <- cohort_realtime_analysis(cohort_analysis_list(), sample_list, lineage)$rel_abundance_filtered_data

      required_col <- which(colnames(alpha_diversity_data)==lineage)
      
      alpha_diversity_data[,required_col] <- as.character(alpha_diversity_data[,required_col])
        
      alpha_diversity_data <- alpha_diversity_data %>%
        pivot_longer(cols = -(required_col), names_to = "Sample_Id", values_to = "Abundance") %>%
        filter(Abundance > 0)
      
      alpha_diversity_data <- alpha_diversity_data %>% 
        group_by(Sample_Id) %>% 
        summarise(Shannon = diversity(Abundance, index = "shannon"),
                  Simpson = diversity(Abundance, index = "simpson"))
      
      alpha_diversity_data <- inner_join(alpha_diversity_data, sample_metadata)
      
      alpha_diversity_data$Group <- factor(alpha_diversity_data$Group)
      
      alpha_diversity_data <- pivot_longer(alpha_diversity_data, cols = c(Shannon, Simpson), 
                                          names_to = "Diversity", values_to = "Value")
      
      alpha_div_p <- compare_means(Value~Group, data = alpha_diversity_data, method = "wilcox",
                                  p.adjust.method = "holm", group.by = "Diversity")
      
      shannon_plot <- alpha_diversity_data %>% filter(Diversity == "Shannon") %>%
        ggplot(aes(x=Group, y=Value, color=Group)) +
        geom_boxplot() +
        scale_color_manual(values = pal_aaas("default")(length(levels(alpha_diversity_data$Group)))) +
        theme_linedraw() +
        labs(y= "Alpha Diversity", x = "") +
        stat_pvalue_manual(subset(alpha_div_p, Diversity=="Shannon"), label = "p.signif", y.position = max(
          subset(alpha_diversity_data, Diversity=="Shannon")$Value) + 0.1, hide.ns = "p.adj", step.increase = 0.1,
          tip.length = 0.02, bracket.size = 0.8, size = 8) +
        guides(color = guide_legend(title = "Shannon", title.position = "top")) +
        theme(
          axis.title.x = element_text(size = 14, face = "bold"),
          axis.title.y = element_text(size = 14, face = "bold"),
          strip.text.x = element_text(size = 14, face = "bold"),
          axis.text.y=element_text(size=14, face = "bold"),
          axis.text.x= element_blank(),
          axis.ticks.x = element_blank(),
          legend.title=element_text(colour="black",size=10, face = "bold"),
          legend.text=element_text(colour="black", size=10, face = "bold"),
          plot.title = element_text(hjust = 0.5,size = rel(2)),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()
        )
      
      simpson_plot <- alpha_diversity_data %>% filter(Diversity == "Simpson") %>%
        ggplot(aes(x=Group, y=Value, color=Group)) +
        geom_boxplot() +
        scale_color_manual(values = pal_aaas("default")(length(levels(alpha_diversity_data$Group)))) + 
        theme_linedraw() +
        labs(y= "Alpha diversity", x = "") +
        stat_pvalue_manual(subset(alpha_div_p, Diversity=="Simpson"), y.position = max(
          subset(alpha_diversity_data, Diversity=="Simpson")$Value) + 0.01, label = "p.signif",
          hide.ns = "p.adj", step.increase = 0.1, tip.length = 0.02, bracket.size = 0.8, size = 8) +
        guides(color = guide_legend(title = "Simpson", title.position = "top")) +
        theme(
          axis.title.x = element_text(size = 14, face = "bold"),
          axis.title.y = element_text(size = 14, face = "bold"),
          strip.text.x = element_text(size = 14, face = "bold"),
          axis.text.y=element_text(size=14, face = "bold"),
          axis.text.x= element_blank(),
          axis.ticks.x = element_blank(),
          legend.title=element_text(colour="black",size=10, face = "bold"),
          legend.text=element_text(colour="black", size=10, face = "bold"),
          plot.title = element_text(hjust = 0.5,size = rel(2)),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()
        )
      
      
      diversity_facet <- as_ggplot(grid.grabExpr(grid.arrange(shannon_plot,
                                          simpson_plot, ncol=2)))

      plot_diversity_box(diversity_facet)

      diversity_facet
    
    }

    else if(route()=="Example")
    {
      req(input$taxa)
      
      file_list <- gsub(".txt", "", list.files("Example/")[grep("\\_final_blast_result.txt$",
                                                                                                list.files("Example/"))])
      lineage <- input$taxa

      sample_metadata <- input_data_reactive()$data
      
      alpha_diversity_data <- example_analysis("Example", file_list, lineage)$rel_abundance_filtered_data

      required_col <- which(colnames(alpha_diversity_data)==lineage)
      
      alpha_diversity_data[,required_col] <- as.character(alpha_diversity_data[,required_col])
        
      alpha_diversity_data <- alpha_diversity_data %>%
        pivot_longer(cols = -(required_col), names_to = "Sample_Id", values_to = "Abundance") %>%
        filter(Abundance > 0)
      
      alpha_diversity_data <- alpha_diversity_data %>% 
        group_by(Sample_Id) %>% 
        summarise(Shannon = diversity(Abundance, index = "shannon"),
                  Simpson = diversity(Abundance, index = "simpson"))
      
      alpha_diversity_data <- inner_join(alpha_diversity_data, sample_metadata)
      
      alpha_diversity_data$Group <- factor(alpha_diversity_data$Group)
      
      alpha_diversity_data <- pivot_longer(alpha_diversity_data, cols = c(Shannon, Simpson), 
                                          names_to = "Diversity", values_to = "Value")
      
      alpha_div_p <- compare_means(Value~Group, data = alpha_diversity_data, method = "wilcox",
                                  p.adjust.method = "holm", group.by = "Diversity")
      
      shannon_plot <- alpha_diversity_data %>% filter(Diversity == "Shannon") %>%
        ggplot(aes(x=Group, y=Value, color=Group)) +
        geom_boxplot() +
        geom_jitter(shape = 16, position = position_jitter(0.2)) +
        scale_color_manual(values = pal_aaas("default")(length(levels(alpha_diversity_data$Group)))) +
        theme_linedraw() +
        labs(y= "Alpha Diversity", x = "") +
        stat_pvalue_manual(subset(alpha_div_p, Diversity=="Shannon"), label = "p.signif", y.position = max(
          subset(alpha_diversity_data, Diversity=="Shannon")$Value) + 0.1, hide.ns = "p.adj", step.increase = 0.1,
          tip.length = 0.02, bracket.size = 0.8, size = 8) +
        guides(color = guide_legend(title = "Shannon", title.position = "top")) +
        theme(
          axis.title.x = element_text(size = 14, face = "bold"),
          axis.title.y = element_text(size = 14, face = "bold"),
          strip.text.x = element_text(size = 14, face = "bold"),
          axis.text.y=element_text(size=14, face = "bold"),
          axis.text.x= element_blank(),
          axis.ticks.x = element_blank(),
          legend.title=element_text(colour="black",size=10, face = "bold"),
          legend.text=element_text(colour="black", size=10, face = "bold"),
          plot.title = element_text(hjust = 0.5,size = rel(2)),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()
        )
      
      simpson_plot <- alpha_diversity_data %>% filter(Diversity == "Simpson") %>%
        ggplot(aes(x=Group, y=Value, color=Group)) +
        geom_boxplot() +
        geom_jitter(shape = 16, position = position_jitter(0.2)) +
        scale_color_manual(values = pal_aaas("default")(length(levels(alpha_diversity_data$Group)))) + 
        theme_linedraw() +
        labs(y= "Alpha diversity", x = "") +
        stat_pvalue_manual(subset(alpha_div_p, Diversity=="Simpson"), y.position = max(
          subset(alpha_diversity_data, Diversity=="Simpson")$Value) + 0.01, label = "p.signif",
          hide.ns = "p.adj", step.increase = 0.1, tip.length = 0.02, bracket.size = 0.8, size = 8) +
        guides(color = guide_legend(title = "Simpson", title.position = "top")) +
        theme(
          axis.title.x = element_text(size = 14, face = "bold"),
          axis.title.y = element_text(size = 14, face = "bold"),
          strip.text.x = element_text(size = 14, face = "bold"),
          axis.text.y=element_text(size=14, face = "bold"),
          axis.text.x= element_blank(),
          axis.ticks.x = element_blank(),
          legend.title=element_text(colour="black",size=10, face = "bold"),
          legend.text=element_text(colour="black", size=10, face = "bold"),
          plot.title = element_text(hjust = 0.5,size = rel(2)),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()
        )
      
      
      diversity_facet <- as_ggplot(grid.grabExpr(grid.arrange(shannon_plot,
                                          simpson_plot, ncol=2)))

      plot_diversity_box(diversity_facet)

      diversity_facet
    
    }
    
    else if(route() == "Offline")
    {
      req(result_dir_val(), input$taxa)

      dir <- result_dir_val()

      lineage <- input$taxa

      sample_metadata <- input_data_reactive()$data
      
      alpha_diversity_data <- cohort_offline_analysis(dir, lineage)$rel_abundance_filtered_data

      required_col <- which(colnames(alpha_diversity_data)==lineage)
      
      alpha_diversity_data[,required_col] <- as.character(alpha_diversity_data[,required_col])
        
      alpha_diversity_data <- alpha_diversity_data %>%
        pivot_longer(cols = -(required_col), names_to = "Sample_Id", values_to = "Abundance") %>%
        filter(Abundance > 0)
      
      alpha_diversity_data <- alpha_diversity_data %>% 
        group_by(Sample_Id) %>% 
        summarise(Shannon = diversity(Abundance, index = "shannon"),
                  Simpson = diversity(Abundance, index = "simpson"))
      
      alpha_diversity_data <- inner_join(alpha_diversity_data, sample_metadata)
      
      alpha_diversity_data$Group <- factor(alpha_diversity_data$Group)
      
      alpha_diversity_data <- pivot_longer(alpha_diversity_data, cols = c(Shannon, Simpson), 
                                          names_to = "Diversity", values_to = "Value")
      
      alpha_div_p <- compare_means(Value~Group, data = alpha_diversity_data, method = "wilcox",
                                  p.adjust.method = "holm", group.by = "Diversity")
      
      shannon_plot <- alpha_diversity_data %>% filter(Diversity == "Shannon") %>%
        ggplot(aes(x=Group, y=Value, color=Group)) +
        geom_boxplot() +
        scale_color_manual(values = pal_aaas("default")(length(levels(alpha_diversity_data$Group)))) +
        theme_linedraw() +
        labs(y= "Alpha Diversity", x = "") +
        stat_pvalue_manual(subset(alpha_div_p, Diversity=="Shannon"), label = "p.signif", y.position = max(
          subset(alpha_diversity_data, Diversity=="Shannon")$Value) + 0.1, hide.ns = "p.adj", step.increase = 0.1,
          tip.length = 0.02, bracket.size = 0.8, size = 8) +
        guides(color = guide_legend(title = "Shannon", title.position = "top")) +
        theme(
          axis.title.x = element_text(size = 14, face = "bold"),
          axis.title.y = element_text(size = 14, face = "bold"),
          strip.text.x = element_text(size = 14, face = "bold"),
          axis.text.y=element_text(size=14, face = "bold"),
          axis.text.x= element_blank(),
          axis.ticks.x = element_blank(),
          legend.title=element_text(colour="black",size=10, face = "bold"),
          legend.text=element_text(colour="black", size=10, face = "bold"),
          plot.title = element_text(hjust = 0.5,size = rel(2)),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()
        )
      
      simpson_plot <- alpha_diversity_data %>% filter(Diversity == "Simpson") %>%
        ggplot(aes(x=Group, y=Value, color=Group)) +
        geom_boxplot() +
        scale_color_manual(values = pal_aaas("default")(length(levels(alpha_diversity_data$Group)))) + 
        theme_linedraw() +
        labs(y= "Alpha diversity", x = "") +
        stat_pvalue_manual(subset(alpha_div_p, Diversity=="Simpson"), y.position = max(
          subset(alpha_diversity_data, Diversity=="Simpson")$Value) + 0.01, label = "p.signif",
          hide.ns = "p.adj", step.increase = 0.1, tip.length = 0.02, bracket.size = 0.8, size = 8) +
        guides(color = guide_legend(title = "Simpson", title.position = "top")) +
        theme(
          axis.title.x = element_text(size = 14, face = "bold"),
          axis.title.y = element_text(size = 14, face = "bold"),
          strip.text.x = element_text(size = 14, face = "bold"),
          axis.text.y=element_text(size=14, face = "bold"),
          axis.text.x= element_blank(),
          axis.ticks.x = element_blank(),
          legend.title=element_text(colour="black",size=10, face = "bold"),
          legend.text=element_text(colour="black", size=10, face = "bold"),
          plot.title = element_text(hjust = 0.5,size = rel(2)),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()
        )
      
      
      diversity_facet <- as_ggplot(grid.grabExpr(grid.arrange(shannon_plot,
                                          simpson_plot, ncol=2)))

      plot_diversity_box(diversity_facet)

      diversity_facet
    
    }
  
  }, height = 500)

  output$plot_pcoa <- renderPlot({

    if(route()=="Realtime")
    {
      req(cohort_analysis_list(), input$taxa, classified_samples_list(), input_data_reactive())

      sample_list <- classified_samples_list()$Barcode

      lineage <- input$taxa
      
      sample_metadata <- input_data_reactive()$data

      matrix <- cohort_realtime_analysis(cohort_analysis_list(), sample_list, lineage)$abundance_matrix

      pcoa_dist <- wcmdscale(vegdist(t(matrix), method = "aitchison", pseudocount=0.5), k=2, eig = TRUE)

      pcoa_df <- pcoa_dist$points[,1:2] %>% as.data.frame()
    
      pcoa_eigenvalues <- pcoa_dist$eig

      pcoa.var <- round(pcoa_eigenvalues/sum(pcoa_eigenvalues)*100, 1)

      pcoa_df$Sample_Id <- rownames(pcoa_df)

      pcoa_df <- inner_join(pcoa_df, sample_metadata)

      pcoa_df$Group <- factor(pcoa_df$Group)
    
      colnames(pcoa_df) <- c("Axis.1", "Axis.2", "Sample_Id", "Group")

      pcoa_plot <- ggplot(data = pcoa_df, aes(x = Axis.1, y = Axis.2, color = Group)) +
        geom_point(size=5) +
        stat_ellipse(aes(colour = Group, fill = Group), level = 0.95, alpha = 0.25, geom = "polygon") +
        scale_color_manual(values = pal_aaas("default")(length(levels(pcoa_df$Group)))) +
        scale_fill_manual(values = pal_aaas("default")(length(levels(pcoa_df$Group)))) +
        theme_linedraw() +
        xlab(paste0("PC1 (", pcoa.var[1], "%", ")")) +
        ylab(paste0("PC2 (", pcoa.var[2], "%", ")")) +
        geom_vline(xintercept = 0, color = "black", linetype = "dashed") +
        geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
        theme(
          axis.text.x = element_text(size = 15, face = "bold"),
          axis.text.y = element_text(size = 15, face = "bold"),
          legend.text = element_text(size = 15, face = "bold"),
          axis.title.x = element_text(size = 15, face = "bold"),
          axis.title.y = element_text(size = 15, face = "bold"),
          legend.title = element_text(size = 15, face = "bold"),
          title = element_text(size = 15, face = "bold")
        ) +
        ggtitle("Principal Coordination Analysis (PCoA) ordination plot using Aitchison Distance")

      plot_pcoa_dot(pcoa_plot)

      pcoa_plot
    
    }

    else if(route()=="Example")
    {
      
      req(input$taxa)

      file_list <- gsub(".txt", "", list.files("Example/")[grep("\\_final_blast_result.txt$",
                                                                                                list.files("Example/"))])
                      
      lineage <- input$taxa

      sample_metadata <- input_data_reactive()$data

      matrix <- example_analysis("Example", file_list, lineage)$abundance_matrix

      pcoa_dist <- wcmdscale(vegdist(t(matrix), method = "aitchison", pseudocount=0.5), k=2, eig = TRUE)

      pcoa_df <- pcoa_dist$points[,1:2] %>% as.data.frame()
    
      pcoa_eigenvalues <- pcoa_dist$eig

      pcoa.var <- round(pcoa_eigenvalues/sum(pcoa_eigenvalues)*100, 1)

      pcoa_df$Sample_Id <- rownames(pcoa_df)

      pcoa_df <- inner_join(pcoa_df, sample_metadata)

      pcoa_df$Group <- factor(pcoa_df$Group)
    
      colnames(pcoa_df) <- c("Axis.1", "Axis.2", "Sample_Id", "Group")

      pcoa_plot <- ggplot(data = pcoa_df, aes(x = Axis.1, y = Axis.2, color = Group)) +
        geom_point(size=5) +
        stat_ellipse(aes(colour = Group, fill = Group), level = 0.95, alpha = 0.25, geom = "polygon") +
        scale_color_manual(values = pal_aaas("default")(length(levels(pcoa_df$Group)))) +
        scale_fill_manual(values = pal_aaas("default")(length(levels(pcoa_df$Group)))) +
        theme_linedraw() +
        xlab(paste0("PC1 (", pcoa.var[1], "%", ")")) +
        ylab(paste0("PC2 (", pcoa.var[2], "%", ")")) +
        geom_vline(xintercept = 0, color = "black", linetype = "dashed") +
        geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
        theme(
          axis.text.x = element_text(size = 15, face = "bold"),
          axis.text.y = element_text(size = 15, face = "bold"),
          legend.text = element_text(size = 15, face = "bold"),
          axis.title.x = element_text(size = 15, face = "bold"),
          axis.title.y = element_text(size = 15, face = "bold"),
          legend.title = element_text(size = 15, face = "bold"),
          title = element_text(size = 15, face = "bold")
        ) +
        ggtitle("Principal Coordination Analysis (PCoA) ordination plot using Aitchison Distance")

      plot_pcoa_dot(pcoa_plot)

      pcoa_plot

    }

    else if(route() == "Offline")
    {
      req(result_dir_val(), input$taxa)

      dir <- result_dir_val()

      lineage <- input$taxa

      sample_metadata <- input_data_reactive()$data

      matrix <- cohort_offline_analysis(dir, lineage)$abundance_matrix

      pcoa_dist <- wcmdscale(vegdist(t(matrix), method = "aitchison", pseudocount=0.5), k=2, eig = TRUE)

      pcoa_df <- pcoa_dist$points[,1:2] %>% as.data.frame()
    
      pcoa_eigenvalues <- pcoa_dist$eig

      pcoa.var <- round(pcoa_eigenvalues/sum(pcoa_eigenvalues)*100, 1)

      pcoa_df$Sample_Id <- rownames(pcoa_df)

      pcoa_df <- inner_join(pcoa_df, sample_metadata)

      pcoa_df$Group <- factor(pcoa_df$Group)
    
      colnames(pcoa_df) <- c("Axis.1", "Axis.2", "Sample_Id", "Group")

      pcoa_plot <- ggplot(data = pcoa_df, aes(x = Axis.1, y = Axis.2, color = Group)) +
        geom_point(size=5) +
        stat_ellipse(aes(colour = Group, fill = Group), level = 0.95, alpha = 0.25, geom = "polygon") +
        scale_color_manual(values = pal_aaas("default")(length(levels(pcoa_df$Group)))) +
        scale_fill_manual(values = pal_aaas("default")(length(levels(pcoa_df$Group)))) +
        theme_linedraw() +
        xlab(paste0("PC1 (", pcoa.var[1], "%", ")")) +
        ylab(paste0("PC2 (", pcoa.var[2], "%", ")")) +
        geom_vline(xintercept = 0, color = "black", linetype = "dashed") +
        geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
        theme(
          axis.text.x = element_text(size = 15, face = "bold"),
          axis.text.y = element_text(size = 15, face = "bold"),
          legend.text = element_text(size = 15, face = "bold"),
          axis.title.x = element_text(size = 15, face = "bold"),
          axis.title.y = element_text(size = 15, face = "bold"),
          legend.title = element_text(size = 15, face = "bold"),
          title = element_text(size = 15, face = "bold")
        ) +
        ggtitle("Principal Coordination Analysis (PCoA) ordination plot using Aitchison Distance")

      plot_pcoa_dot(pcoa_plot)

      pcoa_plot
    
    }
  
  }, height = 500) 

  output$plot_nmds <- renderPlot({
    
    if(route()=="Realtime")
    {
      req(cohort_analysis_list(), input$taxa, classified_samples_list(), input_data_reactive())

      sample_list <- classified_samples_list()$Barcode

      lineage <- input$taxa
      
      sample_metadata <- input_data_reactive()$data

      matrix <- cohort_realtime_analysis(cohort_analysis_list(), sample_list, lineage)$rel_abundance_matrix
      
      nmds_dist <- metaMDS(t(matrix), distance = "bray", trymax = 100)

      nmds_df <- as.data.frame(nmds_dist$points)

      nmds_df$Group <- sample_metadata$Group

      nmds_df$Group <- factor(nmds_df$Group)

      nmds_plot <- ggplot(data = nmds_df, aes(x = MDS1, y = MDS2, color = Group)) +
        geom_point(size=5) +
        stat_ellipse(aes(colour = Group, fill = Group), level = 0.95, alpha = 0.25, geom = "polygon") +
        scale_color_manual(values = pal_aaas("default")(length(levels(nmds_df$Group)))) +
        scale_fill_manual(values = pal_aaas("default")(length(levels(nmds_df$Group)))) +
        theme_linedraw() +
        theme(
          axis.text.x = element_text(size = 15, face = "bold"),
          axis.text.y = element_text(size = 15, face = "bold"),
          legend.text = element_text(size = 15, face = "bold"),
          axis.title.x = element_text(size = 15, face = "bold"),
          axis.title.y = element_text(size = 15, face = "bold"),
          legend.title = element_text(size = 15, face = "bold"),
          title = element_text(size = 15, face = "bold")
        ) +
        geom_vline(xintercept = 0, color = "black", linetype = "dashed") +
        geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
        ggtitle("Non-metric MultiDimensional Scaling (NMDS) ordination plot with Bray-Curtis Distance")
      
      plot_nmds_dot(nmds_plot)

      nmds_plot
    
    }

    else if(route()=="Example")
    {
      req(input$taxa)
      
      file_list <- gsub(".txt", "", list.files("Example/")[grep("\\_final_blast_result.txt$",
                                                                                                list.files("Example/"))])
      lineage <- input$taxa

      sample_metadata <- input_data_reactive()$data

      matrix <- example_analysis("Example", file_list, lineage)$rel_abundance_matrix

      nmds_dist <- metaMDS(t(matrix), distance = "bray", trymax = 100)

      nmds_df <- as.data.frame(nmds_dist$points)

      nmds_df$Group <- sample_metadata$Group

      nmds_df$Group <- factor(nmds_df$Group)

      nmds_plot <- ggplot(data = nmds_df, aes(x = MDS1, y = MDS2, color = Group)) +
        geom_point(size=5) +
        stat_ellipse(aes(colour = Group, fill = Group), level = 0.95, alpha = 0.25, geom = "polygon") +
        scale_color_manual(values = pal_aaas("default")(length(levels(nmds_df$Group)))) +
        scale_fill_manual(values = pal_aaas("default")(length(levels(nmds_df$Group)))) +
        theme_linedraw() +
        theme(
          axis.text.x = element_text(size = 15, face = "bold"),
          axis.text.y = element_text(size = 15, face = "bold"),
          legend.text = element_text(size = 15, face = "bold"),
          axis.title.x = element_text(size = 15, face = "bold"),
          axis.title.y = element_text(size = 15, face = "bold"),
          legend.title = element_text(size = 15, face = "bold"),
          title = element_text(size = 15, face = "bold")
        ) +
        geom_vline(xintercept = 0, color = "black", linetype = "dashed") +
        geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
        ggtitle("Non-metric MultiDimensional Scaling (NMDS) ordination plot with Bray-Curtis Distance")
      
      plot_nmds_dot(nmds_plot)

      nmds_plot

    }

    else if(route() == "Offline")
    {
      req(result_dir_val(), input$taxa)

      dir <- result_dir_val()

      lineage <- input$taxa

      sample_metadata <- input_data_reactive()$data

      matrix <- cohort_offline_analysis(dir, lineage)$rel_abundance_matrix

      nmds_dist <- metaMDS(t(matrix), distance = "bray", trymax = 100)

      nmds_df <- as.data.frame(nmds_dist$points)

      nmds_df$Group <- sample_metadata$Group

      nmds_df$Group <- factor(nmds_df$Group)

      nmds_plot <- ggplot(data = nmds_df, aes(x = MDS1, y = MDS2, color = Group)) +
        geom_point(size=5) +
        stat_ellipse(aes(colour = Group, fill = Group), level = 0.95, alpha = 0.25, geom = "polygon") +
        scale_color_manual(values = pal_aaas("default")(length(levels(nmds_df$Group)))) +
        scale_fill_manual(values = pal_aaas("default")(length(levels(nmds_df$Group)))) +
        theme_linedraw() +
        theme(
          axis.text.x = element_text(size = 15, face = "bold"),
          axis.text.y = element_text(size = 15, face = "bold"),
          legend.text = element_text(size = 15, face = "bold"),
          axis.title.x = element_text(size = 15, face = "bold"),
          axis.title.y = element_text(size = 15, face = "bold"),
          legend.title = element_text(size = 15, face = "bold"),
          title = element_text(size = 15, face = "bold")
        ) +
        geom_vline(xintercept = 0, color = "black", linetype = "dashed") +
        geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
        ggtitle("Non-metric MultiDimensional Scaling (NMDS) ordination plot with Bray-Curtis Distance")
      
      plot_nmds_dot(nmds_plot)

      nmds_plot

    }
    
  }, height = 500)

  output$plot_pca <- renderPlot({
    
    if(route()=="Realtime")
    {
      req(cohort_analysis_list(), input$taxa, classified_samples_list(), input_data_reactive())

      sample_list <- classified_samples_list()$Barcode

      lineage <- input$taxa

      sample_metadata <- input_data_reactive()$data

      matrix <- cohort_realtime_analysis(cohort_analysis_list(), sample_list, lineage)$abundance_matrix

      clr_transformed_abundance_matrix <- clr(matrix+0.5) %>% as.data.frame()

      pca <- prcomp(t(clr_transformed_abundance_matrix), scale. = FALSE, center = TRUE,
              retx = TRUE)

      pca.var <- pca$sdev^2

      pca.var.per <- round(pca.var/sum(pca.var)*100, 1)

      pca_df <- data.frame(Sample_Id = rownames(pca$x), X=pca$x[,1], Y=pca$x[,2])

      pca_df <- inner_join(pca_df, sample_metadata)

      pca_df$Group <- factor(pca_df$Group)

      pca_plot <- ggplot(data = pca_df, aes(x = X, y = Y, color = Group)) +
        geom_point(size=5) +
        stat_ellipse(aes(colour = Group, fill = Group), level = 0.95, alpha = 0.25, geom = "polygon") +
        xlab(paste0("PC1 (", pca.var.per[1], "%", ")")) +
        ylab(paste0("PC2 (", pca.var.per[2], "%", ")")) +
        scale_color_manual(values = pal_aaas("default")(length(levels(pca_df$Group)))) +
        scale_fill_manual(values = pal_aaas("default")(length(levels(pca_df$Group)))) +
        theme_linedraw() +
        theme(
          axis.text.x = element_text(size = 15, face = "bold"),
          axis.text.y = element_text(size = 15, face = "bold"),
          legend.text = element_text(size = 15, face = "bold"),
          axis.title.x = element_text(size = 15, face = "bold"),
          axis.title.y = element_text(size = 15, face = "bold"),
          legend.title = element_text(size = 15, face = "bold"),
          title = element_text(size = 15, face = "bold")
        ) +
        geom_vline(xintercept = 0, color = "black", linetype = "dashed") +
        geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
        ggtitle("Principal Component Analysis (PCA) Plot using CLR Transformed Data")

      plot_pca_dot(pca_plot)
      
      pca_plot
    
    }

    else if(route()=="Example")
    {
      req(input$taxa)

      file_list <- gsub(".txt", "", list.files("Example/")[grep("\\_final_blast_result.txt$",
                                                                                                list.files("Example/"))])                
      lineage <- input$taxa

      sample_metadata <- input_data_reactive()$data

      matrix <- example_analysis("Example", file_list, lineage)$abundance_matrix

      clr_transformed_abundance_matrix <- clr(matrix+0.5) %>% as.data.frame()

      pca <- prcomp(t(clr_transformed_abundance_matrix), scale. = FALSE, center = TRUE,
              retx = TRUE)

      pca.var <- pca$sdev^2

      pca.var.per <- round(pca.var/sum(pca.var)*100, 1)

      pca_df <- data.frame(Sample_Id = rownames(pca$x), X=pca$x[,1], Y=pca$x[,2])

      pca_df <- inner_join(pca_df, sample_metadata)

      pca_df$Group <- factor(pca_df$Group)

      pca_plot <- ggplot(data = pca_df, aes(x = X, y = Y, color = Group)) +
        geom_point(size=5) +
        stat_ellipse(aes(colour = Group, fill = Group), level = 0.95, alpha = 0.25, geom = "polygon") +
        xlab(paste0("PC1 (", pca.var.per[1], "%", ")")) +
        ylab(paste0("PC2 (", pca.var.per[2], "%", ")")) +
        scale_color_manual(values = pal_aaas("default")(length(levels(pca_df$Group)))) +
        scale_fill_manual(values = pal_aaas("default")(length(levels(pca_df$Group)))) +
        theme_linedraw() +
        theme(
          axis.text.x = element_text(size = 15, face = "bold"),
          axis.text.y = element_text(size = 15, face = "bold"),
          legend.text = element_text(size = 15, face = "bold"),
          axis.title.x = element_text(size = 15, face = "bold"),
          axis.title.y = element_text(size = 15, face = "bold"),
          legend.title = element_text(size = 15, face = "bold"),
          title = element_text(size = 15, face = "bold")
        ) +
        geom_vline(xintercept = 0, color = "black", linetype = "dashed") +
        geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
        ggtitle("Principal Component Analysis (PCA) Plot using CLR Transformed Data")

      plot_pca_dot(pca_plot)
      
      pca_plot  
    
    }

    else if(route() == "Offline")
    {
      req(result_dir_val(), input$taxa)

      dir <- result_dir_val()

      lineage <- input$taxa

      sample_metadata <- input_data_reactive()$data

      matrix <- cohort_offline_analysis(dir, lineage)$abundance_matrix

      clr_transformed_abundance_matrix <- clr(matrix+0.5) %>% as.data.frame()

      pca <- prcomp(t(clr_transformed_abundance_matrix), scale. = FALSE, center = TRUE,
              retx = TRUE)

      pca.var <- pca$sdev^2

      pca.var.per <- round(pca.var/sum(pca.var)*100, 1)

      pca_df <- data.frame(Sample_Id = rownames(pca$x), X=pca$x[,1], Y=pca$x[,2])

      pca_df <- inner_join(pca_df, sample_metadata)

      pca_df$Group <- factor(pca_df$Group)

      pca_plot <- ggplot(data = pca_df, aes(x = X, y = Y, color = Group)) +
        geom_point(size=5) +
        stat_ellipse(aes(colour = Group, fill = Group), level = 0.95, alpha = 0.25, geom = "polygon") +
        xlab(paste0("PC1 (", pca.var.per[1], "%", ")")) +
        ylab(paste0("PC2 (", pca.var.per[2], "%", ")")) +
        scale_color_manual(values = pal_aaas("default")(length(levels(pca_df$Group)))) +
        scale_fill_manual(values = pal_aaas("default")(length(levels(pca_df$Group)))) +
        theme_linedraw() +
        theme(
          axis.text.x = element_text(size = 15, face = "bold"),
          axis.text.y = element_text(size = 15, face = "bold"),
          legend.text = element_text(size = 15, face = "bold"),
          axis.title.x = element_text(size = 15, face = "bold"),
          axis.title.y = element_text(size = 15, face = "bold"),
          legend.title = element_text(size = 15, face = "bold"),
          title = element_text(size = 15, face = "bold")
        ) +
        geom_vline(xintercept = 0, color = "black", linetype = "dashed") +
        geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
        ggtitle("Principal Component Analysis (PCA) Plot using CLR Transformed Data")

      plot_pca_dot(pca_plot)
      
      pca_plot

    }

  }, height = 500)

  output$permanova_data <- DT::renderDataTable({
    
    if(route()=="Realtime")
    {
      req(cohort_analysis_list(), input$taxa, classified_samples_list(), input_data_reactive())

      sample_list <- classified_samples_list()$Barcode

      lineage <- input$taxa

      sample_metadata <- input_data_reactive()$data
      
      matrix <- cohort_realtime_analysis(cohort_analysis_list(), sample_list, lineage)$abundance_matrix
      
      perm_dist <- vegdist(t(matrix), method = "aitchison", pseudocount=0.5)

      permanova_res <- pairwise.adonis(perm_dist, as.factor(sample_metadata$Group), p.adjust.m = "BH")

      permanova_res <- permanova_res %>% dplyr::select(c(pairs, R2, p.value, p.adjusted))

      colnames(permanova_res) <- c("Pair", "R2", "P", "Padj")

      table_permanova(permanova_res)

      permanova_res
    
    }

    else if(route()=="Example")
    {
      req(input$taxa)

      file_list <- gsub(".txt", "", list.files("Example/")[grep("\\_final_blast_result.txt$",
                                                                                                list.files("Example/"))])
      lineage <- input$taxa

      sample_metadata <- input_data_reactive()$data
      
      matrix <- example_analysis("Example", file_list, lineage)$abundance_matrix
      
      perm_dist <- vegdist(t(matrix), method = "aitchison", pseudocount=0.5)

      permanova_res <- pairwise.adonis(perm_dist, as.factor(sample_metadata$Group), p.adjust.m = "BH")

      permanova_res <- permanova_res %>% dplyr::select(c(pairs, R2, p.value, p.adjusted))

      colnames(permanova_res) <- c("Pair", "R2", "P", "Padj")

      table_permanova(permanova_res)

      permanova_res

    }

    else if(route() == "Offline")
    {
      req(result_dir_val(), input$taxa)

      dir <- result_dir_val()

      lineage <- input$taxa

      sample_metadata <- input_data_reactive()$data
      
      matrix <- cohort_offline_analysis(dir, lineage)$abundance_matrix
      
      perm_dist <- vegdist(t(matrix), method = "aitchison", pseudocount=0.5)

      permanova_res <- pairwise.adonis(perm_dist, as.factor(sample_metadata$Group), p.adjust.m = "BH")

      permanova_res <- permanova_res %>% dplyr::select(c(pairs, R2, p.value, p.adjusted))

      colnames(permanova_res) <- c("Pair", "R2", "P", "Padj")

      table_permanova(permanova_res)

      permanova_res
    
    }
  
  }, height = 500)

  output$plot_heatmap <- renderPlot({
    
    if(route()=="Realtime")
    {
      req(cohort_analysis_list(), input$taxa, classified_samples_list(), input_data_reactive())

      sample_list <- classified_samples_list()$Barcode

      lineage <- input$taxa

      sample_metadata <- input_data_reactive()$data
      
      matrix <- cohort_realtime_analysis(cohort_analysis_list(), sample_list, lineage)$rel_abundance_filtered_matrix

      heatmap_df <- as.data.frame(log10(matrix+0.00001))

      sample_annotation <- data.frame(sample_metadata$Group)

      sample_metadata$Group <- factor(sample_metadata$Group)

      colnames(sample_annotation) <- "Group"

      sample_annotation$Group <- factor(sample_annotation$Group)

      col_list <- list(Group = setNames(pal_aaas("default")(length(levels(sample_metadata$Group))), 
      levels(sample_metadata$Group)))

      row_annotate <- rowAnnotation(
      df = sample_annotation,
      col = col_list, show_annotation_name = FALSE)

      row_dend <-  hclust(dist(t(heatmap_df)), method = "complete")

      column_dend <- hclust(dist(heatmap_df), method = "complete")

      tmp <- heatmap_df %>% pivot_longer(cols = 1:length(colnames(heatmap_df)), names_to = "Sample",
                                        values_to = "Abundance")

      col_fun <- colorRamp2(c(2, 0, -6), c("#1010fe", "#FFD700", "#FF1212"))

      heatmap_plot <- Heatmap(t(heatmap_df), heatmap_legend_param = list(title = expression("Log"[10]*"Relative Abundance"),
                                              title_gp = gpar(fontsize = 10)),
                            row_names_gp = gpar(fontsize = 10, fontface = "bold"),
                            cluster_rows = color_branches(row_dend),
                            cluster_columns = color_branches(column_dend),
                            show_column_names = FALSE,
                            show_row_names = TRUE,
                            right_annotation = row_annotate,
                            col = col_fun(seq(min(tmp$Abundance), max(tmp$Abundance))))

      heatmap_ggplot <- as_ggplot(grid.grabExpr(print(heatmap_plot)))

      plot_real_heatmap(heatmap_ggplot)

      heatmap_ggplot
    
    }

    else if(route()=="Example")
    {
      req(input$taxa)

      file_list <- gsub(".txt", "", list.files("Example/")[grep("\\_final_blast_result.txt$",
                                                                                                list.files("Example/"))])
      lineage <- input$taxa

      sample_metadata <- input_data_reactive()$data
      
      matrix <- example_analysis("Example", file_list, lineage)$rel_abundance_filtered_matrix

      heatmap_df <- as.data.frame(log10(matrix+0.00001))

      sample_annotation <- data.frame(sample_metadata$Group)

      sample_metadata$Group <- factor(sample_metadata$Group)

      colnames(sample_annotation) <- "Group"

      sample_annotation$Group <- factor(sample_annotation$Group)

      col_list <- list(Group = setNames(pal_aaas("default")(length(levels(sample_metadata$Group))), 
      levels(sample_metadata$Group)))

      row_annotate <- rowAnnotation(
      df = sample_annotation,
      col = col_list, show_annotation_name = FALSE)

      row_dend <-  hclust(dist(t(heatmap_df)), method = "complete")

      column_dend <- hclust(dist(heatmap_df), method = "complete")

      tmp <- heatmap_df %>% pivot_longer(cols = 1:length(colnames(heatmap_df)), names_to = "Sample",
                                        values_to = "Abundance")

      col_fun <- colorRamp2(c(2, 0, -6), c("#1010fe", "#FFD700", "#FF1212"))

      heatmap_plot <- Heatmap(t(heatmap_df), heatmap_legend_param = list(title = expression("Log"[10]*"Relative Abundance"),
                                              title_gp = gpar(fontsize = 10)),
                            row_names_gp = gpar(fontsize = 10, fontface = "bold"),
                            cluster_rows = color_branches(row_dend),
                            cluster_columns = color_branches(column_dend),
                            show_column_names = FALSE,
                            show_row_names = TRUE,
                            right_annotation = row_annotate,
                            col = col_fun(seq(min(tmp$Abundance), max(tmp$Abundance))))

      heatmap_ggplot <- as_ggplot(grid.grabExpr(print(heatmap_plot)))

      plot_real_heatmap(heatmap_ggplot)

      heatmap_ggplot
    
    }

    else if(route() == "Offline")
    {
      req(result_dir_val(), input$taxa)

      dir <- result_dir_val()

      lineage <- input$taxa

      sample_metadata <- input_data_reactive()$data
      
      matrix <- cohort_offline_analysis(dir, lineage)$rel_abundance_filtered_matrix

      heatmap_df <- as.data.frame(log10(matrix+0.00001))

      sample_annotation <- data.frame(sample_metadata$Group)

      sample_metadata$Group <- factor(sample_metadata$Group)

      colnames(sample_annotation) <- "Group"

      sample_annotation$Group <- factor(sample_annotation$Group)

      col_list <- list(Group = setNames(pal_aaas("default")(length(levels(sample_metadata$Group))), 
      levels(sample_metadata$Group)))

      row_annotate <- rowAnnotation(
      df = sample_annotation,
      col = col_list, show_annotation_name = FALSE)

      row_dend <-  hclust(dist(t(heatmap_df)), method = "complete")

      column_dend <- hclust(dist(heatmap_df), method = "complete")

      tmp <- heatmap_df %>% pivot_longer(cols = 1:length(colnames(heatmap_df)), names_to = "Sample",
                                        values_to = "Abundance")

      col_fun <- colorRamp2(c(2, 0, -6), c("#1010fe", "#FFD700", "#FF1212"))

      heatmap_plot <- Heatmap(t(heatmap_df), heatmap_legend_param = list(title = expression("Log"[10]*"Relative Abundance"),
                                              title_gp = gpar(fontsize = 10)),
                            row_names_gp = gpar(fontsize = 10, fontface = "bold"),
                            cluster_rows = color_branches(row_dend),
                            cluster_columns = color_branches(column_dend),
                            show_column_names = FALSE,
                            show_row_names = TRUE,
                            right_annotation = row_annotate,
                            col = col_fun(seq(min(tmp$Abundance), max(tmp$Abundance))))

      heatmap_ggplot <- as_ggplot(grid.grabExpr(print(heatmap_plot)))

      plot_real_heatmap(heatmap_ggplot)

      heatmap_ggplot
    
    }
  
  }, height = 500, width = 1000)

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

  output$download_stacked_barplot <- downloadHandler(
      filename = function() {
        req(input$taxa)
        lineage <- input$taxa
        paste("Stacked_bar_plot_", lineage,"_", Sys.Date(), ".png", sep = "")
      },
      content = function(file) {
        ggsave(file, plot_taxa_stacked(),
               width = 15, height = 10, units = "in", dpi = "retina", bg = "white")
      }
    )
    
    output$download_boxplot <- downloadHandler(
      filename = function() {
        req(input$taxa)
        lineage <- input$taxa
        paste("Alpha_Diversity_plot_", lineage, "_", Sys.Date(), ".png", sep = "")
      },
      content = function(file) {
        ggsave(file, plot_diversity_box(),
               width = 13.69, height = 8.27, units = "in", dpi = "retina", bg = "white")
        
      }
    )
    
    output$download_pcoa_plot <- downloadHandler(
      filename = function() {
        req(input$taxa)
        lineage <- input$taxa
        paste("PCoA_plot_", lineage, "_", Sys.Date(), ".png", sep = "")
      },
      content = function(file) {
        ggsave(file, plot_pcoa_dot(),
               width = 13.69, height = 8.27, units = "in", dpi = "retina", bg = "white")
      }
    )
    
    output$download_nmds_plot <- downloadHandler(
      filename = function() {
        req(input$taxa)
        lineage <- input$taxa
        paste("NMDS_plot_", lineage, "_", Sys.Date(), ".png", sep = "")
      },
      content = function(file) {
        
        ggsave(file, plot_nmds_dot(),
               width = 13.69, height = 8.27, units = "in", dpi = "retina", bg = "white")
        
      }
    )
    
    output$download_pca_plot <- downloadHandler(
      filename = function() {
        req(input$taxa)
        lineage <- input$taxa
        paste("PCA_plot_", lineage, "_", Sys.Date(), ".png", sep = "")
      },
      content = function(file) {
        ggsave(file, plot_pca_dot(),
               width = 13.69, height = 8.27, units = "in", dpi = "retina", bg = "white")
        }
    )
    
    output$download_heatmap <- downloadHandler(
      filename = function() {
        req(input$taxa)
        lineage <- input$taxa
        paste("HeatMap_", lineage, "_", Sys.Date(), ".png", sep = "")
      },
      content = function(file) {
        ggsave(file, plot_real_heatmap(),
               width = 13.69, height = 8.27, units = "in", dpi = "retina", bg = "white")
        
      }
    )

    output$download_permanova_csv <- downloadHandler(
      filename = function() {
        req(input$taxa)
        lineage <- input$taxa
        paste0("PERMANOVA_Result_", lineage, "_",Sys.Date(), ".tsv")
        },
      content = function(file) {
        write.csv(table_permanova(), file, row.names = FALSE, quote = FALSE, sep = "\t")
      }
    )


}

# Run the app
shinyApp(ui = ui, server = server)