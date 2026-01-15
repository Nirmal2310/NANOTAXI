source("packages.R")

plan(multisession, workers = 16)

ui <- navbarPage(title = div(class="titleimg",img(src="Nanotaxi.png", height="100%", width="12%"), "",
              style="position: absolute; top: 3px; left: 10px; background-color:white"),
              id = "main_navbar",
  tags$head(
    tags$style(HTML('.navbar-nav {
                            height: 25px;
                            padding-left: 120px;
                            backgroun-color:white}'),
               HTML('.navbar-brand {width: 20px; font-size:35px; text-align:left; background-color:white}'),
               HTML('.navbar-default {
                    background-color: white}'))
  ),
  useShinyjs(),
  source("UI/ui-tab-intro.R", local = TRUE)$value,
  source("UI/ui-tab-input.R", local = TRUE)$value,
  source("UI/ui-tab-cohort-analysis.R", local = TRUE)$value,
  source("UI/ui-tab-help.R", local = TRUE)$value,
  source("UI/ui-tab-conditions.R", local = TRUE)$value
)

conda_path <- system('if [ $(which conda | grep "condabin") ]; then conda_path=$(which conda | sed "s/\\/condabin.*$//g"); else conda_path=$(which conda | sed "s/\\/bin.*$//g"); fi && echo $conda_path', intern = TRUE)
 
Sys.setenv(RETICULATE_PYTHON = paste0(conda_path, "/envs/minknow_api/bin/python"))

server <- function(input, output, session) {
  observe({
    if (input$data_file_type == "upload") {
      insertTab(
        inputId = "main_navbar",
        source("UI/ui-tab-real-time-analysis.R", local = TRUE)$value,
        target = "cohort_tab",
        position = "before")
        }
        else {
          removeTab(inputId = "main_navbar", target = "realtime_tab")
        }
  })
  
  state <- reactiveVal()

  route <- reactiveVal()
  
  is_running <- reactiveVal(FALSE)

  cohort_run <- reactiveVal(FALSE)
  
  timer_10s <- reactiveTimer(10000)

  cohort_delay_done <- reactiveVal(FALSE)

  status_checked <- reactiveVal(FALSE)

  trigger <- reactiveVal(0)

  cohort_trigger <- reactiveVal(0)
  
  pores <- reactiveVal()
  
  passed_counts <- reactiveVal()

  sample_list <- reactiveVal()

  classified_samples_list <- reactiveVal()

  classified_list <- reactiveVal()

  cohort_analysis_list <- reactiveVal()

  cohort_sample_list <- reactiveVal()

  hist_list <- reactiveVal()

  qual_list <- reactiveVal()

  mean_df <- reactiveVal()

  processed_reads_df <- reactiveVal()

  taxa_count_table <- reactiveVal()

  plot_read_histogram <- reactiveVal()

  plot_q_histogram <- reactiveVal()

  plot_taxa_bar <- reactiveVal()

  rarefaction_combined <- reactiveVal(data.frame(
    Sample_Id = character(0),
    Classified_reads = numeric(0),
    Unique_taxa = numeric(0)
  ))

  diversity_combined <- reactiveVal(data.frame(
    Sample_Id = character(0),
    Classified_reads = numeric(0),
    Diversity = character(0),
    Value = numeric(0)
  ))

  plot_rarefaction_curve <- reactiveVal()

  plot_diversity_curve <- reactiveVal()

  result_dir_val <- reactiveVal()

  reactive_rel_abundance_matrix <- reactiveVal()

  reactive_counts_data <- reactiveVal()

  reactive_counts_matrix <- reactiveVal()

  realtime_daa_data <- reactiveVal()

  realtime_daa_metadata <- reactiveVal()

  realtime_daa_prev_cutoff <- reactiveVal()

  realtime_daa_counts_cutoff <- reactiveVal()

  realtime_taxa_group <- reactiveVal()

  abundance_val <- reactiveVal()
    
  plot_taxa_stacked <- reactiveVal()

  plot_diversity_box <- reactiveVal()

  plot_pcoa_dot <- reactiveVal()

  plot_nmds_dot <- reactiveVal()

  plot_pca_dot <- reactiveVal()

  table_permanova <- reactiveVal()

  plot_real_heatmap <- reactiveVal()
  
  plot_daa_volcano <- reactiveVal()

  table_ancombc_daa <- reactiveVal()

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

    trigger(1)
    
    route("Realtime")
    
    updateActionButton(session, "start_analysis", 
                       label = if (is_running()) "Stop Analysis" else "Start Analysis")
    
    if (is_running()) {
      
      status_checked(FALSE)
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

    req(trigger() > 0)
    
    req(is_running(), status_checked(), state() == "Sequencing")

    data_base <- input$realtime_database

    min_len <- input$realtime_min

    max_len <- input$realtime_max

    q_score <- input$realtime_q_score

    coverage <- input$realtime_cov

    identity <- input$realtime_iden

    threads <- input$realtime_threads

    conf_score <- input$realtime_conf
      
    taxa <- input$realtime_tax

    input_pipeline <- input$realtime_pipeline

    batch_size <- as.integer(input$chunk_size)

    if (realtime_task$status() != "running") {

      isolate({
        
        realtime_task$invoke(min_len, max_len, q_score, taxa, conf_score, coverage, identity, data_base, threads, input_pipeline, batch_size)
      
      })

    }

  })

  observeEvent(realtime_task$result(), {
    
    tryCatch({
      
      req(realtime_task$result(), input$taxon_counts_cutoff)
    
      result <- realtime_task$result()

      sample_list(result$sample_list)

      classified_samples_list(result$classified_samples_list)
      
      mean_df(result$mean_read_length_df)

      processed_reads_df(result$processed_reads_df)

      hist_list(result$hist_data_list)

      qual_list(result$quality_data_list)

      classified_list(result$classification_data_list)

      sample_names <- sample_list()$Barcode

      sample_id <- vector()

      classified_reads <- vector()

      unique_taxa <- vector()

      abundance_data_list <- list()

      lineage <-  input$taxon_select

      if(input$taxon_select=="Species") {
        
        for(i in seq_along(sample_names)) {
          
          sample_id[i] <- sample_names[i]
          
          abundance_data_list[[i]] <- classified_list()[[i]] %>%
                                      group_by(.data[[input$taxon_select]]) %>%
                                      summarise(Counts = sum(Counts), .groups = "drop") %>%
                                      filter(.data[[input$taxon_select]] != "Unclassified" & Counts >= input$taxon_counts_cutoff)

          unique_taxa[i] <- nrow(abundance_data_list[[i]])
          
          classified_reads[i] <- sum(abundance_data_list[[i]]$Counts)

          colnames(abundance_data_list[[i]]) <- c(input$taxon_select, sample_names[i])
        
        }

        current_snapshot <- data.frame(Sample_Id = sample_id, Classified_reads = classified_reads, Unique_taxa = unique_taxa)

        old_data <- rarefaction_combined()
        
        new_data <- rbind(old_data, current_snapshot)

        rarefaction_combined(new_data)

        old_diversity_data <- diversity_combined()

        counts_data <- abundance_data_list %>% purrr::reduce(full_join, by=lineage)

        counts_data[is.na(counts_data)] <- 0

        counts_data <- as.data.frame(counts_data)

        counts_data_long <- counts_data %>% pivot_longer(cols = -!!sym(lineage), names_to = "Sample", values_to = "Counts")

        sample_total_counts_df <- counts_data_long %>% group_by(Sample) %>% summarise(Total_Counts = sum(Counts))

        colnames(sample_total_counts_df) <- c("Sample_Id", "Counts")

        counts_data_wide <- counts_data_long %>% pivot_wider(names_from = !!sym(lineage), values_from = Counts, values_fill = 0) %>% 
          as.data.frame()
        
        counts_data_wide[is.na(counts_data_wide)] <- 0

        counts_data_wide <- counts_data_wide %>% dplyr::select(Sample, everything())

        rownames(counts_data_wide) <- counts_data_wide$Sample

        counts_data_wide <- counts_data_wide[,-1]

        counts_data_wide <- round(counts_data_wide, digits = 0)

        diversity_data <- counts_data_wide %>% 
                          as.data.frame() %>% rownames_to_column("Sample_Id") %>%
                          pivot_longer(cols = -Sample_Id, names_to = lineage, values_to = "Counts") %>% 
                          group_by(Sample_Id) %>% 
                          summarise(
                            Shannon = diversity(Counts, index = "shannon"),
                            Simpson = diversity(Counts, index = "simpson")
                          )

        diversity_data <- diversity_data %>% pivot_longer(cols = c(Shannon, Simpson),
                                              names_to = "Diversity", values_to = "Value")

        diversity_data <- inner_join(diversity_data, sample_total_counts_df, by = "Sample_Id")

        diversity_data <- diversity_data %>% dplyr::select(Sample_Id, Counts, Diversity, Value)

        colnames(diversity_data) <- c("Sample_Id", "Classified_reads", "Diversity", "Value")

        new_diveristy_data <- rbind(old_diversity_data, diversity_data)

        diversity_combined(new_diveristy_data)
      
      }

      refresh_time <- as.integer(input$refresh_rate)*1000

      delay(refresh_time, {
        trigger(trigger() + 1)
      })
    
    }, error = function(e) {
      
      showNotification(paste("Realtime Run Error: ", e$message), type = "error")

      refresh_time <- as.integer(input$refresh_rate)*1000

      delay(refresh_time, {
        trigger(trigger() + 1)
      })
    
    })
  
  })

  observeEvent(trigger(), {
    
    req(trigger() >= 6, cohort_delay_done() == FALSE)

    print(paste0("Cohort Delay Done at ", Sys.time(), " after ", trigger(), " triggers."))
    
    cohort_delay_done(TRUE)

    cohort_trigger(1)
  })

  observe({

    req(is_running(), cohort_delay_done(), status_checked(), state()=="Sequencing", classified_list(), classified_samples_list())
      
    cohort_analysis_list(classified_list())

    cohort_sample_list(classified_samples_list())

    invalidateLater(120000, session)
  
  })

  daa_reactive_conditions <- reactive({
    
    list(cohort_analysis_list(), cohort_sample_list(), input$taxa,
          input$prevalence_cutoff, input$abundance_cutoff, input$counts_cutoff)

  })

  observeEvent(ignoreInit = TRUE, daa_reactive_conditions(), {
      
      req(route()=="Realtime", is_running(), cohort_delay_done(), input$realtime_control, status_checked(), state()=="Sequencing", input_data_reactive())

      req(input$taxa, input$prevalence_cutoff, input$counts_cutoff, cohort_analysis_list(), cohort_sample_list(), input$abundance_cutoff)

      sample_list <- cohort_sample_list()$Barcode

      lineage <- input$taxa

      sample_metadata <- input_data_reactive()$data

      control_group <- input$realtime_control

      prevalence_cutoff <- input$prevalence_cutoff

      counts_cutoff <- input$counts_cutoff

      abundance_cutoff <- input$abundance_cutoff

      counts_matrix <- cohort_realtime_analysis(cohort_analysis_list(), sample_list, lineage, prevalence_cutoff, abundance_cutoff)$counts_matrix

      if (daa_ancombc_run$status() != "running") {

        isolate({
          
          daa_ancombc_run$invoke(counts_matrix, lineage, sample_metadata, control_group, prevalence_cutoff, counts_cutoff)
        
        })

      }

    }
  )

  observeEvent(daa_ancombc_run$result(), {

    tryCatch({

      req(daa_ancombc_run$result())

      req(nrow(daa_ancombc_run$result()$res_dunn)>0)

      result <- daa_ancombc_run$result()

      res_dunn <- result$res_dunn

      metadata <- result$metadata

      lineage <- result$lineage

      prev_cutoff <- result$prevalence_cutoff

      counts_cutoff <- result$counts_cutoff

      realtime_daa_data(res_dunn)

      realtime_daa_metadata(metadata)

      realtime_daa_counts_cutoff(counts_cutoff)

      realtime_daa_prev_cutoff(prev_cutoff)

      realtime_taxa_group(lineage)
    
    }, error = function(e) {

      showNotification(paste("Real-time ANCOM-BC2 Error: ", e$message), type = "error")

    })

  })

  observe({
    
    req(cohort_run(), route()=="Offline")

    work_dir <- getwd()

    install_dir <- paste0(work_dir, "/Installation")

    pipeline_dir <- paste0(work_dir, "/Pipelines")

    if(input$pipeline=="BLASTn")
    {
      if(input$setup)
      {
        print("Setting up BLASTn for 16S Data Analysis....")

        system(paste0('bash ', install_dir,'/blast_install.sh'))

        print("Analyzing Data....")

        system(paste0("bash ", pipeline_dir, "/blast_run.sh -p ", input$fastqdir, " -k ", input$kitname,
        " -t ", input$threads, " -m ", input$min, " -M ", input$max, " -i ", input$iden
        , " -c ", input$cov, " -q ", input$q_score, " -n ", input$database))

        result_dir <- paste0(input$fastqdir, "/Blast_Results/", input$database,"/")

        result_dir_val(result_dir)

      }
      else
      {
        print("Analyzing Data....")

        system(paste0("bash ", pipeline_dir, "/blast_run.sh -p ", input$fastqdir, " -k ", input$kitname,
                      " -t ", input$threads, " -m ", input$min, " -M ", input$max, " -i ", input$iden
                      , " -c ", input$cov, " -q ", input$q_score, " -n ", input$database))

        result_dir <- paste0(input$fastqdir, "/Blast_Results/", input$database,"/")

        result_dir_val(result_dir)
      
      }
    }
    
    else if(input$pipeline == "Kraken2") 
    {
      if(input$setup) 
      {
        print("Installation Already Completed.")
        
        print("Analyzing Data....")
        
        system(paste0("bash ", pipeline_dir, "/kraken_run.sh -p ", input$fastqdir,
                      " -k ", input$kitname, " -t ",
                      input$threads, " -m ", input$min, " -M ",
                      input$max, " -r ", input$tax, " -c ", input$conf
                      , " -q ", input$q_score, " -n ", input$database))
        result_dir <- paste0(input$fastqdir, "/Kraken2_Results/", input$database,"/")
        
        result_dir_val(result_dir)
      }
      else
      {
        print("Analyzing Data....")
        
        system(paste0("bash ", pipeline_dir, "/kraken_run.sh -p ", input$fastqdir,
                      " -k ", input$kitname, " -t ",
                      input$threads, " -m ", input$min, " -M ",
                      input$max, " -r ", input$tax, " -c ", input$conf
                      , " -q ", input$q_score, " -n ", input$database))
        result_dir <- paste0(input$fastqdir, "/Kraken2_Results/", input$database,"/")
        
        result_dir_val(result_dir)
      }
    }
    
    else if(input$pipeline == "Minimap2")
    {
      if(input$setup)
      {
        print("Installation Already Completed.")
        
        print("Analyzing Data....")
        
        pipeline_path <- paste0(getwd(), "/Pipelines")
        
        system(paste0("bash ", pipeline_dir, "/minimap2_run.sh -p ", input$fastqdir, " -k ", input$kitname,
                      " -s ", pipeline_path," -t ", input$threads, " -m ", input$min, " -M ", input$max, " -i ", input$iden
                      , " -c ", input$cov, " -q ", input$q_score, " -n ", input$database))
        result_dir <- paste0(input$fastqdir,"/Minimap2_Results/", input$database,"/")
        
        result_dir_val(result_dir)
      }
      else
      {
        print("Analyzing Data....")
        
        pipeline_path <- paste0(getwd(), "/Pipelines")
        
        system(paste0("bash ", pipeline_dir, "/minimap2_run.sh -p ", input$fastqdir, " -k ", input$kitname,
                      " -s ", pipeline_path," -t ", input$threads, " -m ", input$min, " -M ", input$max, " -i ", input$iden
                      , " -c ", input$cov, " -q ", input$q_score, " -n ", input$database))
        result_dir <- paste0(input$fastqdir,"/Minimap2_Results/", input$database,"/")
        
        result_dir_val(result_dir)
      }
    }

    else if(input$pipeline == "EMU")
    {
      if(input$setup)
      {
        print("Setting up EMU Pipeline for 16S Data Analysis....")

        system(paste0("bash ", install_dir, "/emu_install.sh"))

        print("Analyzing Data....")

        system(paste0("bash ", pipeline_dir, "/emu_run.sh -p ", input$fastqdir, " -k ", input$kitname, " -t ",
        input$threads, " -m ", input$min, " -M ",
        input$max, " -q ", input$q_score, " -n ", input$database))

        result_dir <- paste0(input$fastqdir,"/EMU_Results/", input$database,"/")

        result_dir_val(result_dir)

      }
      else
      {
        print("Analyzing Data....")

        system(paste0("bash ", pipeline_dir, "/emu_run.sh -p ", input$fastqdir, " -k ", input$kitname, " -t ",
        input$threads, " -m ", input$min, " -M ",
        input$max, " -q ", input$q_score, " -n ", input$database))
      
        result_dir <- paste0(input$fastqdir,"/EMU_Results/", input$database,"/")

        result_dir_val(result_dir)

      }
    }
  
  })

  observe({
    
    req(cohort_run(), route()=="Offline", result_dir_val(), input$taxa, input$prevalence_cutoff, input$abundance_cutoff)

    lineage <- input$taxa

    prevalence_cutoff <- input$prevalence_cutoff

    abundance_cutoff <- input$abundance_cutoff

    dir <- result_dir_val()

    result_list <- cohort_offline_analysis(dir, lineage, prevalence_cutoff, abundance_cutoff)

    reactive_rel_abundance_matrix(result_list$rel_abundance_renormalized_matrix)

    reactive_counts_data(result_list$counts_data)

    reactive_counts_matrix(result_list$counts_matrix)

  })

  observe({
    
    req(route()=="Example", input$taxa, input$prevalence_cutoff, input$abundance_cutoff)

    file_list <- gsub(".txt", "", list.files("Example/")[grep("\\_final_emu_result.txt$",
                                                                                                list.files("Example/"))])
      
    lineage <- input$taxa

    sample_metadata <- input_data_reactive()$data

    prevalence_cutoff <- input$prevalence_cutoff

    abundance_cutoff <- input$abundance_cutoff

    result_list <- example_analysis("Example", file_list, lineage, prevalence_cutoff, abundance_cutoff)

    reactive_rel_abundance_matrix(result_list$rel_abundance_renormalized_matrix)

    reactive_counts_data(result_list$counts_data)

    reactive_counts_matrix(result_list$counts_matrix)
  
  })

  observe({

    req(nrow(reactive_rel_abundance_matrix())>0)

    taxa_names_data <- sort(rownames(reactive_rel_abundance_matrix()))

    choices_list <- setNames(taxa_names_data, taxa_names_data)

    current_selection <- isolate(input$toi)

    updateSelectizeInput(session, 'toi', choices = choices_list, server = TRUE, selected = if (is.null(current_selection)) NULL else intersect(current_selection, taxa_names_data))

  })

  check_sequencing_status <- function() {
    
    path <- paste0(getwd(),"/Scripts")
      
    sequencing_state_check <- paste0(path,"/get_sequencing_state.py")
    
    state(py_run_file(sequencing_state_check)$status)
  }
  
  run_analysis_scripts <- function() {
    
    script_path <- paste0(getwd(),"/Scripts")
    
    pore_check <- paste0(script_path,"/get_active_pores.py")
    
    pores(py_run_file(pore_check)$pore_counts)

  }

  realtime_task <- ExtendedTask$new(function(min, max, q_score, taxa, conf, coverage, identity, data_base, threads, pipeline, chunk_size) {
    
    future_promise({
      
      script_path <- paste0(getwd(), "/Scripts")

      pipeline_path <- paste0(getwd(), "/Pipelines")
      
      output_path_check <- paste0(script_path, "/get_minknow_output_dir.py")

      reads_path <- paste0(py_run_file(output_path_check)$data_directory, "/fastq_pass")

      kit_name <- py_run_file(paste0(script_path, "/get_kit_name.py"))$kit_name

      if(!dir.exists(reads_path)) {
        
        return(NULL)
      
      } else {

        run_log <- system(paste0("bash ",pipeline_path,"/main_run.sh -d ", reads_path," -k ", kit_name, " -s ", pipeline_path, " -m ", min, " -M ", max,
                      " -t ", threads, " -v ", conf, " -q ", q_score, " -n ", data_base, " -i ", identity, 
                      " -c ", coverage, " -r ", taxa, " -p ", pipeline, " -b ", chunk_size," 2>&1 | tee -a ", reads_path, "/realtime_", pipeline,"_run.log"), intern = TRUE)

        length_list <- gsub("/.*$", "", list.files(reads_path, pattern = "average_length.txt", recursive = TRUE))

        run_status <- run_log[1]

        if(run_status == "No new files found.") {
          
          return(NULL)
        
        } else if (run_status == "New files found, now processing.") {

          sample_list <- data.frame(Barcode = length_list)
        
          mean_read_length_df <- data.frame()

          processed_reads_df <- data.frame()

          quality_data_list <- list()

          hist_data_list <- list()

          for(i in 1:length(length_list))
          {
            
            tmp <- read.table(file = paste0(reads_path, "/", length_list[i], "/", data_base, "/", length_list[i], "_average_length.txt"), sep = "\t", header = FALSE)

            mean_read_length_df <- rbind(mean_read_length_df, tmp)

            tmp <- read.table(file = paste0(reads_path, "/", length_list[i], "/", data_base, "/", length_list[i], "_processed_reads.txt"), sep = "\t", header = FALSE)

            processed_reads_df <- rbind(processed_reads_df, tmp)

            quality_data_list[[i]] <- read.table(file = paste0(reads_path, "/", length_list[i], "/", data_base, "/", length_list[i], "_quality.txt"), sep = "\t", header = FALSE)

            colnames(quality_data_list[[i]]) <- "Phred"

            hist_data_list[[i]] <- read.table(file = paste0(reads_path, "/", length_list[i], "/", data_base, "/", length_list[i], "_hist.txt"), sep = "\t", header = FALSE)

            colnames(hist_data_list[[i]]) <- c("Bin", "Counts")
          
          }
          
          colnames(mean_read_length_df) <- c("Barcode", "Mean_Length")

          colnames(processed_reads_df) <- c("Barcode", "Processed_Reads")

          classification_data_list <- list()

          if(pipeline=="Kraken2") {
            
            classified_samples <- gsub(paste0("^.*", data_base, "/"), "", Sys.glob(file.path(reads_path, "barcode*", data_base, "*final_kraken2_result.txt"))) %>%
                                  gsub("_final_kraken2_result.txt", "", x = .)

          } else if(pipeline=="Minimap2") {
            
            classified_samples <- gsub(paste0("^.*", data_base, "/"), "", Sys.glob(file.path(reads_path, "barcode*", data_base, "*final_minimap2_result.txt"))) %>%
                                  gsub("_final_minimap2_result.txt", "", x = .)

          } else if(pipeline=="BLASTn") {
            
            classified_samples <- gsub(paste0("^.*", data_base, "/"), "", Sys.glob(file.path(reads_path, "barcode*", data_base, "*final_blast_result.txt"))) %>%
                                  gsub("_final_blast_result.txt", "", x = .)

          } else if(pipeline=="EMU") {
            
            classified_samples <- gsub(paste0("^.*", data_base, "/"), "", Sys.glob(file.path(reads_path, "barcode*", data_base, "*final_emu_result.txt"))) %>%
                                  gsub("_final_emu_result.txt", "", x = .)

          }

          classified_samples_list <- data.frame(Barcode = classified_samples)

          for (i in 1:length(classified_samples))
          {
              
              if(pipeline=="Kraken2") {
                
                classification_data_list[[i]] <- read.delim(file = paste0(reads_path, "/", classified_samples[i], "/", 
                data_base, "/", classified_samples[i], "_final_kraken2_result.txt"), header = FALSE, sep = "\t")
              
              } else if(pipeline=="Minimap2") {
                
                classification_data_list[[i]] <- read.delim(file = paste0(reads_path, "/", classified_samples[i], "/", 
                data_base, "/", classified_samples[i], "_final_minimap2_result.txt"), header = FALSE, sep = "\t")
              
              } else if(pipeline=="BLASTn") {
                
                classification_data_list[[i]] <- read.delim(file = paste0(reads_path, "/", classified_samples[i], "/", 
                data_base, "/", classified_samples[i], "_final_blast_result.txt"), header = FALSE, sep = "\t")
              
              } else if(pipeline=="EMU") {
                
                classification_data_list[[i]] <- read.delim(file = paste0(reads_path, "/", classified_samples[i], "/", 
                data_base, "/", classified_samples[i], "_final_emu_result.txt"), header = FALSE, sep = "\t")
              
              }
              
              colnames(classification_data_list[[i]]) <- c("TAX_ID", "Counts", "Kingdom", "Phylum",
                                                  "Class", "Order", "Family", "Genus", "Species")
              
              classification_data_list[[i]][classification_data_list[[i]]==""] <- "Unclassified"
          }

          list(
            sample_list = sample_list,
            mean_read_length_df = mean_read_length_df,
            processed_reads_df = processed_reads_df,
            hist_data_list = hist_data_list,
            quality_data_list = quality_data_list,
            classified_samples_list = classified_samples_list,
            classification_data_list = classification_data_list
          )

        } else {

          stop("Realtime Task failed with status: ", run_log)

        }
      }
    }, seed = TRUE)
  })

  example_analysis <- function(path, file_list, lineage, prevalence_cutoff, abundance_cutoff) 
  {
    samples_header <- gsub("_final_emu_result","",file_list)
    
    sample_data_list <- list()
    
    abundance_data_list <- list()
    
    rel_abundance_data_list <- list()
    
    for(i in 1:length(file_list))
    {
      sample_data_list[[i]] <- read.delim(paste0(path,"/",file_list[i],".txt"), header = FALSE)

      sample_data_list[[i]] <- sample_data_list[[i]][,1:9]

      colnames(sample_data_list[[i]]) <- c("TAX_ID", "Counts", "Kingdom", "Phylum",
                                            "Class", "Order", "Family", "Genus",
                                            "Species")

      sample_data_list[[i]][sample_data_list[[i]]==""] <- "Unclassified"

      col_num <- which(colnames(sample_data_list[[i]])==lineage)

      abundance_data_list[[i]] <- sample_data_list[[i]] %>% group_by(.data[[lineage]]) %>% 
        summarise(Counts = sum(Counts))

      abundance_data_list[[i]] <- abundance_data_list[[i]] %>% filter(.data[[lineage]] != "Unclassified")

      rel_abundance_data_list[[i]] <- abundance_data_list[[i]] %>% mutate(Freq = (Counts/sum(Counts)))

      rel_abundance_data_list[[i]] <- rel_abundance_data_list[[i]][,-2]

      colnames(rel_abundance_data_list[[i]]) <- c(lineage, samples_header[i])

      colnames(abundance_data_list[[i]]) <- c(lineage, samples_header[i])

    }

    counts_data <- abundance_data_list %>% purrr::reduce(full_join, by=lineage)

    counts_data[is.na(counts_data)] <- 0

    counts_data <- as.data.frame(counts_data)

    counts_matrix <- as.matrix(counts_data[,-1])

    rownames(counts_matrix) <- counts_data[, which(colnames(counts_data)==lineage)]

    rel_abundance_data <- rel_abundance_data_list %>% purrr::reduce(full_join, by=lineage)

    rel_abundance_data <- as.data.frame(rel_abundance_data)

    rel_abundance_data[is.na(rel_abundance_data)] <- 0

    rel_abundance_matrix <- as.matrix(rel_abundance_data[,-1])

    rownames(rel_abundance_matrix) <- rel_abundance_data[,which(colnames(rel_abundance_data)==lineage)]

    pre_cutoff <- round((length(samples_header)*prevalence_cutoff)/100,0)

    mean_cutoff <- abundance_cutoff/100

    rel_abundance_filtered_matrix <- rel_abundance_matrix[rowSums(rel_abundance_matrix>0)>=pre_cutoff & rowMeans(rel_abundance_matrix)>=mean_cutoff,]

    rel_abundance_filtered_matrix[rel_abundance_filtered_matrix==0] <- 1e-10
    
    rel_abundance_renormalized_matrix <- apply(rel_abundance_filtered_matrix, 2 , function(x) x/sum(x))

    return(list(
    'rel_abundance_renormalized_matrix' = rel_abundance_renormalized_matrix,
    'counts_matrix' = counts_matrix,
    'counts_data' = counts_data))
  }
  
  cohort_realtime_analysis <- function(classification_list, sample_list, lineage, prevalence_cutoff, abundance_cutoff) {
      
      data <- classification_list

      abundance_data_list <- list()

      rel_abundance_data_list <- list()

      for(i in 1:length(sample_list))
        {

          col_num <- which(colnames(data[[i]])==lineage)
          
          abundance_data_list[[i]] <- data[[i]] %>% group_by(.data[[lineage]]) %>% 
          summarise(Counts = sum(Counts))

          abundance_data_list[[i]] <- abundance_data_list[[i]] %>% filter(.data[[lineage]] != "Unclassified")

          rel_abundance_data_list[[i]] <- abundance_data_list[[i]] %>% mutate(Freq = (Counts/sum(Counts)))

          rel_abundance_data_list[[i]] <- rel_abundance_data_list[[i]][,-2]

          colnames(rel_abundance_data_list[[i]]) <- c(lineage, sample_list[i])

          colnames(abundance_data_list[[i]]) <- c(lineage, sample_list[i])
        
      }
      
      counts_data <- abundance_data_list %>% purrr::reduce(full_join, by=lineage)

      counts_data[is.na(counts_data)] <- 0

      counts_data <- as.data.frame(counts_data)

      counts_matrix <- as.matrix(counts_data[,-1])

      rownames(counts_matrix) <- counts_data[, which(colnames(counts_data)==lineage)]

      rel_abundance_data <- rel_abundance_data_list %>% purrr::reduce(full_join, by=lineage)

      rel_abundance_data <- as.data.frame(rel_abundance_data)

      rel_abundance_data[is.na(rel_abundance_data)] <- 0

      rel_abundance_matrix <- as.matrix(rel_abundance_data[,-1])

      rownames(rel_abundance_matrix) <- rel_abundance_data[,which(colnames(rel_abundance_data)==lineage)]

      pre_cutoff <- round((length(sample_list)*prevalence_cutoff)/100,0)

      mean_cutoff <- abundance_cutoff/100

      rel_abundance_filtered_matrix <- rel_abundance_matrix[rowSums(rel_abundance_matrix>0)>=pre_cutoff & rowMeans(rel_abundance_matrix)>=mean_cutoff,]

      rel_abundance_filtered_matrix[rel_abundance_filtered_matrix==0] <- 1e-10
      
      rel_abundance_renormalized_matrix <- apply(rel_abundance_filtered_matrix, 2 , function(x) x/sum(x))

      return(list(
        'rel_abundance_renormalized_matrix' = rel_abundance_renormalized_matrix,
        'counts_data' = counts_data,
        'counts_matrix' = counts_matrix))
  }

  cohort_offline_analysis <- function(result_dir, lineage, prevalence_cutoff, abundance_cutoff) {
    
    file_list <- gsub(".txt", "", list.files(result_dir))[grep("\\_final.*result.txt", list.files(result_dir))]

    samples_header <- gsub("_final.*result", "", file_list)

    sample_data_list <- list()

    abundance_data_list <- list()

    rel_abundance_data_list <- list()

    for (i in 1:length(file_list))
    {
      sample_data_list[[i]] <- read.delim(file = paste0(result_dir, "/", file_list[i], ".txt"), header = FALSE)

      sample_data_list[[i]] <- sample_data_list[[i]][,1:9]

      colnames(sample_data_list[[i]]) <- c("TAX_ID", "Counts", "Kingdom", "Phylum",
                                                                "Class", "Order", "Family", "Genus", "Species")
                            
      sample_data_list[[i]][sample_data_list[[i]]==""] <- "Unclassified"
      
      col_num <- which(colnames(sample_data_list[[i]])==lineage)

      abundance_data_list[[i]] <- sample_data_list[[i]] %>% group_by(.data[[lineage]]) %>% 
        summarise(Counts = sum(Counts))

      abundance_data_list[[i]] <- abundance_data_list[[i]] %>% filter(.data[[lineage]] != "Unclassified")

      rel_abundance_data_list[[i]] <- abundance_data_list[[i]] %>% mutate(Freq = (Counts/sum(Counts)))

      rel_abundance_data_list[[i]] <- rel_abundance_data_list[[i]][,-2]

      colnames(rel_abundance_data_list[[i]]) <- c(lineage, samples_header[i])

      colnames(abundance_data_list[[i]]) <- c(lineage, samples_header[i])
    
    }

    counts_data <- abundance_data_list %>% purrr::reduce(full_join, by=lineage)

    counts_data[is.na(counts_data)] <- 0

    counts_data <- as.data.frame(counts_data)

    counts_matrix <- as.matrix(counts_data[,-1])

    rownames(counts_matrix) <- counts_data[, which(colnames(counts_data)==lineage)]

    rel_abundance_data <- rel_abundance_data_list %>% purrr::reduce(full_join, by=lineage)

    rel_abundance_data <- as.data.frame(rel_abundance_data)

    rel_abundance_data[is.na(rel_abundance_data)] <- 0

    rel_abundance_matrix <- as.matrix(rel_abundance_data[,-1])

    rownames(rel_abundance_matrix) <- rel_abundance_data[,which(colnames(rel_abundance_data)==lineage)]

    pre_cutoff <- round((length(samples_header)*prevalence_cutoff)/100,0)

    mean_cutoff <- abundance_cutoff/100

    rel_abundance_filtered_matrix <- rel_abundance_matrix[rowSums(rel_abundance_matrix>0)>=pre_cutoff & rowMeans(rel_abundance_matrix)>=mean_cutoff, , drop = FALSE]

    if(nrow(rel_abundance_filtered_matrix)>1) {
      
      rel_abundance_filtered_matrix <- rel_abundance_filtered_matrix + 1e-10

      rel_abundance_renormalized_matrix <- apply(rel_abundance_filtered_matrix, 2, function(x) x/sum(x))
    
    } else {

      rel_abundance_renormalized_matrix <- rel_abundance_filtered_matrix

    }
    
    return(list(
      'rel_abundance_renormalized_matrix' = rel_abundance_renormalized_matrix,
      'counts_data' = counts_data,
      'counts_matrix' = counts_matrix))
  }

  stacked_barplot_function <- function(stacked_df, lineage, top_n)
  {
    stacked_plot_data <- as.data.frame(stacked_df)
    
    stacked_plot_data[[lineage]] <- rownames(stacked_df)

    rownames(stacked_plot_data) <- NULL

    stacked_plot_data <- stacked_plot_data %>% dplyr::select(c(all_of(lineage), everything()))

    required_col <- which(colnames(stacked_plot_data)==lineage)

    stacked_plot_data[,required_col] <- as.character(stacked_plot_data[,required_col])

    stacked_plot_data[,-1] <- stacked_plot_data[,-1]*100

    stacked_plot_data <- stacked_plot_data %>% 
      pivot_longer(cols = -all_of(required_col), names_to = "Sample_Id", values_to = "Abundance") %>% 
      filter(Abundance > 0)
    
    stacked_plot_data <- stacked_plot_data %>% group_by(Sample_Id) %>% 
      arrange(desc(Abundance)) %>% dplyr::slice(1:ifelse(n()<top_n,n(),top_n))

    stacked_plot_data$Sample_Id <- gsub("barcode", "", stacked_plot_data$Sample_Id)

    temp <- stacked_plot_data %>% group_by(Sample_Id) %>% summarise(Abundance = 100-sum(Abundance)) %>% mutate(Species="Others") %>% 
      dplyr::select(Species, Sample_Id, Abundance)

    colnames(temp) <- colnames(stacked_plot_data)

    stacked_plot_data <- rbind(stacked_plot_data, temp) %>% arrange(Sample_Id)

    total_abundance <- stacked_plot_data %>%
      group_by(!!sym(lineage)) %>%
      summarise(Total_Abundance = sum(Abundance)) %>%
      arrange(Total_Abundance)

    stacked_plot_data[[required_col]] <- factor(stacked_plot_data[[required_col]], 
                                levels = c("Others", total_abundance[[required_col]][total_abundance[[required_col]] != "Others"]))

    fill_colors <- c("Others" = "#D3D3D3",setNames(viridis_pal(option = "D")(
      length(levels(stacked_plot_data[[required_col]]))-1), 
      levels(stacked_plot_data[[required_col]])[levels(stacked_plot_data[[required_col]]) != "Others"]))
    
    stacked_barplot <- ggplot(stacked_plot_data, aes(x = Sample_Id, y = Abundance, fill = !!sym(lineage))) +
    geom_bar(stat = "identity", position = "stack", color = "black", linewidth=0.2) +
    theme_linedraw() +
    theme(
      axis.text.x = element_text(size = 10, face = "bold", colour = "#07446a"),
      axis.text.y = element_text(size = 12, face = "bold", colour = "#07446a"),
      axis.title.x = element_text(size = 12, face = "bold", colour = "#07446a"),
      axis.title.y = element_text(size = 12, face = "bold", colour = "#07446a"),
      legend.title = element_text(size = 12, face = "bold", colour = "#07446a"),
      legend.text = element_text(face = "bold"),
      legend.position = "right",
      title = element_text(size = 14, face = "bold", colour = "#07446a"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing=unit(0,"lines"),
      strip.background=element_rect(color="grey30", fill="grey90"),
      panel.border=element_rect(color="grey90")
    ) +
    guides(fill = guide_legend(ncol = 4)) +
    labs(x = "BARCODE", y="RELATIVE ABUNDANCE (%)", fill = lineage) + 
    scale_fill_manual(values = fill_colors) +
    scale_y_continuous(expand = c(0, 0), limits = c(0,102))

    return(stacked_barplot)
  }

  stacked_subsampled_barplot_function <- function(stacked_df, lineage, species_of_interest) {

    stacked_plot_data <- as.data.frame(stacked_df)
    
    stacked_plot_data[[lineage]] <- rownames(stacked_df)

    rownames(stacked_plot_data) <- NULL

    stacked_plot_data <- stacked_plot_data %>% dplyr::select(c(all_of(lineage), everything()))

    required_col <- which(colnames(stacked_plot_data)==lineage)

    stacked_plot_data[,required_col] <- as.character(stacked_plot_data[,required_col])

    stacked_plot_data[,-1] <- stacked_plot_data[,-1]*100

    stacked_plot_data <- stacked_plot_data %>% 
      pivot_longer(cols = -all_of(required_col), names_to = "Sample_Id", values_to = "Abundance") %>% 
      filter(Abundance > 0)

    stacked_plot_data <- stacked_plot_data %>% group_by(Sample_Id) %>% 
      dplyr::filter(!!sym(lineage) %in% species_of_interest)

    stacked_plot_data$Sample_Id <- gsub("barcode", "", stacked_plot_data$Sample_Id)

    temp <- stacked_plot_data %>% group_by(Sample_Id) %>% summarise(Abundance = 100-sum(Abundance)) %>% mutate(Species="Others") %>% 
      dplyr::select(Species, Sample_Id, Abundance)

    colnames(temp) <- colnames(stacked_plot_data)

    stacked_plot_data <- rbind(stacked_plot_data, temp) %>% arrange(Sample_Id)

    total_abundance <- stacked_plot_data %>%
      group_by(!!sym(lineage)) %>%
      summarise(Total_Abundance = sum(Abundance)) %>%
      arrange(Total_Abundance)

    stacked_plot_data[[required_col]] <- factor(stacked_plot_data[[required_col]], 
                                levels = c("Others", total_abundance[[required_col]][total_abundance[[required_col]] != "Others"]))

    fill_colors <- c("Others" = "#D3D3D3",setNames(viridis_pal(option = "D")(
      length(levels(stacked_plot_data[[required_col]]))-1), 
      levels(stacked_plot_data[[required_col]])[levels(stacked_plot_data[[required_col]]) != "Others"]))
    
    stacked_barplot <- ggplot(stacked_plot_data, aes(x = Sample_Id, y = Abundance, fill = !!sym(lineage))) +
    geom_bar(stat = "identity", position = "stack", color = "black", linewidth=0.2) +
    theme_linedraw() +
    theme(
      axis.text.x = element_text(size = 10, face = "bold", colour = "#07446a"),
      axis.text.y = element_text(size = 12, face = "bold", colour = "#07446a"),
      axis.title.x = element_text(size = 12, face = "bold", colour = "#07446a"),
      axis.title.y = element_text(size = 12, face = "bold", colour = "#07446a"),
      legend.title = element_text(size = 12, face = "bold", colour = "#07446a"),
      legend.text = element_text(face = "bold"),
      legend.position = "right",
      title = element_text(size = 14, face = "bold", colour = "#07446a"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    guides(fill = guide_legend(ncol = 4)) +
    labs(x = "BARCODE", y="RELATIVE ABUNDANCE (%)", fill = lineage) + 
    scale_fill_manual(values = fill_colors) +
    scale_y_continuous(expand = c(0, 0), limits = c(0,102))

    return(stacked_barplot)
  
  }

  diversity_boxplot_function <- function(counts_data_wide, lineage, sample_metadata, rarefaction_cutoff)
  {
    
    alpha_diversity_data <- vegan::rrarefy(counts_data_wide, rarefaction_cutoff) %>% 
      as.data.frame() %>% rownames_to_column("Sample_Id") %>% 
      pivot_longer(cols = -Sample_Id, names_to = lineage, values_to = "Counts") %>% 
      group_by(Sample_Id) %>% 
      summarise(
        Shannon = diversity(Counts, index = "shannon"),
        Simpson = diversity(Counts, index = "simpson")
      )
    
    alpha_diversity_data <- inner_join(alpha_diversity_data, sample_metadata)
    
    alpha_diversity_data$Group <- factor(alpha_diversity_data$Group)
    
    alpha_diversity_data <- pivot_longer(alpha_diversity_data, cols = c(Shannon, Simpson), 
                                        names_to = "Diversity", values_to = "Value")
    
    alpha_div_p <- compare_means(Value~Group, data = alpha_diversity_data, method = "wilcox",
                                p.adjust.method = "BH", group.by = "Diversity")

    get_legend <- function(myggplot){
      tmp <- ggplot_gtable(ggplot_build(myggplot))
      leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
      legend <- tmp$grobs[[leg]]
      return(legend)
    }
    
    shannon_plot <- alpha_diversity_data %>% filter(Diversity == "Shannon") %>%
      ggplot(aes(x=Group, y=Value, color=Group)) +
      geom_boxplot() +
      geom_jitter(shape = 16, position = position_jitter(0.2)) +
      scale_color_manual(values = pal_aaas("default")(length(levels(alpha_diversity_data$Group)))) +
      theme_linedraw() +
      labs(y= "Alpha Diversity", x = "") +
      stat_pvalue_manual(subset(alpha_div_p, Diversity=="Shannon"), label = "p.signif", y.position = max(
        subset(alpha_diversity_data, Diversity=="Shannon")$Value) + 0.1, hide.ns = "p.adj", step.increase = 0.1,
        tip.length = 0.02, bracket.size = 0.8, size = 8, color = "#5B5DC7") +
      guides(color = guide_legend(title = "Shannon", title.position = "top")) +
      theme(
        axis.title.x = element_text(size = 14, face = "bold", colour = "#5B5DC7"),
        axis.title.y = element_text(size = 14, face = "bold", colour = "#5B5DC7", vjust=+2),
        strip.text.x = element_text(size = 14, face = "bold", colour = "#5B5DC7"),
        axis.text.y=element_text(size=14, face = "bold", colour = "#5B5DC7"),
        axis.text.x= element_blank(),
        axis.ticks.x = element_blank(),
        legend.title=element_blank(),
        legend.text=element_text(colour="#5B5DC7", size=12, face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 18, face = "bold", colour = "#5B5DC7", hjust = 0.5)
      ) +
      ggtitle("Shannon")
    
    simpson_plot <- alpha_diversity_data %>% filter(Diversity == "Simpson") %>%
      ggplot(aes(x=Group, y=Value, color=Group)) +
      geom_boxplot() +
      geom_jitter(shape = 16, position = position_jitter(0.2)) +
      scale_color_manual(values = pal_aaas("default")(length(levels(alpha_diversity_data$Group)))) + 
      theme_linedraw() +
      labs(y= "Alpha diversity", x = "") +
      stat_pvalue_manual(subset(alpha_div_p, Diversity=="Simpson"), y.position = max(
        subset(alpha_diversity_data, Diversity=="Simpson")$Value) + 0.01, label = "p.signif",
        hide.ns = "p.adj", step.increase = 0.1, tip.length = 0.02, bracket.size = 0.8, size = 8, color = "#5B5DC7") +
      guides(color = guide_legend(title = "Simpson", title.position = "top")) +
      theme(
        axis.title.x = element_text(size = 14, face = "bold", colour = "#5B5DC7"),
        axis.title.y = element_text(size = 14, face = "bold", colour = "#5B5DC7", vjust=+2),
        strip.text.x = element_text(size = 14, face = "bold", colour = "#5B5DC7"),
        axis.text.y=element_text(size=14, face = "bold", colour = "#5B5DC7"),
        axis.text.x= element_blank(),
        axis.ticks.x = element_blank(),
        legend.position="none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 18, face = "bold", colour = "#5B5DC7", hjust = 0.5)
      ) +
      ggtitle("Simpson")
    
    legend <- get_legend(shannon_plot)

    shannon_plot <- shannon_plot + theme(legend.position="none")
    
    diversity_facet <- as_ggplot(grid.grabExpr(grid.arrange(shannon_plot,
                                        simpson_plot, legend, ncol=3, widths=c(2.2, 2.2, 1.0)))) +
      labs(caption = paste0("Alpha Diversity metrices calculated on rarified ", lineage, " data with the 
                              cutoff library size of ", scales::comma(rarefaction_cutoff), " reads.<br>The P-value 
                              is calculated using Kruskal-Walis Test with Benjamini-Hochberg Correction.")) +
                                          theme(plot.caption = element_markdown(
                                            color = "#0F6E73", size = 15,
                                            margin = margin(20,0,5,0), face = "bold",
                                            hjust = 0.5
                                          ))
    return(diversity_facet)
  }

  diversity_pcoa_function <- function(mat, lineage, sample_metadata, prevalence_cutoff, abundance_cutoff)
  {
    
    pcoa_data <- mat %>% as.matrix()

    transformed_data <- t(compositions::clr(t(pcoa_data)))

    pcoa_dist <- wcmdscale(vegdist(t(transformed_data), method = "euclidean"), k=2, eig = TRUE)

    pcoa_df <- pcoa_dist$points[,1:2] %>% as.data.frame()
  
    pcoa_eigenvalues <- pcoa_dist$eig

    pcoa.var <- round(pcoa_eigenvalues/sum(pcoa_eigenvalues)*100, 1)

    pcoa_df$Sample_Id <- rownames(pcoa_df)

    pcoa_df <- inner_join(pcoa_df, sample_metadata)

    pcoa_df$Group <- factor(pcoa_df$Group)
  
    colnames(pcoa_df) <- c("Axis.1", "Axis.2", "Sample_Id", "Group")

    pcoa_plot <- ggplot(data = pcoa_df, aes(x = Axis.1, y = Axis.2, fill = Group)) +
      geom_point(size=5, shape = 21, color = "black") +
      stat_ellipse(aes(colour = Group, fill = Group), level = 0.95, alpha = 0.25, geom = "polygon") +
      scale_color_manual(values = pal_aaas("default")(length(levels(pcoa_df$Group)))) +
      scale_fill_manual(values = pal_aaas("default")(length(levels(pcoa_df$Group)))) +
      theme_linedraw() +
      xlab(paste0("PC1 (", pcoa.var[1], "%", ")")) +
      ylab(paste0("PC2 (", pcoa.var[2], "%", ")")) +
      geom_vline(xintercept = 0, color = "black", linetype = "dashed") +
      geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
      ggtitle("Principal Coordination Analysis (PCoA) ordination plot using Aitchison Distance") +
      labs(caption = paste0("The TSS Normalized Data was filtered to keep the ", lineage, " with ", prevalence_cutoff, 
                            "% prevalence and mean relative abundance >= ", abundance_cutoff,"%.<br>CLR transformation
                            was applied on the filtered data and Euclidean Distance was calculated and plotted.")) +
      theme(
        axis.text.x = element_text(size = 15, face = "bold", colour = "#5B5DC7"),
        axis.text.y = element_text(size = 15, face = "bold", colour = "#5B5DC7"),
        legend.text = element_text(size = 15, face = "bold", colour = "#5B5DC7"),
        axis.title.x = element_text(size = 15, face = "bold", colour = "#5B5DC7"),
        axis.title.y = element_text(size = 15, face = "bold", colour = "#5B5DC7"),
        legend.title = element_text(size = 15, face = "bold", colour = "#5B5DC7"),
        plot.title = element_text(size = 15, face = "bold", colour = "#5B5DC7", hjust = 0.5),
        plot.caption = element_markdown(
          color = "#0F6E73", size = 15,
          margin = margin(20, 0, 10, 0), face = "bold",
          hjust = 0.5
        )
      )
    
    return(pcoa_plot)
  }

  diversity_nmds_function <- function(mat, lineage, sample_metadata, prevalence_cutoff, abundance_cutoff)
  {
    nmds_dist <- metaMDS(t(mat), distance = "bray", trymax = 1000)

    nmds_df <- as.data.frame(nmds_dist$points)

    passed_df <- data.frame(Sample_Id = rownames(nmds_df))

    passed_df <- inner_join(passed_df, sample_metadata, by = "Sample_Id")

    nmds_df$Group <- passed_df$Group

    nmds_df$Group <- factor(nmds_df$Group)

    nmds_plot <- ggplot(data = nmds_df, aes(x = MDS1, y = MDS2, fill = Group)) +
      geom_point(size=5, shape = 21, color = "black") +
      stat_ellipse(aes(colour = Group, fill = Group), level = 0.95, alpha = 0.25, geom = "polygon") +
      scale_color_manual(values = pal_aaas("default")(length(levels(nmds_df$Group)))) +
      scale_fill_manual(values = pal_aaas("default")(length(levels(nmds_df$Group)))) +
      theme_linedraw() +
      ggtitle("Non-metric MultiDimensional Scaling (NMDS) ordination plot with Bray-Curtis Distance") +
      labs(caption = paste0("The TSS Normalized Data was filtered to keep ", lineage, " with ", prevalence_cutoff,
                            "% prevalence and mean relative abundance >= ", abundance_cutoff,"%.<br>Bray-Curtis
                            Distance was calculated on the filtered data and plotted.")) +
      theme(
        axis.text.x = element_text(size = 15, face = "bold", colour = "#5B5DC7"),
        axis.text.y = element_text(size = 15, face = "bold", colour = "#5B5DC7"),
        legend.text = element_text(size = 15, face = "bold", colour = "#5B5DC7"),
        axis.title.x = element_text(size = 15, face = "bold", colour = "#5B5DC7"),
        axis.title.y = element_text(size = 15, face = "bold", colour = "#5B5DC7", vjust = -1),
        legend.title = element_text(size = 15, face = "bold", colour = "#5B5DC7"),
        plot.title = element_text(size = 15, face = "bold", colour = "#5B5DC7", hjust = 0.5),
        plot.caption = element_markdown(
          color = "#0F6E73", size = 15,
          margin = margin(20, 0, 10, 0), face = "bold",
          hjust = 0.5
        )
      ) +
      geom_vline(xintercept = 0, color = "black", linetype = "dashed") +
      geom_hline(yintercept = 0, color = "black", linetype = "dashed")

    return(nmds_plot)
  }

  diversity_pca_function <- function(mat, lineage, sample_metadata, prevalence_cutoff, abundance_cutoff, top_taxa)
  {
    
    pca_data <- mat %>% as.matrix()

    transformed_data <- t(clr(t(pca_data)))

    pca_out <- FactoMineR::PCA(transformed_data, scale.unit = TRUE, ncp = 5, graph = FALSE)
    
    pca_data <- pca_out$var$coord %>% as.data.frame() %>% rownames_to_column("Sample_Id")

    pca_data <- inner_join(pca_data, sample_metadata)

    indv_data <- pca_out$ind$coord %>% head(n=top_taxa)

    sample_metadata$Group <- factor(sample_metadata$Group)

    pca_plot <- ggplot(data = pca_data, aes(x = Dim.1, y = Dim.2, fill = Group)) +
      geom_point(size=5, shape = 21, color = "black") +
      stat_ellipse(aes(colour = Group, fill = Group), level = 0.95, alpha = 0.4, geom = "polygon") +
      xlab(paste0("PC1 (", round(pca_out$eig[1,2],1), "%", ")")) +
      ylab(paste0("PC2 (", round(pca_out$eig[2,2],1), "%", ")")) +
      geom_segment(data = indv_data, aes(x=0, y=0, xend = Dim.1, yend = Dim.2), inherit.aes = FALSE,
                  arrow = arrow(type = "closed", length = unit(0.05, "inches")), color = "#0773C1") +
      geom_text_repel(data = indv_data, aes(Dim.1, Dim.2), inherit.aes = FALSE, label = rownames(indv_data), color = "#0773C1", max.overlaps = top_taxa) +
      scale_color_manual(values = pal_aaas("default")(length(levels(sample_metadata$Group)))) +
      scale_fill_manual(values = pal_aaas("default")(length(levels(sample_metadata$Group)))) +
      theme_linedraw() +
      ggtitle("Principal Component Analysis (PCA) Bi-Plot using CLR Transformed Data") +
      labs(caption = paste0("The TSS Normalized Data was filtered to keep ", lineage, " taxon with ", prevalence_cutoff,
                            "% prevalence and mean relative abundance >= ", abundance_cutoff,"%.<br>CLR transformation
                            was applied on the filtered data and biplot was plotted, showing loadings of top ", top_taxa,
                            " taxon.")) +
      theme(
        axis.text.x = element_text(size = 15, face = "bold", colour = "#5B5DC7"),
        axis.text.y = element_text(size = 15, face = "bold", colour = "#5B5DC7"),
        legend.text = element_text(size = 15, face = "bold", colour = "#5B5DC7"),
        axis.title.x = element_text(size = 15, face = "bold", colour = "#5B5DC7"),
        axis.title.y = element_text(size = 15, face = "bold", colour = "#5B5DC7"),
        legend.title = element_text(size = 15, face = "bold", colour = "#5B5DC7"),
        plot.title = element_text(size = 15, face = "bold", colour = "#5B5DC7", hjust = 0.5),
        plot.caption = element_markdown(
          color = "#0F6E73", size = 15,
          margin = margin(20, 0, 10, 0), face = "bold",
          hjust = 0.5
        )
      ) +
      geom_vline(xintercept = 0, color = "black", linetype = "dashed") +
      geom_hline(yintercept = 0, color = "black", linetype = "dashed")

    return(pca_plot)
  }

  diversity_permanova_function <- function(mat, lineage, sample_metadata, control_group)
  {
    filtered_data <- mat %>% as.matrix()

    passed_df <- data.frame(Sample_Id = colnames(filtered_data))

    passed_df <- inner_join(passed_df, sample_metadata, by = "Sample_Id")

    transformed_data <- t(compositions::clr(t(filtered_data)))

    perm_dist <- vegdist(t(transformed_data), method = "euclidean")

    permanova_res <- pairwise.adonis(perm_dist, as.factor(passed_df$Group), p.adjust.m = "BH")

    permanova_res <- permanova_res %>% dplyr::select(c(pairs, R2, p.value, p.adjusted))

    colnames(permanova_res) <- c("Pair", "R2", "P", "Padj")

    permanova_res <- permanova_res[grepl(control_group, permanova_res$Pair), ]

    return(permanova_res)
  }

  diversity_heatmap_function <- function(mat, lineage, sample_metadata)
  {
    
    data <- mat %>% as.matrix()

    passed_df <- data.frame(Sample_Id = colnames(data))
    
    transformed_data <- t(clr(t(data)))

    heatmap_df <- as.data.frame(transformed_data)

    sample_metadata <- inner_join(passed_df, sample_metadata, by = "Sample_Id")

    sample_metadata$Group <- factor(sample_metadata$Group)

    sample_annotation <- data.frame(sample_metadata$Group)

    colnames(sample_annotation) <- "Group"

    sample_annotation$Group <- factor(sample_annotation$Group)

    col_list <- list(Group = setNames(pal_aaas("default")(length(levels(sample_metadata$Group))), 
                                 levels(sample_metadata$Group)))
    
    col_annotate <- columnAnnotation(
    df = sample_annotation,
    col = col_list,
    show_annotation_name = FALSE,
    annotation_legend_param = list(
      Group = list(
        title = "Groups",
        title_gp = gpar(fontsize = 12, fontface = "bold"),
        labels_gp = gpar(fontsize = 10),
        direction = "horizontal"
        )
      )
    )
    row_dend <- hclust(dist(heatmap_df, method = "euclidean"), method = "complete")

    column_dend <- hclust(dist(t(heatmap_df), method = "euclidean"), method = "complete")

    tmp <- heatmap_df %>% 
      pivot_longer(cols = everything(), 
              names_to = "Sample",
              values_to = "Abundance")
    
    heatmap_matrix <- as.matrix(heatmap_df)

    breaks <- seq(
      -ceiling(max(abs(heatmap_matrix))), ceiling(max(abs(heatmap_matrix))),
      length.out = ifelse(max(abs(heatmap_matrix))>5, 2*ceiling(max(abs(heatmap_matrix))), 10 ) 
    )
    
    col_fun <- colorRampPalette(c("darkblue", "blue", "white", "red", "darkred"))(length(breaks)-1)
    
    prevalence <- rowSums(data > 1e-9)
    
    heatmap_plot <- Heatmap(
      heatmap_matrix,
      name = "CLR",
      
      
      row_names_gp = gpar(fontsize = 8),
      row_names_max_width = unit(5, "cm"),
      row_names_side = "right",
      
      
      cluster_rows = color_branches(row_dend, k = 4),
      cluster_columns = color_branches(column_dend, k = 4),
      
      
      show_column_names = TRUE,
      column_names_gp = gpar(fontsize = 10, fontface = "bold"),
      column_title = NULL,
      column_title_gp = gpar(fontsize = 14, fontface = "bold"),
      column_dend_height = unit(15, "mm"),
    
      bottom_annotation = col_annotate,
    
      right_annotation = rowAnnotation(
        Prevalence = anno_barplot(prevalence,
        gp = gpar(fill = "#008B45FF", color = "white"),
        ylim = c(1,25)),
        annotation_name_rot = 0,
        annotation_name_side = "top",
        annotation_name_align = TRUE,
        annotation_name_gp = gpar(fontsize = 12, fontface = "bold"),
        border = TRUE,
        width = unit(2.5, "cm")
      ),
      
      
      col = col_fun,
      
      
      rect_gp = gpar(col = "white", lwd = 0.15),
      
      
      heatmap_legend_param = list(
        title_gp = gpar(fontsize = 12, fontface = "bold"),
        labels_gp = gpar(fontsize = 10),
        legend_height = unit(4, "cm"),
        grid_width = unit(0.5, "cm")
      ),
      
      
      clustering_distance_rows = "euclidean",
      clustering_distance_columns = "euclidean",
      clustering_method_rows = "complete",
      clustering_method_columns = "complete"
    )
    
    

    heatmap_ggplot <- as_ggplot(grid.grabExpr(print(heatmap_plot)))

    return(heatmap_ggplot)
  }

  daa_ancombc_run <- ExtendedTask$new(function(counts_matrix, lineage, sample_metadata, control_group, prevalence_cutoff, counts_cutoff) {

    future_promise({

      daa_matrix <- counts_matrix

      counts_cutoff <- counts_cutoff
      
      metadata <- sample_metadata

      rownames(metadata) <- metadata$Sample_Id

      control <- control_group

      metadata$Group <- as.factor(metadata$Group)

      metadata$Group <- relevel(metadata$Group, ref = control)

      prev_cutoff <- prevalence_cutoff/100

      passed_samples <- data.frame(Sample_Id = colnames(daa_matrix))

      passed_samples <- inner_join(passed_samples, sample_metadata, by = "Sample_Id")

      passed_samples <- passed_samples %>% group_by(Group) %>% summarise(Counts = n())

      if (nrow(passed_samples) > 0 && all(passed_samples$Counts >=2, na.rm=TRUE)) {
        
        output <- ancombc2(data = daa_matrix, meta_data = metadata,
          taxa_are_rows = TRUE,
          fix_formula = "Group", rand_formula = NULL,
          p_adj_method = "BH", pseudo_sens = TRUE,
          prv_cut = prev_cutoff, lib_cut = counts_cutoff, s0_perc = 0.05,
          group = "Group", struc_zero = TRUE, neg_lb = FALSE,
          alpha = 0.05, n_cl = 2, verbose = FALSE,
          global = FALSE, pairwise = FALSE, 
          dunnet = TRUE, trend = FALSE,
          iter_control = list(tol = 1e-5, max_iter = 20, 
          verbose = FALSE),
          em_control = list(tol = 1e-5, max_iter = 100),
          lme_control = lme4::lmerControl(), 
          mdfdr_control = list(fwer_ctrl_method = "BH", B = 100), 
          trend_control = NULL)

        res_dunn <- output$res_dunn

        return(list(
          'res_dunn' = res_dunn,
          'metadata' = metadata,
          'prevalence_cutoff' = prevalence_cutoff,
          'counts_cutoff' = counts_cutoff,
          'lineage' = lineage
        ))
      
      } else {
        
        return(NULL)
      
      }

    }, seed = TRUE)

  })

  daa_volcano_plot <- function(counts_matrix, lineage, sample_metadata, control_group, prevalence_cutoff, counts_cutoff) {
    
    daa_matrix <- counts_matrix

    passed_samples <- data.frame(Sample_Id = colnames(daa_matrix))

    passed_samples <- inner_join(passed_samples, sample_metadata, by = "Sample_Id")

    passed_samples <- passed_samples %>% group_by(Group) %>% summarise(Counts = n())

    if (nrow(passed_samples) > 0 && all(passed_samples$Counts >=2, na.rm=TRUE)) {
      
      counts_cutoff <- counts_cutoff
    
      metadata <- sample_metadata

      rownames(metadata) <- metadata$Sample_Id

      control <- control_group

      metadata$Group <- as.factor(metadata$Group)

      metadata$Group <- relevel(metadata$Group, ref = control)

      prev_cutoff <- prevalence_cutoff/100

      output <- ancombc2(data = daa_matrix, meta_data = metadata,
          taxa_are_rows = TRUE,
          fix_formula = "Group", rand_formula = NULL,
          p_adj_method = "BH", pseudo_sens = TRUE,
          prv_cut = prev_cutoff, lib_cut = counts_cutoff, s0_perc = 0.05,
          group = "Group", struc_zero = TRUE, neg_lb = FALSE,
          alpha = 0.05, n_cl = 2, verbose = FALSE,
          global = FALSE, pairwise = FALSE, 
          dunnet = TRUE, trend = FALSE,
          iter_control = list(tol = 1e-5, max_iter = 20, 
          verbose = FALSE),
          em_control = list(tol = 1e-5, max_iter = 100),
          lme_control = lme4::lmerControl(), 
          mdfdr_control = list(fwer_ctrl_method = "BH", B = 100), 
          trend_control = NULL)
      
      res_dunn <- output$res_dunn

      total_groups <- length(levels(metadata$Group))

      req_start_values <- c("passed_ss_Group", "diff_Group", "lfc_Group", "q_Group")

      out_col_values <- c("Sensitive", "Significance", "LFC", "P_adj")

      req_columns <- outer(levels(metadata$Group)[2:total_groups], req_start_values, function(x,y) paste0(y, x)) %>% 
        as.vector()

      filtered_dunn <- res_dunn[,c(1,which(colnames(res_dunn) %in% req_columns))]

      long_data_list <- list()

      for(i in 1:length(req_start_values)) {
        long_data_list[[i]] <- filtered_dunn %>% dplyr::select(c(taxon, starts_with(req_start_values[i]))) %>% 
          pivot_longer(cols = -taxon, names_to = "Comparison", values_to = out_col_values[i])

        long_data_list[[i]]$Comparison <- sub(req_start_values[i], "", long_data_list[[i]]$Comparison)
      }

      final_data <- purrr::reduce(long_data_list, left_join, by=c("taxon", "Comparison"))

      final_data <- final_data %>% dplyr::select(taxon, everything())

      colnames(final_data)[1] <- lineage

      final_data <- final_data %>% filter(Sensitive==TRUE)

      final_data$Comparison <- paste0(control, " - ", final_data$Comparison)

      final_data$Name <- ifelse(final_data$LFC < -1 & final_data$P_adj < 0.05, "Downregulated",
                                ifelse(final_data$LFC > 1 & final_data$P_adj < 0.05, "Upregulated", "Not Significant"))
      
      final_data$Name <- factor(final_data$Name, levels = c("Upregulated", "Downregulated", "Not Significant"))

      color_group <- c("Upregulated" = "#F6807F", "Downregulated" = "#9EB5F0", "Not Significant" = "#A7A7A7")

      n_comp <- length(unique(final_data$Comparison))

      total_plots <- ifelse(n_comp %% 2 == 0, n_comp, n_comp+1)

      if(total_plots==2) {
        volcano_plot <- ggplot(final_data, aes(x = LFC, y = -log10(P_adj), color = Name)) +
          geom_point(aes(color = Name), alpha = 0.6, size = 5) +
          facet_wrap(~Comparison, nrow = 1, ncol = 2, scales = "free") +
          scale_color_manual(values = color_group) +
          geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
          geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "red") +
          labs(
            x = "Log2 Fold Change",
            y = "-Log10(Adjusted P-value)",
            caption = paste0("ANCOM-BC2 is applied on ", lineage, " counts matrix with the prevalence cutoff of ", prevalence_cutoff, "% and counts cutoff of ", counts_cutoff, ".")
          ) +
          theme_classic() +
          geom_label_repel(aes(label = ifelse(Name %in% c("Upregulated", "Downregulated"), !!sym(lineage), "")), size = 5, color = "#2b71c2", box.padding = 0.5) +
          theme(
            axis.text.x = element_text(size = 15, face = "bold", colour = "#5B5DC7"),
            axis.text.y = element_text(size = 15, face = "bold", colour = "#5B5DC7"),
            legend.text = element_text(size = 15, face = "bold", colour = "#5B5DC7"),
            axis.title.x = element_text(size = 15, face = "bold", colour = "#5B5DC7"),
            axis.title.y = element_text(size = 15, face = "bold", colour = "#5B5DC7"),
            legend.title = element_text(size = 15, face = "bold", colour = "#5B5DC7"),
            axis.line = element_line(colour = "black", linewidth = 0.5, linetype = "solid" ),
            strip.text.x = element_text(size = 20, face = "bold", colour = "#5B5DC7"),
            strip.background = element_blank(),
            plot.caption = element_markdown(
              color = "#0F6E73", size = 15,
              margin = margin(20, 0, 10, 0), face = "bold",
              hjust = 0.5
            )
          )
      } else if(total_plots>2) {
        volcano_plot <- ggplot(final_data, aes(x = LFC, y = -log10(P_adj), color = Name)) +
          geom_point(aes(color = Name), alpha = 0.6, size = 5) +
          facet_wrap(~Comparison, nrow = total_plots/2, ncol = total_plots/2, scales = "free") +
          scale_color_manual(values = color_group) +
          geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
          geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "red") +
          labs(
            x = "Log2 Fold Change",
            y = "-Log10(Adjusted P-value)",
            caption = paste0("ANCOM-BC2 is applied on ", lineage, " counts matrix with the prevalence cutoff of ", prevalence_cutoff, "% and counts cutoff of ", counts_cutoff, ".")
          ) +
          theme_classic() +
          geom_label_repel(aes(label = ifelse(Name %in% c("Upregulated", "Downregulated"), !!sym(lineage), "")), size = 5, color = "#2b71c2", box.padding = 0.5) +
          theme(
            axis.text.x = element_text(size = 15, face = "bold", colour = "#5B5DC7"),
            axis.text.y = element_text(size = 15, face = "bold", colour = "#5B5DC7"),
            legend.text = element_text(size = 15, face = "bold", colour = "#5B5DC7"),
            axis.title.x = element_text(size = 15, face = "bold", colour = "#5B5DC7"),
            axis.title.y = element_text(size = 15, face = "bold", colour = "#5B5DC7"),
            legend.title = element_text(size = 15, face = "bold", colour = "#5B5DC7"),
            axis.line = element_line(colour = "black", linewidth = 0.5, linetype = "solid" ),
            strip.text.x = element_text(size = 20, face = "bold", colour = "#5B5DC7"),
            strip.background = element_blank(),
            plot.caption = element_markdown(
              color = "#0F6E73", size = 15,
              margin = margin(20, 0, 10, 0), face = "bold",
              hjust = 0.5
            )
          )
      }
      return(list(
        'volcano_plot' = volcano_plot,
        'daa_result' = final_data
      ))
    } else {
      
      return(NULL)
    
    }


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

  output$sampleinfo <- DT::renderDataTable({
    
    temp <- input_data_reactive()$data
    
    temp
    
  },
  options = list(paging = TRUE,
                 pageLength = 10,
                 scrollX = TRUE,
                 scrollY = TRUE,
                 autoWidth = FALSE,
                 server = TRUE,
                 rownames = FALSE,
                 initComplete = JS(
                    "function(settings, json) {",
                    "$(this.api().table().header()).css({'background-color': '#6C7AE0', 'color': '#FFFFFF', 'font-weight': 'bold'});",
                    "$(this.api().table().body()).find('tr.odd').css({'background-color': '#FFFFFF', 'color': '#000000'})",
                    "$(this.api().table().body()).find('tr.even').css({'background-color': '#F8F6FF', 'color': '#000000'})",
                    "}"
                  ),
                  drawCallback = JS(
                    "function(settings) {",
                    "$(this.api().table().body()).find('tr.odd').css({'background-color': '#FFFFFF', 'color': '#000000'});",
                    "$(this.api().table().body()).find('tr.even').css({'background-color': '#F8F6FF', 'color': '#000000'});",
                    "}"
                  )
                 )
                 )

  output$analysisoutput <- DT::renderDataTable({
  
    print("Data Engineering Output")

    if(route()=="Realtime")
    {
      
      req(cohort_analysis_list(), cohort_sample_list(), input$taxa, input_data_reactive(), input$prevalence_cutoff, input$abundance_cutoff)

      sample_list <- cohort_sample_list()$Barcode

      lineage <- input$taxa

      prevalence_cutoff <- input$prevalence_cutoff

      abundance_cutoff <- input$abundance_cutoff

      abundance_data <- cohort_realtime_analysis(cohort_analysis_list(), sample_list, lineage, prevalence_cutoff, abundance_cutoff)$counts_data

      abundance_val(abundance_data)

      abundance_data
    
    }
    
    else if(route()=="Example")
    {
      req(reactive_counts_data(), input$taxa)
                      
      lineage <- input$taxa

      abundance_data <- reactive_counts_data()

      abundance_val(abundance_data)

      abundance_data
    
    }

    else if(route()=="Offline")
    {
      req(reactive_counts_data())

      abundance_data <- reactive_counts_data()

      abundance_val(abundance_data)

      abundance_data
    
    }
  
  },
  options = list(paging = TRUE,
                 pageLength = 10,
                 scrollX = TRUE,
                 scrollY = TRUE,
                 autoWidth = TRUE,
                 server = TRUE,
                 rownames = FALSE,
                 initComplete = JS(
                    "function(settings, json) {",
                    "$(this.api().table().header()).css({'background-color': '#6C7AE0', 'color': '#FFFFFF', 'font-weight': 'bold'});",
                    "$(this.api().table().body()).find('tr.odd').css({'background-color': '#FFFFFF', 'color': '#000000'})",
                    "$(this.api().table().body()).find('tr.even').css({'background-color': '#F8F6FF', 'color': '#000000'})",
                    "}"
                  ),
                  drawCallback = JS(
                    "function(settings) {",
                    "$(this.api().table().body()).find('tr.odd').css({'background-color': '#FFFFFF', 'color': '#000000'});",
                    "$(this.api().table().body()).find('tr.even').css({'background-color': '#F8F6FF', 'color': '#000000'});",
                    "}"
                  )
                 )
                 
    )
    
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
        paste0(lineage,"_Counts_Data_", Sys.Date(), ".csv")
        },
      content = function(file) {
        write.csv(abundance_val(), file, row.names = FALSE, quote = FALSE)
      }
  )

  output$state <- renderText({
    state()
  })
  
  output$pores <- renderText({
    pores()
  })
  
  
  output$passed_reads <- renderText({
    
    req(input$barcode_select, processed_reads_df(), sample_list())

    if(!input$barcode_select %in% sample_list()$Barcode) {

      mean_length <- 0

    } else {
      
      selected_count <- processed_reads_df()[processed_reads_df()$Barcode == input$barcode_select, "Processed_Reads"]
    
      selected_count
    
    }
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

      taxa <- input$taxon_select

      df <- df %>% dplyr::select(c(!!sym(input$taxon_select), Counts)) %>% 
      dplyr::filter(!!sym(taxa) !="Unclassified") %>% dplyr::summarise(Total_Counts = sum(Counts))
    
      classified_reads <- df$Total_Counts

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

  output$plot_classification_plot <- renderPlotly({

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
      mutate(Abundance = (Counts/sum(Counts))*100)
      
      if(input$switch) {
        df <- df %>% arrange(desc(Abundance)) %>% head(n=num)
        df[[taxa]] = reorder(df[[taxa]], -df$Abundance)
        classification_plot <- ggplot(df, aes(x = !!sym(taxa), Abundance)) +
        geom_bar(stat = "identity", fill = "#557EAE") +
        theme_minimal() +
        ylab("Relative Abundance (%)") +
        theme(
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 15, face = "bold", colour = "#5B5DC7"),
          axis.line = element_line(colour = "#557EAE", linewidth = 1, linetype = "solid" ),
          strip.text.x = element_text(size = 15, face = "bold", colour = "#5B5DC7"),
          axis.text.y= element_text(size=15, face = "bold", colour = "#5B5DC7"),
          axis.text.x= element_text(size=15, face = "bold", colour = "#5B5DC7", angle = 90, vjust = 1, hjust=1),
          legend.position = "none",
          title = element_text(size = 15, face = "bold")
        )
      } else {
        df <- df %>% arrange(desc(Counts)) %>% head(n=num)
        df[[taxa]] = reorder(df[[taxa]], -df$Counts)
        classification_plot <- ggplot(df, aes(x = !!sym(taxa), Counts)) +
        geom_bar(stat = "identity", fill = "#557EAE") +
        theme_minimal() +
        ylab("Read Counts") +
        theme(
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 15, face = "bold", colour = "#5B5DC7"),
          axis.line = element_line(colour = "#557EAE", linewidth = 1, linetype = "solid" ),
          strip.text.x = element_text(size = 15, face = "bold", colour = "#5B5DC7"),
          axis.text.y= element_text(size=15, face = "bold", colour = "#5B5DC7"),
          axis.text.x= element_text(size=15, face = "bold", colour = "#5B5DC7", angle = 90, vjust = 1, hjust=1),
          legend.position = "none",
          title = element_text(size = 15, face = "bold")
        )
      }

      plot_taxa_bar(classification_plot)

      ggplotly(classification_plot, height = 600)

    }

  })

  output$plot_rarefaction_plot <- renderPlotly({
    
    df <- rarefaction_combined()

    req(nrow(df) > 0, input$barcode_select)

    df <- unique(df)

    if(!input$barcode_select %in% df$Sample_Id) {

      return(NULL)

    } else {

      barcode_df <- df %>% filter(Sample_Id == input$barcode_select)

      barcode_df <- rbind(data.frame(Sample_Id = input$barcode_select, Classified_reads = 0, Unique_taxa = 0), barcode_df)

      max_reads <- max(barcode_df$Classified_reads) + 10

      max_taxa <- max(barcode_df$Unique_taxa) + 2

      rarefaction_plot <- ggplot(barcode_df, aes(Classified_reads, Unique_taxa, color = "#5869B9")) +
                            geom_line(linewidth =  0.8, alpha = 0.7) +
                            theme_linedraw() +
                            ylab("Unique Species Count") +
                            xlab("Classified Reads") +
                            theme(
                            axis.title.x = element_text(size = 15, face = "bold", colour = "#5B5DC7"),
                            axis.title.y = element_text(size = 15, face = "bold", colour = "#5B5DC7"),
                            axis.line = element_line(colour = "#557EAE", linewidth = 1, linetype = "solid" ),
                            strip.text.x = element_text(size = 15, face = "bold", colour = "#5B5DC7"),
                            axis.text.y= element_text(size=15, face = "bold", colour = "#5B5DC7"),
                            axis.text.x= element_text(size=15, face = "bold", colour = "#5B5DC7"),
                            legend.position = "none"
                          ) +
                          scale_x_continuous(expand = expansion(mult = 0), limits = c(0, max_reads)) +
                          scale_y_continuous(expand = expansion(mult = 0), limits = c(0, max_taxa))

      plot_rarefaction_curve(rarefaction_plot)

      ggplotly(rarefaction_plot)

    }
  
  })

  output$plot_diversity_plot <- renderPlot({

    df <- diversity_combined()

    req(nrow(df) > 0, input$barcode_select)

    df <- unique(df)

    if(!input$barcode_select %in% df$Sample_Id) {

      return(NULL)

    } else {

      barcode_df <- df %>% filter(Sample_Id == input$barcode_select)

      origin_dataframe <- data.frame(Sample_Id = c(input$barcode_select, input$barcode_select), Classified_reads = c(0, 0), Diversity = c("Shannon", "Simpson"), Value = c(0, 0))

      barcode_df <- rbind(origin_dataframe, barcode_df)

      barcode_df <- unique(barcode_df)

      if(nrow(barcode_df)>=4) {
        
        shannon_df <- barcode_df %>% filter(Diversity == "Shannon")

        max_reads <- max(shannon_df$Classified_reads) + 10

        shannon_max_value <- max(shannon_df$Value) + 0.1

        shannon_diversity_plot <- ggplot(shannon_df, aes(Classified_reads, Value, color = "#5869B9")) +
                                  geom_line(linewidth =  0.8, alpha = 0.7) +
                                  theme_linedraw() +
                                  ylab("Diversity") +
                                  xlab("Classified Reads") +
                                  theme(
                                    axis.title.x = element_text(size = 15, face = "bold", colour = "#5B5DC7"),
                                    axis.title.y = element_text(size = 15, face = "bold", colour = "#5B5DC7"),
                                    axis.line = element_line(colour = "#557EAE", linewidth = 1, linetype = "solid" ),
                                    strip.text.x = element_text(size = 15, face = "bold", colour = "#5B5DC7"),
                                    axis.text.y= element_text(size=15, face = "bold", colour = "#5B5DC7"),
                                    axis.text.x= element_text(size=15, face = "bold", colour = "#5B5DC7"),
                                    legend.position = "none",
                                    plot.title = element_text(size = 18, face = "bold", colour = "#5B5DC7", hjust = 0.5)
                                ) +
                                scale_x_continuous(expand = expansion(mult = 0), limits = c(0, max_reads)) +
                                scale_y_continuous(expand = expansion(mult = 0), limits = c(0, shannon_max_value)) +
                                ggtitle("Shannon")

        simpson_df <- barcode_df %>% filter(Diversity == "Simpson")

        simpson_max_value <- max(simpson_df$Value) + 0.02

        simpson_diversity_plot <- barcode_df %>% filter(Diversity == "Simpson") %>%
                                  ggplot(aes(Classified_reads, Value, color = "#5869B9")) +
                                  geom_line(linewidth =  0.8, alpha = 0.7) +
                                  theme_linedraw() +
                                  ylab("Diversity") +
                                  xlab("Classified Reads") +
                                  theme(
                                  axis.title.x = element_text(size = 15, face = "bold", colour = "#5B5DC7"),
                                  axis.title.y = element_text(size = 15, face = "bold", colour = "#5B5DC7"),
                                  axis.line = element_line(colour = "#557EAE", linewidth = 1, linetype = "solid" ),
                                  strip.text.x = element_text(size = 15, face = "bold", colour = "#5B5DC7"),
                                  axis.text.y= element_text(size=15, face = "bold", colour = "#5B5DC7"),
                                  axis.text.x= element_text(size=15, face = "bold", colour = "#5B5DC7"),
                                  legend.position = "none",
                                  plot.title = element_text(size = 18, face = "bold", colour = "#5B5DC7", hjust = 0.5)
                                ) +
                                scale_x_continuous(expand = expansion(mult = 0), limits = c(0, max_reads)) +
                                scale_y_continuous(expand = expansion(mult = 0), limits = c(0, simpson_max_value)) +
                                ggtitle("Simpson")

        # diversity_facet <- as_ggplot(grid.grabExpr(grid.arrange(shannon_diversity_plot, simpson_diversity_plot, ncol=2)))

        diversity_facet <- plot_grid(shannon_diversity_plot, simpson_diversity_plot, ncol = 2, align = 'v', axis = 'lrtb')

        plot_diversity_curve(diversity_facet)

        diversity_facet
      }
    
    }

  }, height = 600)

  output$plot_read_length <- renderPlotly({

    req(input$barcode_select, sample_list(), hist_list())

    if(!input$barcode_select %in% sample_list()$Barcode) {

      return(NULL)

    }

    else {
      
      df <- hist_list()[[which(sample_list()$Barcode == input$barcode_select)]]

      df <- df %>% group_by(Bin) %>% summarise(Counts = sum(Counts)) %>% filter(Bin>=1000 & Bin<=2000) %>% 
      mutate(Selected = ifelse(Bin>=1400 & Bin<=1800, "Yes", "No"))

      cols <- c("Yes" = "#557EAE", "No" = "#ea7cac")

      histogram_plot <- ggplot(df, aes(Bin, Counts, fill = Selected)) +
      geom_bar(stat = "identity", color = "black") +
      scale_fill_manual(values = cols) +
      scale_x_continuous(expand=c(0,0)) +
      scale_y_continuous(expand=c(0,0)) +
      theme_minimal() +
      theme(
        axis.title.x = element_text(size = 14, face = "bold", colour = "#5B5DC7"),
        axis.title.y = element_text(size = 14, face = "bold", colour = "#5B5DC7"),
        strip.text.x = element_text(size = 14, face = "bold", colour = "#5B5DC7"),
        axis.text.y=element_text(size=14, face = "bold", colour = "#5B5DC7"),
        axis.text.x=element_text(size=14, face = "bold", colour = "#5B5DC7"),
        axis.line = element_line(colour = "#557EAE", linewidth = 0.5, linetype = "solid"),
        legend.title=element_blank(),
        legend.text=element_blank(),
        plot.title = element_text(hjust = 0.5, size = rel(2)),
        legend.position = "none"
      ) +
      xlab("Read Length") +
      ylab("Read Counts")
      
      plot_read_histogram(histogram_plot)
      
      ggplotly(histogram_plot)
      
    }
  
  })

  output$plot_q_score_plot <- renderPlotly({

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

      cols <- c("Yes" = "#557EAE", "No" = "#ea7cac")

      q_score_plot <- df %>% ggplot(aes(Q_Score, Frequency, fill = Selected)) +
        geom_bar(stat = "identity", color = "black") +
        theme_minimal() +
        xlab("Phred Score") +
        ylab("Read Counts") +
        theme(
          axis.title.x = element_text(size = 14, face = "bold", colour = "#5B5DC7"),
          axis.title.y = element_text(size = 14, face = "bold", colour = "#5B5DC7"),
          strip.text.x = element_text(size = 14, face = "bold", colour = "#5B5DC7"),
          axis.text.y= element_text(size=14, face = "bold", colour = "#5B5DC7"),
          axis.text.x= element_text(size=14, face = "bold", colour = "#5B5DC7"),
          axis.line = element_line(colour = "#557EAE", linewidth = 0.5, linetype = "solid" ),
          strip.background = element_blank(),
          legend.position = "none"
        ) +
        scale_fill_manual(values = cols)
      
      plot_q_histogram(q_score_plot)

      ggplotly(q_score_plot)

    }

  })

  output$taxa_table <- renderDT({

    req(input$barcode_select, classified_list(), input$taxon_select, classified_samples_list())

    if(!input$barcode_select %in% classified_samples_list()$Barcode) {

      return(NULL)

    }

    else {
      
      df <- classified_list()[[which(classified_samples_list()$Barcode == input$barcode_select)]]

      taxa <- input$taxon_select

      df <- df  %>% dplyr::select(c(!!sym(taxa), Counts)) %>% 
      filter(!!sym(taxa) != "Unclassified") %>% 
      group_by(!!sym(taxa)) %>% 
      summarise(Counts = sum(Counts)) %>% 
      mutate(Abundance = (Counts/sum(Counts))*100) %>%
      arrange(desc(Counts)) %>%
      select(!!sym(taxa), Counts, Abundance)

      taxa_count_table(df)

      df$Abundance = df$Abundance/100

      datatable(df, escape = FALSE, options = list(
        pageLength = 10,
        autoWidth = TRUE,
        rownames = FALSE,
        columnDefs = list(list(className = "dt-left", targets = 0:2)))) %>%
        formatPercentage("Abundance", 1) %>%
        formatStyle("Abundance",
        background = styleColorBar(df$Abundance, "steelblue",  -90),
        backgroundSize = '98% 88%',
        backgroundRepeat = 'no-repeat',
        backgroundPosition = "center")
    }

  }, ,height = 350, options = list(paging = TRUE,
                                   pageLength = 10,
                                   scrollX = TRUE,
                                   scrollY = TRUE,
                                   autoWidth = FALSE,
                                   server = TRUE,
                                   rownames = FALSE,
                                   initComplete = JS("function(settings, json) {",
                                                     "$(this.api().table().header()).css({'background-color': '#6C7AE0', 'color': '#FFFFFF', 'font-weight': 'bold'});",
                                                     "$(this.api().table().body()).find('tr.odd').css({'background-color': '#FFFFFF', 'color': '#000000'})",
                                                     "$(this.api().table().body()).find('tr.even').css({'background-color': '#F8F6FF', 'color': '#000000'})",
                                                     "}"),
                                   drawCallback = JS(
                                     "function(settings) {",
                                     "$(this.api().table().body()).find('tr.odd').css({'background-color': '#FFFFFF', 'color': '#000000'});",
                                     "$(this.api().table().body()).find('tr.even').css({'background-color': '#F8F6FF', 'color': '#000000'});",
                                     "}")
                 ))

  output$plot_stacked_barplot <- renderPlotly({

    if(route()=="Realtime")
    {
      
      req(cohort_analysis_list(), input$taxa, input$top_taxa, cohort_sample_list(), input$prevalence_cutoff, input$abundance_cutoff)

      sample_list <- cohort_sample_list()$Barcode

      top_n <- input$top_taxa

      lineage <- input$taxa

      prevalence_cutoff <- input$prevalence_cutoff

      abundance_cutoff <- input$abundance_cutoff

      stacked_df <- cohort_realtime_analysis(cohort_analysis_list(), sample_list, lineage, prevalence_cutoff, abundance_cutoff)$rel_abundance_renormalized_matrix

      if(is.null(input$toi) || length(input$toi) == 0) {
        stacked_barplot <- stacked_barplot_function(stacked_df, lineage, top_n)
      } else {
        required_species <- input$toi
        stacked_barplot <- stacked_subsampled_barplot_function(stacked_df, lineage, required_species)
      }
      
      plot_taxa_stacked(stacked_barplot)

      ggplotly(stacked_barplot)
    
    }

    else if(route()=="Example")
    {
      req(reactive_rel_abundance_matrix(), input$taxa, input$top_taxa, input_data_reactive())

      lineage <- input$taxa

      top_n <- input$top_taxa
      
      stacked_df <- reactive_rel_abundance_matrix()

      if(is.null(input$toi) || length(input$toi) == 0) {
        stacked_barplot <- stacked_barplot_function(stacked_df, lineage, top_n)
      } else {
        required_species <- input$toi
        stacked_barplot <- stacked_subsampled_barplot_function(stacked_df, lineage, required_species)
      }

      plot_taxa_stacked(stacked_barplot)

      ggplotly(stacked_barplot)
    
    }

    else if(route()=="Offline")
    {
      req(reactive_rel_abundance_matrix(), input$taxa, input$top_taxa, input_data_reactive())

      lineage <- input$taxa

      top_n <- input$top_taxa

      stacked_df <- reactive_rel_abundance_matrix()

      if(is.null(input$toi) || length(input$toi) == 0) {
        stacked_barplot <- stacked_barplot_function(stacked_df, lineage, top_n)
      } else {
        required_species <- input$toi
        stacked_barplot <- stacked_subsampled_barplot_function(stacked_df, lineage, required_species)
      }

      plot_taxa_stacked(stacked_barplot)

      ggplotly(stacked_barplot, tooltip = "text")
    }

  })

  output$plot_boxplot <- renderPlot({
    
    if(route()=="Realtime")
    {

      req(cohort_analysis_list(), cohort_sample_list(), input$taxa, input_data_reactive(), input$prevalence_cutoff, input$abundance_cutoff)

      sample_list <- cohort_sample_list()$Barcode

      lineage <- input$taxa

      sample_metadata <- input_data_reactive()$data

      prevalence_cutoff <- input$prevalence_cutoff

      abundance_cutoff <- input$abundance_cutoff
      
      counts_data <- cohort_realtime_analysis(cohort_analysis_list(), sample_list, lineage, prevalence_cutoff, abundance_cutoff)$counts_data

      counts_data_long <- counts_data %>% pivot_longer(cols = -!!sym(lineage), names_to = "Sample", values_to = "Counts")

      rarefaction_cutoff <- counts_data_long %>% group_by(Sample) %>% 
        summarise(Total = sum(Counts), Singletons = sum(Counts==1), GC = 100*(1-(Singletons/Total))) %>% 
        filter(GC>=95) %>% summarise(Min = min(Total)) %>% pull(Min)

      counts_data_wide <- counts_data_long %>% pivot_wider(names_from = !!sym(lineage), values_from = Counts, values_fill = 0) %>% 
        as.data.frame()

      counts_data_wide <- counts_data_wide %>% dplyr::select(Sample, everything())

      rownames(counts_data_wide) <- counts_data_wide$Sample

      counts_data_wide <- counts_data_wide[,-1]

      counts_data_wide <- round(counts_data_wide, digits = 0)

      diversity_facet <- diversity_boxplot_function(counts_data_wide, lineage, sample_metadata, rarefaction_cutoff)

      plot_diversity_box(diversity_facet)

      diversity_facet
    
    }

    else if(route()=="Example")
    {
      req(reactive_counts_data(), input$taxa, input_data_reactive())
      
      lineage <- input$taxa

      sample_metadata <- input_data_reactive()$data
      
      counts_data <- reactive_counts_data()

      counts_data_long <- counts_data %>% pivot_longer(cols = -!!sym(lineage), names_to = "Sample", values_to = "Counts")

      rarefaction_cutoff <- counts_data_long %>% group_by(Sample) %>% 
        summarise(Total = sum(Counts), Singletons = sum(Counts==1), GC = 100*(1-(Singletons/Total))) %>% 
        filter(GC>=95) %>% summarise(Min = min(Total)) %>% pull(Min)

      counts_data_wide <- counts_data_long %>% pivot_wider(names_from = !!sym(lineage), values_from = Counts, values_fill = 0) %>% 
        as.data.frame()

      counts_data_wide <- counts_data_wide %>% dplyr::select(Sample, everything())

      rownames(counts_data_wide) <- counts_data_wide$Sample

      counts_data_wide <- counts_data_wide[,-1]

      counts_data_wide <- round(counts_data_wide, digits = 0)

      diversity_facet <- diversity_boxplot_function(counts_data_wide, lineage, sample_metadata, rarefaction_cutoff)

      plot_diversity_box(diversity_facet)

      diversity_facet
    
    }
    
    else if(route() == "Offline")
    {
      req(reactive_counts_data(), input_data_reactive(), input$taxa, input$prevalence_cutoff, input$abundance_cutoff)

      lineage <- input$taxa

      sample_metadata <- input_data_reactive()$data
      
      counts_data <- reactive_counts_data()

      counts_data_long <- counts_data %>% pivot_longer(cols = -!!sym(lineage), names_to = "Sample", values_to = "Counts")

      rarefaction_cutoff <- counts_data_long %>% group_by(Sample) %>% 
        summarise(Total = sum(Counts), Singletons = sum(Counts==1), GC = 100*(1-(Singletons/Total))) %>% 
        filter(GC>=95) %>% summarise(Min = min(Total)) %>% pull(Min)

      counts_data_wide <- counts_data_long %>% pivot_wider(names_from = !!sym(lineage), values_from = Counts, values_fill = 0) %>% 
        as.data.frame()

      counts_data_wide <- counts_data_wide %>% dplyr::select(Sample, everything())

      rownames(counts_data_wide) <- counts_data_wide$Sample

      counts_data_wide <- counts_data_wide[,-1]

      counts_data_wide <- round(counts_data_wide, digits = 0)

      diversity_facet <- diversity_boxplot_function(counts_data_wide, lineage, sample_metadata, rarefaction_cutoff)

      plot_diversity_box(diversity_facet)

      diversity_facet
    
    }
  
  }, height = 500)

  output$plot_pcoa <- renderPlot({

    if(route()=="Realtime")
    {
      req(cohort_analysis_list(), cohort_sample_list(), input$taxa, input_data_reactive(), input$prevalence_cutoff, input$abundance_cutoff)

      lineage <- input$taxa
      
      sample_metadata <- input_data_reactive()$data

      prevalence_cutoff <- input$prevalence_cutoff

      abundance_cutoff <- input$abundance_cutoff

      sample_list <- cohort_sample_list()$Barcode

      matrix <- cohort_realtime_analysis(cohort_analysis_list(), sample_list, lineage, prevalence_cutoff, abundance_cutoff)$rel_abundance_renormalized_matrix

      pcoa_plot <- diversity_pcoa_function(matrix, lineage, sample_metadata, prevalence_cutoff, abundance_cutoff)

      plot_pcoa_dot(pcoa_plot)

      pcoa_plot
    
    }

    else if(route()=="Example")
    {
      
      req(reactive_rel_abundance_matrix(), input_data_reactive(), input$taxa, input$prevalence_cutoff, input$abundance_cutoff)
                      
      lineage <- input$taxa

      sample_metadata <- input_data_reactive()$data

      matrix <- reactive_rel_abundance_matrix()

      prevalence_cutoff <- input$prevalence_cutoff

      abundance_cutoff <- input$abundance_cutoff

      pcoa_plot <- diversity_pcoa_function(matrix, lineage, sample_metadata, prevalence_cutoff, abundance_cutoff)

      plot_pcoa_dot(pcoa_plot)

      pcoa_plot

    }

    else if(route() == "Offline")
    {
      req(reactive_rel_abundance_matrix(), input_data_reactive(), input$taxa, input$prevalence_cutoff, input$abundance_cutoff)

      lineage <- input$taxa

      sample_metadata <- input_data_reactive()$data

      prevalence_cutoff <- input$prevalence_cutoff

      abundance_cutoff <- input$abundance_cutoff

      matrix <- reactive_rel_abundance_matrix()

      pcoa_plot <- diversity_pcoa_function(matrix, lineage, sample_metadata, prevalence_cutoff, abundance_cutoff)

      plot_pcoa_dot(pcoa_plot)

      pcoa_plot
    
    }
  
  }, height = 600)

  output$plot_nmds <- renderPlot({
    
    if(route()=="Realtime")
    {
      req(cohort_analysis_list(), cohort_sample_list(), input$taxa, input_data_reactive(), input$prevalence_cutoff, input$abundance_cutoff)

      lineage <- input$taxa
      
      sample_metadata <- input_data_reactive()$data

      prevalence_cutoff <- input$prevalence_cutoff

      abundance_cutoff <- input$abundance_cutoff

      sample_list <- cohort_sample_list()$Barcode

      matrix <- cohort_realtime_analysis(cohort_analysis_list(), sample_list, lineage, prevalence_cutoff, abundance_cutoff)$rel_abundance_renormalized_matrix
    
      nmds_plot <- diversity_nmds_function(matrix, lineage, sample_metadata, prevalence_cutoff, abundance_cutoff)
      
      plot_nmds_dot(nmds_plot)

      nmds_plot
    
    }

    else if(route()=="Example")
    {
      req(reactive_rel_abundance_matrix(), input_data_reactive(), input$taxa, input$prevalence_cutoff, input$abundance_cutoff)
      
      lineage <- input$taxa

      sample_metadata <- input_data_reactive()$data

      prevalence_cutoff <- input$prevalence_cutoff

      abundance_cutoff <- input$abundance_cutoff

      matrix <- reactive_rel_abundance_matrix()

      nmds_plot <- diversity_nmds_function(matrix, lineage, sample_metadata, prevalence_cutoff, abundance_cutoff)
      
      plot_nmds_dot(nmds_plot)

      nmds_plot

    }

    else if(route() == "Offline")
    {
      req(reactive_rel_abundance_matrix(), input_data_reactive(), input$taxa, input$prevalence_cutoff, input$abundance_cutoff)

      lineage <- input$taxa

      sample_metadata <- input_data_reactive()$data

      prevalence_cutoff <- input$prevalence_cutoff

      abundance_cutoff <- input$prevalence_cutoff

      matrix <- reactive_rel_abundance_matrix()

      nmds_plot <- diversity_nmds_function(matrix, lineage, sample_metadata, prevalence_cutoff, abundance_cutoff)
      
      plot_nmds_dot(nmds_plot)

      nmds_plot

    }
    
  }, height = 600)

  output$plot_pca <- renderPlot({
    
    if(route()=="Realtime")
    {
      req(cohort_analysis_list(), cohort_sample_list(), input$taxa, input_data_reactive(), input$prevalence_cutoff, input$abundance_cutoff, input$biplot_taxa)

      lineage <- input$taxa

      sample_metadata <- input_data_reactive()$data

      prevalence_cutoff <- input$prevalence_cutoff

      abundance_cutoff <- input$abundance_cutoff

      top_taxa <- input$biplot_taxa

      sample_list <- cohort_sample_list()$Barcode

      matrix <- cohort_realtime_analysis(cohort_analysis_list(), sample_list, lineage, prevalence_cutoff, abundance_cutoff)$rel_abundance_renormalized_matrix

      pca_plot <- diversity_pca_function(matrix, lineage, sample_metadata, prevalence_cutoff, abundance_cutoff, top_taxa)

      plot_pca_dot(pca_plot)
      
      pca_plot
    
    }

    else if(route()=="Example")
    {
      req(reactive_rel_abundance_matrix(), input_data_reactive(), input$taxa, input$prevalence_cutoff, input$abundance_cutoff, input$biplot_taxa)

      lineage <- input$taxa

      sample_metadata <- input_data_reactive()$data

      prevalence_cutoff <- input$prevalence_cutoff

      abundance_cutoff <- input$abundance_cutoff

      top_taxa <- input$biplot_taxa

      matrix <- reactive_rel_abundance_matrix()

      pca_plot <- diversity_pca_function(matrix, lineage, sample_metadata, prevalence_cutoff, abundance_cutoff, top_taxa)

      plot_pca_dot(pca_plot)
      
      pca_plot  
    
    }

    else if(route() == "Offline")
    {
      req(reactive_rel_abundance_matrix(), input_data_reactive(), input$taxa, input$prevalence_cutoff, input$abundance_cutoff, input$biplot_taxa)

      lineage <- input$taxa

      sample_metadata <- input_data_reactive()$data

      prevalence_cutoff <- input$prevalence_cutoff

      abundance_cutoff <- input$abundance_cutoff

      top_taxa <- input$biplot_taxa

      matrix <- reactive_rel_abundance_matrix()

      pca_plot <- diversity_pca_function(matrix, lineage, sample_metadata, prevalence_cutoff, abundance_cutoff, top_taxa)

      plot_pca_dot(pca_plot)
      
      pca_plot

    }

  }, height = 600)

  output$permanova_data <- DT::renderDataTable({
    
    if(route()=="Realtime")
    {
      req(cohort_analysis_list(), cohort_sample_list(), input$taxa, input_data_reactive(), input$realtime_control, input$prevalence_cutoff, input$abundance_cutoff)

      lineage <- input$taxa

      sample_metadata <- input_data_reactive()$data

      sample_list <- cohort_sample_list()$Barcode

      prevalence_cutoff <- input$prevalence_cutoff

      abundance_cutoff <- input$abundance_cutoff
      
      matrix <- cohort_realtime_analysis(cohort_analysis_list(), sample_list, lineage, prevalence_cutoff, abundance_cutoff)$rel_abundance_renormalized_matrix

      control_group <- input$realtime_control

      permanova_res <- diversity_permanova_function(matrix, lineage, sample_metadata, control_group)

      table_permanova(permanova_res)

      permanova_res
    
    }

    else if(route()=="Example")
    {
      req(reactive_rel_abundance_matrix(), input_data_reactive(), input$taxa, input$prevalence_cutoff, input$abundance_cutoff, input$control)

      lineage <- input$taxa

      sample_metadata <- input_data_reactive()$data

      prevalence_cutoff <- input$prevalence_cutoff

      abundance_cutoff <- input$abundance_cutoff

      control_group <- input$control
      
      matrix <- reactive_rel_abundance_matrix()
      
      permanova_res <- diversity_permanova_function(matrix, lineage, sample_metadata, control_group)

      table_permanova(permanova_res)

      permanova_res

    }

    else if(route() == "Offline")
    {
      req(reactive_rel_abundance_matrix(), input_data_reactive(), input$taxa, input$control)

      lineage <- input$taxa

      sample_metadata <- input_data_reactive()$data

      control_group <- input$control
      
      matrix <- reactive_rel_abundance_matrix()
      
      permanova_res <- diversity_permanova_function(matrix, lineage, sample_metadata, control_group)

      table_permanova(permanova_res)

      permanova_res
    
    }
  },
  options = list(paging = TRUE,
                 pageLength = 10,
                 scrollX = TRUE,
                 scrollY = TRUE,
                 autoWidth = FALSE,
                 server = TRUE,
                 rownames = FALSE,
                 initComplete = JS(
                    "function(settings, json) {",
                    "$(this.api().table().header()).css({'background-color': '#6C7AE0', 'color': '#FFFFFF', 'font-weight': 'bold'});",
                    "$(this.api().table().body()).find('tr.odd').css({'background-color': '#FFFFFF', 'color': '#000000'})",
                    "$(this.api().table().body()).find('tr.even').css({'background-color': '#F8F6FF', 'color': '#000000'})",
                    "}"
                  ),
                  drawCallback = JS(
                    "function(settings) {",
                    "$(this.api().table().body()).find('tr.odd').css({'background-color': '#FFFFFF', 'color': '#000000'});",
                    "$(this.api().table().body()).find('tr.even').css({'background-color': '#F8F6FF', 'color': '#000000'});",
                    "}"
                  )
                 )
  )

  output$plot_heatmap <- renderPlot({
    
    if(route()=="Realtime")
    {
      req(cohort_analysis_list(), cohort_sample_list(), input$taxa, input_data_reactive(), input$prevalence_cutoff, input$abundance_cutoff)

      lineage <- input$taxa

      sample_metadata <- input_data_reactive()$data

      sample_list <- classified_samples_list()$Barcode

      prevalence_cutoff <- input$prevalence_cutoff

      abundance_cutoff <- input$abundance_cutoff
      
      matrix <- cohort_realtime_analysis(cohort_analysis_list(), sample_list, lineage, prevalence_cutoff, abundance_cutoff)$rel_abundance_renormalized_matrix

      heatmap_ggplot <- diversity_heatmap_function(matrix, lineage, sample_metadata)

      plot_real_heatmap(heatmap_ggplot)

      heatmap_ggplot
    
    }

    else if(route()=="Example")
    {
      req(reactive_rel_abundance_matrix(), input_data_reactive(), input$taxa, input$prevalence_cutoff, input$abundance_cutoff)

      lineage <- input$taxa

      sample_metadata <- input_data_reactive()$data

      prevalence_cutoff <- input$prevalence_cutoff

      abundance_cutoff <- input$abundance_cutoff
      
      matrix <- reactive_rel_abundance_matrix()

      heatmap_ggplot <- diversity_heatmap_function(matrix, lineage, sample_metadata)

      plot_real_heatmap(heatmap_ggplot)

      heatmap_ggplot
    
    }

    else if(route() == "Offline")
    {
      req(reactive_rel_abundance_matrix(), input$taxa, input_data_reactive())

      lineage <- input$taxa

      sample_metadata <- input_data_reactive()$data
      
      matrix <- reactive_rel_abundance_matrix()

      heatmap_ggplot <- diversity_heatmap_function(matrix, lineage, sample_metadata)

      plot_real_heatmap(heatmap_ggplot)

      heatmap_ggplot
    
    }
  
  }, height = 1200)

  output$plot_volcano <- renderPlot({
    
    if(route()=="Realtime")
    {
      req(realtime_daa_data(), realtime_daa_metadata(), realtime_taxa_group(), realtime_daa_counts_cutoff(), realtime_daa_prev_cutoff(), input$realtime_control)

      res_dunn <- realtime_daa_data()

      metadata <- realtime_daa_metadata()

      lineage <- realtime_taxa_group()

      prevalence_cutoff <- realtime_daa_prev_cutoff()

      counts_cutoff <- realtime_daa_counts_cutoff()

      control <- input$realtime_control

      total_groups <- length(levels(metadata$Group))

      req_start_values <- c("passed_ss_Group", "diff_Group", "lfc_Group", "q_Group")

      out_col_values <- c("Sensitive", "Significance", "LFC", "P_adj")

      req_columns <- outer(levels(metadata$Group)[2:total_groups], req_start_values, function(x,y) paste0(y, x)) %>% 
        as.vector()

      filtered_dunn <- res_dunn[,c(1,which(colnames(res_dunn) %in% req_columns))]

      long_data_list <- list()

      for(i in 1:length(req_start_values)) {
        long_data_list[[i]] <- filtered_dunn %>% dplyr::select(c(taxon, starts_with(req_start_values[i]))) %>% 
          pivot_longer(cols = -taxon, names_to = "Comparison", values_to = out_col_values[i])

        long_data_list[[i]]$Comparison <- sub(req_start_values[i], "", long_data_list[[i]]$Comparison)
      }

      final_data <- purrr::reduce(long_data_list, left_join, by=c("taxon", "Comparison"))

      final_data <- final_data %>% dplyr::select(taxon, everything())

      colnames(final_data)[1] <- lineage

      final_data <- final_data %>% filter(Sensitive==TRUE)

      final_data$Comparison <- paste0(control, " - ", final_data$Comparison)

      final_data$Name <- ifelse(final_data$LFC < -1 & final_data$P_adj < 0.05, "Downregulated",
                                ifelse(final_data$LFC > 1 & final_data$P_adj < 0.05, "Upregulated", "Not Significant"))
      
      final_data$Name <- factor(final_data$Name, levels = c("Upregulated", "Downregulated", "Not Significant"))

      color_group <- c("Upregulated" = "#F6807F", "Downregulated" = "#9EB5F0", "Not Significant" = "#A7A7A7")

      n_comp <- length(unique(final_data$Comparison))

      total_plots <- ifelse(n_comp %% 2 == 0, n_comp, n_comp+1)

      if(total_plots==2) {
        volcano_plot <- ggplot(final_data, aes(x = LFC, y = -log10(P_adj), color = Name)) +
          geom_point(aes(color = Name), alpha = 0.6, size = 5) +
          facet_wrap(~Comparison, nrow = 1, ncol = 2, scales = "free") +
          scale_color_manual(values = color_group) +
          geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
          geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "red") +
          labs(
            x = "Log2 Fold Change",
            y = "-Log10(Adjusted P-value)",
            caption = paste0("ANCOM-BC2 is applied on ", lineage, " counts matrix with the prevalence cutoff of ", prevalence_cutoff, "% and counts cutoff of ", counts_cutoff, ".")
          ) +
          theme_classic() +
          geom_label_repel(aes(label = ifelse(Name %in% c("Upregulated", "Downregulated"), !!sym(lineage), "")), size = 5, color = "#2b71c2", box.padding = 0.5) +
          theme(
            axis.text.x = element_text(size = 15, face = "bold", colour = "#5B5DC7"),
            axis.text.y = element_text(size = 15, face = "bold", colour = "#5B5DC7"),
            legend.text = element_text(size = 15, face = "bold", colour = "#5B5DC7"),
            axis.title.x = element_text(size = 15, face = "bold", colour = "#5B5DC7"),
            axis.title.y = element_text(size = 15, face = "bold", colour = "#5B5DC7"),
            legend.title = element_text(size = 15, face = "bold", colour = "#5B5DC7"),
            axis.line = element_line(colour = "black", linewidth = 0.5, linetype = "solid" ),
            strip.text.x = element_text(size = 20, face = "bold", colour = "#5B5DC7"),
            strip.background = element_blank(),
            plot.caption = element_markdown(
              color = "#0F6E73", size = 15,
              margin = margin(20, 0, 10, 0), face = "bold",
              hjust = 0.5
            )
          )
      } else if(total_plots>2) {
        volcano_plot <- ggplot(final_data, aes(x = LFC, y = -log10(P_adj), color = Name)) +
          geom_point(aes(color = Name), alpha = 0.6, size = 5) +
          facet_wrap(~Comparison, nrow = total_plots/2, ncol = total_plots/2, scales = "free") +
          scale_color_manual(values = color_group) +
          geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
          geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "red") +
          labs(
            x = "Log2 Fold Change",
            y = "-Log10(Adjusted P-value)",
            caption = paste0("ANCOM-BC2 is applied on ", lineage, " counts matrix with the prevalence cutoff of ", prevalence_cutoff, "% and counts cutoff of ", counts_cutoff, ".")
          ) +
          theme_classic() +
          geom_label_repel(aes(label = ifelse(Name %in% c("Upregulated", "Downregulated"), !!sym(lineage), "")), size = 5, color = "#2b71c2", box.padding = 0.5) +
          theme(
            axis.text.x = element_text(size = 15, face = "bold", colour = "#5B5DC7"),
            axis.text.y = element_text(size = 15, face = "bold", colour = "#5B5DC7"),
            legend.text = element_text(size = 15, face = "bold", colour = "#5B5DC7"),
            axis.title.x = element_text(size = 15, face = "bold", colour = "#5B5DC7"),
            axis.title.y = element_text(size = 15, face = "bold", colour = "#5B5DC7"),
            legend.title = element_text(size = 15, face = "bold", colour = "#5B5DC7"),
            axis.line = element_line(colour = "black", linewidth = 0.5, linetype = "solid" ),
            strip.text.x = element_text(size = 20, face = "bold", colour = "#5B5DC7"),
            strip.background = element_blank(),
            plot.caption = element_markdown(
              color = "#0F6E73", size = 15,
              margin = margin(20, 0, 10, 0), face = "bold",
              hjust = 0.5
            )
          ) 
      }

      plot_daa_volcano(volcano_plot)

      table_ancombc_daa(final_data)

      volcano_plot
    
    }
    else if(route()=="Example")
    {
      req(reactive_counts_matrix(), input_data_reactive(), input$taxa, input$prevalence_cutoff, input$counts_cutoff, input$abundance_cutoff, input$control)

      lineage <- input$taxa

      sample_metadata <- input_data_reactive()$data

      prevalence_cutoff <- input$prevalence_cutoff

      abundance_cutoff <- input$abundance_cutoff

      counts_cutoff <- input$counts_cutoff

      control <- input$control

      matrix <- reactive_counts_matrix()

      req(daa_volcano_plot(matrix, lineage, sample_metadata, control, prevalence_cutoff, counts_cutoff))

      volcano_plot <- daa_volcano_plot(matrix, lineage, sample_metadata, control, prevalence_cutoff, counts_cutoff)$volcano_plot

      daa_result <- daa_volcano_plot(matrix, lineage, sample_metadata, control, prevalence_cutoff, counts_cutoff)$daa_result

      plot_daa_volcano(volcano_plot)

      table_ancombc_daa(daa_result)

      volcano_plot
    }
    else if(route()=="Offline")
    {
      req(reactive_counts_matrix(), input_data_reactive(), input$taxa, input$prevalence_cutoff, input$counts_cutoff, input$control)
      
      lineage <- input$taxa
      
      sample_metadata <- input_data_reactive()$data
      
      prevalence_cutoff <- input$prevalence_cutoff
      
      counts_cutoff <- input$counts_cutoff
      
      control <- input$control
      
      matrix <- reactive_counts_matrix()

      req(daa_volcano_plot(matrix, lineage, sample_metadata, control, prevalence_cutoff, counts_cutoff))
      
      volcano_plot <- daa_volcano_plot(matrix, lineage, sample_metadata, control, prevalence_cutoff, counts_cutoff)$volcano_plot

      daa_result <- daa_volcano_plot(matrix, lineage, sample_metadata, control, prevalence_cutoff, counts_cutoff)$daa_result
      
      plot_daa_volcano(volcano_plot)

      table_ancombc_daa(daa_result)
      
      volcano_plot
    }
  }, height = 600)

  output$download_readlength_plot <- downloadHandler(
    req(input$barcode_select),
    filename = function() {
      paste0("Read_Length_Histogram_", input$barcode_select, "_", Sys.Date(), "_", format(Sys.time(), "%H-%M-%S"),".pdf")
    },
    content = function(file) {
      ggsave(file, plot_read_histogram(), width = 13.69, height = 8.27,
      units = "in", dpi = 600, bg = "white", device = "pdf")
    }
  )

  output$download_qscore_plot <- downloadHandler(
    req(input$barcode_select),
    filename = function() {
      paste0("Phred_Score_Histogram_", input$barcode_select, "_", Sys.Date(), "_", format(Sys.time(), "%H-%M-%S"),".pdf")
    },
    content = function(file) {
      ggsave(file, plot_q_histogram(), width = 13.69, height = 8.27,
      units = "in", dpi = 600, bg = "white", device = "pdf")
    }
  )

  output$download_classification_plot <- downloadHandler(
    req(input$barcode_select, input$taxon_select),
    filename = function() {
      paste0(input$taxon_select,"_Classification_", input$barcode_select, "_", Sys.Date(), "_", format(Sys.time(), "%H-%M-%S"),".pdf")
    },
    content = function(file) {
      ggsave(file, plot_taxa_bar(), width = 13.69, height = 8.27,
      units = "in", dpi = 600, bg = "white", device = "pdf")
    }
  )

  output$download_taxa_table <- downloadHandler(
    req(input$barcode_select, input$taxon_select),
    filename = function() {
      paste0(input$taxon_select, "_Counts_", input$barcode_select, "_", Sys.Date(), "_", format(Sys.time(), "%H-%M-%S"), ".csv")
    },
    content = function(file) {
      write.table(taxa_count_table(), file, row.names = FALSE, sep = "\t", quote = FALSE, col.names = TRUE)
    }
  )

  output$download_rarefaction_curve <- downloadHandler(
    req(input$taxon_select),
    filename = function() {
      paste0("Realtime_", input$taxon_select, "_Rarefaction_curve_", Sys.Date(), "_", format(Sys.time(), "%H-%M-%S"), ".pdf")
    },
    content = function(file) {
      ggsave(file, plot_rarefaction_curve(), width = 13.69, height = 8.27,
      units = "in", dpi = 600, bg = "white", device = "pdf")
    }
  )

  output$download_diversity_curve <- downloadHandler(
    req(input$taxon_select),
    filename = function() {
      paste0("Realtime_", input$taxon_select, "_Diversity_curve_", Sys.Date(), "_", format(Sys.time(), "%H-%M-%S"), ".pdf")
    },
    content = function(file) {
      ggsave(file, plot_diversity_curve(), width = 13.69, height = 8.27,
      units = "in", dpi = 600, bg = "white", device = "pdf")
    }
  )

  output$download_stacked_barplot <- downloadHandler(
      filename = function() {
        req(input$taxa)
        lineage <- input$taxa
        paste("Stacked_bar_plot_", lineage,"_", Sys.Date(), ".pdf", sep = "")
      },
      content = function(file) {
        ggsave(file, plot_taxa_stacked(),
               width = 18, height = 9, units = "in", dpi = "retina", bg = "white", device = "pdf")
      }
    )
    
    output$download_boxplot <- downloadHandler(
      filename = function() {
        req(input$taxa)
        lineage <- input$taxa
        paste("Alpha_Diversity_plot_", lineage, "_", Sys.Date(), ".pdf", sep = "")
      },
      content = function(file) {
        ggsave(file, plot_diversity_box(),
               width = 13.69, height = 8.27, units = "in", dpi = 600, bg = "white", device = "pdf")
        
      }
    )
    
    output$download_pcoa_plot <- downloadHandler(
      filename = function() {
        req(input$taxa)
        lineage <- input$taxa
        paste("PCoA_plot_", lineage, "_", Sys.Date(), ".pdf", sep = "")
      },
      content = function(file) {
        ggsave(file, plot_pcoa_dot(),
               width = 13.69, height = 8.27, units = "in", dpi = 600, bg = "white", device = "pdf")
      }
    )
    
    output$download_nmds_plot <- downloadHandler(
      filename = function() {
        req(input$taxa)
        lineage <- input$taxa
        paste("NMDS_plot_", lineage, "_", Sys.Date(), ".pdf", sep = "")
      },
      content = function(file) {
        
        ggsave(file, plot_nmds_dot(),
               width = 13.69, height = 8.27, units = "in", dpi = 600, bg = "white", device = "pdf")
        
      }
    )
    
    output$download_pca_plot <- downloadHandler(
      filename = function() {
        req(input$taxa)
        lineage <- input$taxa
        paste("PCA_plot_", lineage, "_", Sys.Date(), ".pdf", sep = "")
      },
      content = function(file) {
        ggsave(file, plot_pca_dot(),
               width = 13.69, height = 8.27, units = "in", dpi = 600, bg = "white", device = "pdf")
        }
    )
    
    output$download_heatmap <- downloadHandler(
      filename = function() {
        req(input$taxa)
        lineage <- input$taxa
        paste("HeatMap_", lineage, "_", Sys.Date(), ".pdf", sep = "")
      },
      content = function(file) {
        ggsave(file, plot_real_heatmap(),
               width = 10, height = 12, units = "in", dpi = 600, bg = "white", device = "pdf")
        
      }
    )

    output$download_permanova_csv <- downloadHandler(
      filename = function() {
        req(input$taxa)
        lineage <- input$taxa
        paste0("PERMANOVA_Result_", lineage, "_",Sys.Date(), ".csv")
        },
      content = function(file) {
        write.csv(table_permanova(), file, row.names = FALSE, quote = FALSE)
      }
    )
  
    output$download_daa <- downloadHandler(
      filename = function() {
        req(input$taxa)
        lineage <- input$taxa
        paste("DAA_Volcano_", lineage, "_", Sys.Date(), ".pdf", sep="")
      },
      content = function(file) {
        ggsave(file, plot_daa_volcano(),
                width = 23.69, height = 18.27, units = "in", dpi = 600, bg = "white", device = "pdf")
      }
    )

    output$download_daa_csv <- downloadHandler(
      filename = function() {
        req(input$taxa)
        lineage <- input$taxa
        paste("DAA_ANCOMBC2_", lineage, "_final_", Sys.Date(), ".csv", sep = "")
      },
      content = function(file) {
        write.csv(table_ancombc_daa(), file, row.names = FALSE, quote = FALSE)
      }
    )


}


shinyApp(ui = ui, server = server)