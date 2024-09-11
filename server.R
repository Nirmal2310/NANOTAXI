options(shiny.maxRequestSize = 100*1024^2)

source("packages.R")

print(sessionInfo())

# conda_path <- system("if [ $(which conda | grep 'condabin') ]; then conda_path=$(which conda | sed 's/\/condabin.*$//g'); else conda_path=$(which conda | sed 's/\/bin.*$//g'); fi && echo $conda_path", intern = TRUE)

# Sys.setenv(RETICULATE_PYTHON = paste0(conda_path,"/envs/minknow_api/bin/python"))

# print(paste0(conda_path,"/envs/minknow_api/bin/python"))

server <- function(input, output, session) {
    
    source("server-input.R", local = TRUE)
    
    source("server-cohort-analysis.R", local = TRUE)

    mode <- reactiveVal(NULL)
    
    ########### Reactive Values For Real Time Analysis ###########

    state <- reactiveVal()
  
    is_running <- reactiveVal(FALSE)
    
    timer_10s <- reactiveTimer(10000)
    
    timer_30s <- reactiveTimer(30000)

    initial_delay_done <- reactiveVal(FALSE)

    status_checked <- reactiveVal(FALSE)
    
    pores <- reactiveVal()
    
    passed_counts <- reactiveVal()

    sample_list <- reactiveVal()

    classified_samples_list <- reactiveVal()

    hist_list <- reactiveVal()

    qual_list <- reactiveVal()

    mean_df <- reactiveVal()

    taxa_count_table <- reactiveVal()

    classified_list <- reactiveVal()

    plot_read_histogram <- reactiveVal()

    plot_q_histogram <- reactiveVal()

    plot_taxa_bar <- reactiveVal()

    ########### Reactive Values For Realtime Cohort Analysis ########

    real_rel_abundance_val <- reactiveVal()

    real_abundance_matrix_val <- reactiveVal()

    real_rel_abundance_filtered_val <- reactiveVal()

    real_rel_abundance_filtered_matrix_val <- reactiveVal()


    ######## Reactive Values For Offline Cohort Analysis ###########
    
    rel_abundance_val <- reactiveVal()
    
    plot_taxa_stacked <- reactiveVal()

    plot_diversity_box <- reactiveVal()

    plot_pcoa_dot <- reactiveVal()

    plot_nmds_dot <- reactiveVal()

    plot_pca_dot <- reactiveVal()

    table_permanova <- reactiveVal()

    plot_real_heatmap <- reactiveVal()

}