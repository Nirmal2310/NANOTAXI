observeEvent(input$StartAnalysis, {
        sample_data <- input_data_reactive()
        list <- sample_data$data[-1, 1]
        ref <- input$fastafile
        threads <- input$threads
        memory <- input$memory
        comp <- input$comp
        cont <- input$cont
    if (input$Setup) {
        system("bash Installation/env_install.sh")
        system("while read sample; \
        do bash Pipeline/rgi_main.sh -s \ 
        $sample -r $ref -t $threads -m $memory \
        -c $comp -d $cont; done < list")
    }else {
        system("while read sample; \
        do bash Pipeline/rgi_main.sh -s \ 
        $sample -r $ref -t $threads -m $memory \
        -c $comp -d $cont; done < list")
    }
    dir.create("Results")
    system("mv *_out/*_consolidated_final_arg_counts.txt Results/")
    system("mv *_out/*_family_info.txt Results/")
})
