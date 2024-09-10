library(shiny)
library(shinyjs)

# UI
ui <- fluidPage(
  useShinyjs(),  # Initialize shinyjs
  actionButton("start_analysis", "Start Analysis"),
  textOutput("status")
)

# Server
server <- function(input, output, session) {
  is_running <- reactiveVal(FALSE)    # Reactive value to track if the process is running
  initial_delay_done <- reactiveVal(FALSE)  # Reactive value to track if the initial delay is done
  
  # Function to check for new files in the target directory
  check_files <- function() {
    # Code to check for files in the directory
    print(paste("Checking files at", Sys.time()))
  }
  
  # Start analysis when button is clicked
  observeEvent(input$start_analysis, {
    is_running(!is_running())
    updateActionButton(session, "start_analysis", 
                       label = if (is_running()) "Stop Analysis" else "Start Analysis")
    
    if (is_running()) {
      # Reset the initial delay status
      initial_delay_done(FALSE)
      
      print(paste0("Initial Delay Started at", Sys.time()))
      
      # Use a delayed execution using shinyjs::delay for the initial 2-minute delay
      delay(3000, {  # 2-minute delay
        initial_delay_done(TRUE)  # Mark initial delay as done
      })
    }
  })
  
  # Observe block to handle file checking
  observe({
    # Only proceed if the process is running and initial delay is done
    req(is_running(), initial_delay_done())
    
    # Perform the file check
    check_files()
    
    # Schedule the next check after 2 minutes
    invalidateLater(30000, session)
  })
}

# Run the application
shinyApp(ui = ui, server = server)