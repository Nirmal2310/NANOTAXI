options(shiny.maxRequestSize = 100*1024^2)
source("packages.R")
print(sessionInfo())
server <- function(input, output, session) {
    source("server-input.R", local = TRUE)
    source("server-cohort-analysis.R", local = TRUE)
}