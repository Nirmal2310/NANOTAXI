#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

source("packages.R")

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
  source("ui-tab-intro.R", local = TRUE)$value,
  source("ui-tab-input.R", local = TRUE)$value,
  source("ui-tab-cohort-analysis.R", local = TRUE)$value,
  source("ui-tab-help.R", local = TRUE)$value,
  source("ui-tab-conditions.R", local = TRUE)$value
 )