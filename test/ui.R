## http://stackoverflow.com/questions/31641704/plotting-multiple-symbols-with-a-reactive-statement-with-chartseries

library(shiny)

shinyUI(fluidPage(
  titlePanel("StockComp"),
  
  sidebarLayout(
    sidebarPanel(
      helpText("Select two stocks and a time frame to compare. 
        Information will be collected from google finance."),
      
      textInput("symb1", "1st Stock Symbol", "WMH"),
      textInput("symb2", "2nd Stock Symbol", "LAD"),
      
      dateRangeInput("dates", 
                     "Date range",
                     start = "2013-01-01", 
                     end = as.character(Sys.Date())),
      
      actionButton("get", "Compare Stocks"),
      
      br(),
      br(),
      
      checkboxInput("log", "Plot y axis on log scale", 
                    value = FALSE)
      
    ),
    
    mainPanel(plotOutput("plot"))
  )
))