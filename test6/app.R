suppressMessages(library('shiny'))
suppressMessages(library('quantmod'))
suppressMessages(library('plyr'))
suppressMessages(library('dplyr'))
suppressMessages(library('magrittr'))
suppressMessages(library('stringr'))
suppressMessages(library('stringi'))
suppressMessages(library('lubridate'))

# Download data for a stock if needed, and return the data
require_symbol <- function(symbol, envir = parent.frame()) {
  if(is.null(envir[[symbol]]))
    envir[[symbol]] = getSymbols(symbol, auto.assign = FALSE)
  envir[[symbol]]
}

# Define UI for application that draws a histogram
ui <- shinyUI(fluidPage(
  
  # Application title
  titlePanel("Old Faithful Geyser Data"),
  
  # Sidebar with a slider input for number of bins 
  sidebarLayout(
    sidebarPanel(
      textInput(inputId = "stock1", label = "Counter 01"), 
      selectInput(inputId = "chart_type", label = "Graph", 
                  choices = c("Candle-Sticks" = "candlesticks", 
                              "Match-Sticks" = "matchsticks", 
                              "Bar-Charts" = "bars", "Line-Charts" = "line")),
      sliderInput(inputId = "time_num", label = "Number of Time Interval", 
                  min = 1, max = 24, step = 1, value = 6), 
      selectInput(inputId = "time_unit", 
                  label = "Unit of Time Interval", 
                  choices = c("Day" = "days", "Week" = "weeks", 
                              "Month" = "months", "Year" = "years"), 
                  selected = "Months")),
    # Show a plot of the generated distribution
    mainPanel(
      conditionalPanel(condition = "input$stock1", br(), 
                       div(plotOutput(outputId = "plot1")))
      )
    )
  )
)

# Define server logic required to draw a histogram
server <- shinyServer(function(input, output) {
  # Create an environment for storing data
  symbol_env <- new.env()
  
  # Make a chart for a symbol, with the settings from the inputs
  make_chart <- function(symbol) {
    symbol_data <- require_symbol(symbol, symbol_env)
    #TA_STR <- paste0()
    chartSeries(symbol_data,
                name       = symbol,
                type       = input$chart_type,
                subset     = paste("last", input$time_num, input$time_unit),
                theme      = "white")
  }
  
  output$plot1 <- renderPlot( make_chart(input$stock1))
})

# Run the application 
shinyApp(ui = ui, server = server)

