rm(list = ls())
library(shiny)
library(xts)
library(DT)
library(quantmod)
library(shinysky)
library(PerformanceAnalytics)

x <- list()
Data <- data.frame()
clo <- data.frame()

ui <- fluidPage(
  titlePanel("test"),
  tags$style(type="text/css", ".shiny-output-error { visibility: hidden; }", 
             ".shiny-output-error:before { visibility: hidden; }"),
  sidebarLayout(
    sidebarPanel(
      helpText("text"),
      select2Input("txt", "stock", choices = c("AAPL", "GOOG", "MSFT"), 
                   selected = c("AAPL", "GOOG")),
      dateInput("dates", "Date:", value = "2016-01-01"),
      actionButton("go", "submit") 
    ),
    mainPanel(
      tabsetPanel(type = "tab", 
                  tabPanel("Plot", plotOutput("plot")), 
                  tabPanel("summary", dataTableOutput("table")),
                  tabPanel("close", dataTableOutput("table1")))
    )
  )
)

server <- function(input, output) {
  
  stockData <- new.env()
  dataInput <- reactive({ 
    if(input$go == 0){return()} #confirming button click
    isolate({
      input$go
      getSymbols(input$txt, src = "yahoo", env = stockData, 
                 from = as.Date(input$dates))    
      Data <- data.frame()
      
      validate(need(input$txt != "" & !is.null(input$txt) & is.vector(input$txt), label = "stock"))
      for (i in 1:length(input$txt)) {
        x[[i]] <- get(input$txt[i], pos = stockData)  # get data from stockData environment  
        Data <- cbind(Data, diff(log(Cl(x[[i]]))))
      }
      Data
    })
  })
  
  Last_Close <- reactive({
    if(input$go==0){return()} #confirming button click
    isolate({
      input$go
      validate(need(input$txt != "" & !is.null(input$txt) & is.vector(input$txt), label = "stock"))
      for (i in 1:length(input$txt)) {
        x[[i]] <- get(input$txt[i], pos=stockData)  # get data from stockData environment  
        clo <- cbind(clo,Cl(x[[i]]))
      }
      clo
    })
  })
  
  output$plot <- renderPlot(chart.Correlation(dataInput()))
  output$table <- DT::renderDataTable(datatable(as.data.frame(Cl(dataInput()))))
  output$table1 <- DT::renderDataTable(datatable(as.data.frame(Cl(Last_Close()))))
}

shinyApp(ui, server)