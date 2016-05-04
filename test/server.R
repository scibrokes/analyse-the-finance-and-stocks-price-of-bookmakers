## http://stackoverflow.com/questions/31641704/plotting-multiple-symbols-with-a-reactive-statement-with-chartseries

library(quantmod)
library(magrittr)
library(dplyr)

shinyServer(function(input, output) {
  
  dataInput <- reactive({   
    data <- llply(c(input$symb1, input$symb2), function(x){
      getSymbols(x, src = "google", from = input$dates[1], to = input$dates[2], auto.assign = FALSE)})
    names(data) <- c(input$symb1, input$symb2)
    return (data)
  })
  
  output$plot <- renderPlot({
    data <- dataInput()
    names(data) <- c(input$symb1, input$symb2)
    llply(data, function(x) chartSeries(x, TA=c(addVo(),addBBands()), theme = chartTheme("white"),
                type = "line", log.scale = input$log))
  })
})