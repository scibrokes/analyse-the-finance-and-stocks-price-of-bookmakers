## http://stackoverflow.com/questions/31641704/plotting-multiple-symbols-with-a-reactive-statement-with-chartseries

library(quantmod)

shinyServer(function(input, output) {
  
  dataInput <- reactive({   
    data1 <- getSymbols(input$symb1, src = "google",      #Seperated the data into two seperate data sets and set auto.assign=FALSE
                        from = input$dates[1],
                        to = input$dates[2],
                        auto.assign = FALSE)
    data2 <- getSymbols(input$symb2, src = "google",     #Seperated the data into two seperate data sets and set auto.assign=FALSE
                        from = input$dates[1],
                        to = input$dates[2],
                        auto.assign = FALSE)
    return (list(data1,data2))                          #Stored the data sets in a single list 
  })
  
  output$plot <- renderPlot({
    chartSeries(dataInput()[[1]], TA=paste0("addTA(",input$symb1,",on=1)"),theme = chartTheme("white"),    #added the TA argument with the paste helper function
                type = "line", log.scale = input$log)
    chartSeries(dataInput()[[2]], TA=paste0("addTA(",input$symb2,",on=1)"),theme = chartTheme("white"), 
                type = "line", log.scale = input$log)
  }) 
})