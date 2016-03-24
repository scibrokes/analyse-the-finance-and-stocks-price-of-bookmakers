# server.R

suppressMessages(library('BBmisc'))
pkgs <- c('shiny', 'shinythemes', 'shinyjs', 'shinyBS', 'shinydashboard', 'shinyAce', 'quantmod', 
          'TTR', 'plyr', 'dplyr', 'stringr', 'purrr', 'googleCharts', 'lubridate', 'googleVis', 'zoo')
rm(pkgs)
suppressAll(source('helpers.R'))
tickers <- c('BET', 'BPTY', 'WMH', 'SPO', '888', 'PTEC', 'PCGE', 'TTR', 'GVC', 'WEB', 
             'BOX', 'RNK', 'LAD', 'STR') %>% sort

## =================================================================================================
shinyServer(function(input, output, session) {
  
  observe({
    if(input$selectall == 0) return(NULL)
    else if(input$selectall%%2 == 0){
      updateCheckboxGroupInput(session, 'counter', 'Select counter(s)', choices = tickers)
    } else{
      updateCheckboxGroupInput(session, 'counter', 'Select counter(s)', choices = tickers, selected = tickers)
    }
  })

  dataInput <- reactive({
    data <- llply(seq(input$counter), function(i) {
      getSymbols(Symbols = tickers[i], src = 'google', 
                 from = input$dates[1], to = input$dates[2], auto.assign = FALSE)})
    
    #'@ data01 <- ifelse(input$symb01==TRUE, getSymbols(Symbols = tickers[1], src = 'google', 
    #'@                        from = input$dates[1], to = input$dates[2], auto.assign = FALSE), NULL)
    #'@ data02 <- ifelse(input$symb02==TRUE, getSymbols(Symbols = tickers[2], src = 'google', 
    #'@                        from = input$dates[1], to = input$dates[2], auto.assign = FALSE), NULL)
    #'@ data03 <- ifelse(input$symb03==TRUE, getSymbols(Symbols = tickers[3], src = 'google', 
    #'@                        from = input$dates[1], to = input$dates[2], auto.assign = FALSE), NULL)
    #'@ data04 <- ifelse(input$symb04==TRUE, getSymbols(Symbols = tickers[4], src = 'google', 
    #'@                        from = input$dates[1], to = input$dates[2], auto.assign = FALSE), NULL)
    #'@ data05 <- ifelse(input$symb05==TRUE, getSymbols(Symbols = tickers[5], src = 'google', 
    #'@                        from = input$dates[1], to = input$dates[2], auto.assign = FALSE), NULL)
    #'@ data06 <- ifelse(input$symb06==TRUE, getSymbols(Symbols = tickers[6], src = 'google', 
    #'@                        from = input$dates[1], to = input$dates[2], auto.assign = FALSE), NULL)
    #'@ data07 <- ifelse(input$symb07==TRUE, getSymbols(Symbols = tickers[7], src = 'google', 
    #'@                        from = input$dates[1], to = input$dates[2], auto.assign = FALSE), NULL)
    #'@ data08 <- ifelse(input$symb08==TRUE, getSymbols(Symbols = tickers[8], src = 'google', 
    #'@                        from = input$dates[1], to = input$dates[2], auto.assign = FALSE), NULL)
    #'@ data09 <- ifelse(input$symb09==TRUE, getSymbols(Symbols = tickers[9], src = 'google', 
    #'@                        from = input$dates[1], to = input$dates[2], auto.assign = FALSE), NULL)
    #'@ data10 <- ifelse(input$symb10==TRUE, getSymbols(Symbols = tickers[10], src = 'google', 
    #'@                        from = input$dates[1], to = input$dates[2], auto.assign = FALSE), NULL)
    #'@ data11 <- ifelse(input$symb11==TRUE, getSymbols(Symbols = tickers[11], src = 'google', 
    #'@                        from = input$dates[1], to = input$dates[2], auto.assign = FALSE), NULL)
    #'@ data12 <- ifelse(input$symb12==TRUE, getSymbols(Symbols = tickers[12], src = 'google', 
    #'@                        from = input$dates[1], to = input$dates[2], auto.assign = FALSE), NULL)
    #'@ data13 <- ifelse(input$symb13==TRUE, getSymbols(Symbols = tickers[13], src = 'google', 
    #'@                        from = input$dates[1], to = input$dates[2], auto.assign = FALSE), NULL)
    #'@ data14 <- ifelse(input$symb14==TRUE, getSymbols(Symbols = tickers[14], src = 'google', 
    #'@                        from = input$dates[1], to = input$dates[2], auto.assign = FALSE), NULL)
    #'@ data15 <- ifelse(input$symb15==TRUE, getSymbols(Symbols = tickers[15], src = 'google', 
    #'@                        from = input$dates[1], to = input$dates[2], auto.assign = FALSE), NULL)
    #'@ data16 <- ifelse(input$symb16==TRUE, getSymbols(Symbols = tickers[16], src = 'google', 
    #'@                        from = input$dates[1], to = input$dates[2], auto.assign = FALSE), NULL)
    #'@ data17 <- ifelse(input$symb17==TRUE, getSymbols(Symbols = tickers[17], src = 'google', 
    #'@                        from = input$dates[1], to = input$dates[2], auto.assign = FALSE), NULL)
    #'@ data18 <- ifelse(input$symb18==TRUE, getSymbols(Symbols = tickers[18], src = 'google', 
    #'@                        from = input$dates[1], to = input$dates[2], auto.assign = FALSE), NULL)
    #'@ return (list(data01, data02, data03, data04, data05, data06, data07, data08, data09, 
    #'@              data10, data11, data12, data13, data14, data15, data16, data17, data18))
    return(data)
  })
  
  output$plot <- renderPlot({   
    data <- dataInput()
    if (input$adjust) data <- adjust(dataInput())
    
    llply(seq(dataInput()), function(i) {
      chartSeries(dataInput()[[i]], name=tickers[i], 
                  TA=c(addVo(), addBBands()), #paste0("addTA(", tickers[i], ",on=1)"), 
                  theme = chartTheme("white"), type = "line", log.scale = input$log)})
    
    ## added the TA argument with the paste helper function
    #'@ if(!is.null(dataInput()[[1]]))
    #'@   chartSeries(dataInput()[[1]], TA=paste0("addTA(", tickers[1], ",on=1)"), 
    #'@               theme = chartTheme("white"), type = "line", log.scale = input$log)
    #'@ if(!is.null(dataInput()[[2]]))
    #'@   chartSeries(dataInput()[[2]], TA=paste0("addTA(", tickers[2], ",on=1)"), 
    #'@               theme = chartTheme("white"), type = "line", log.scale = input$log)
    #'@ if(!is.null(dataInput()[[3]]))
    #'@   chartSeries(dataInput()[[3]], TA=paste0("addTA(", tickers[3], ",on=1)"), 
    #'@               theme = chartTheme("white"), type = "line", log.scale = input$log)
    #'@ if(!is.null(dataInput()[[4]]))
    #'@   chartSeries(dataInput()[[4]], TA=paste0("addTA(", tickers[4], ",on=1)"), 
    #'@               theme = chartTheme("white"), type = "line", log.scale = input$log)
    #'@ if(!is.null(dataInput()[[5]]))
    #'@   chartSeries(dataInput()[[5]], TA=paste0("addTA(", tickers[5], ",on=1)"), 
    #'@               theme = chartTheme("white"), type = "line", log.scale = input$log)
    #'@ if(!is.null(dataInput()[[6]]))
    #'@   chartSeries(dataInput()[[6]], TA=paste0("addTA(", tickers[6], ",on=1)"), 
    #'@               theme = chartTheme("white"), type = "line", log.scale = input$log)
    #'@ if(!is.null(dataInput()[[7]]))
    #'@   chartSeries(dataInput()[[7]], TA=paste0("addTA(", tickers[7], ",on=1)"), 
    #'@               theme = chartTheme("white"), type = "line", log.scale = input$log)
    #'@ if(!is.null(dataInput()[[8]]))
    #'@   chartSeries(dataInput()[[8]], TA=paste0("addTA(", tickers[8], ",on=1)"), 
    #'@               theme = chartTheme("white"), type = "line", log.scale = input$log)
    #'@ if(!is.null(dataInput()[[9]]))
    #'@   chartSeries(dataInput()[[9]], TA=paste0("addTA(", tickers[9], ",on=1)"), 
    #'@               theme = chartTheme("white"), type = "line", log.scale = input$log)
    #'@ if(!is.null(dataInput()[[10]]))
    #'@   chartSeries(dataInput()[[10]], TA=paste0("addTA(", tickers[10], ",on=1)"), 
    #'@               theme = chartTheme("white"), type = "line", log.scale = input$log)
    #'@ if(!is.null(dataInput()[[11]]))
    #'@   chartSeries(dataInput()[[11]], TA=paste0("addTA(", tickers[11], ",on=1)"), 
    #'@               theme = chartTheme("white"), type = "line", log.scale = input$log)
    #'@ if(!is.null(dataInput()[[12]]))
    #'@   chartSeries(dataInput()[[12]], TA=paste0("addTA(", tickers[12], ",on=1)"), 
    #'@               theme = chartTheme("white"), type = "line", log.scale = input$log)
    #'@ if(!is.null(dataInput()[[13]]))
    #'@   chartSeries(dataInput()[[13]], TA=paste0("addTA(", tickers[13], ",on=1)"), 
    #'@               theme = chartTheme("white"), type = "line", log.scale = input$log)
    #'@ if(!is.null(dataInput()[[14]]))
    #'@   chartSeries(dataInput()[[14]], TA=paste0("addTA(", tickers[14], ",on=1)"), 
    #'@               theme = chartTheme("white"), type = "line", log.scale = input$log)
    #'@ if(!is.null(dataInput()[[15]]))
    #'@   chartSeries(dataInput()[[15]], TA=paste0("addTA(", tickers[15], ",on=1)"), 
    #'@               theme = chartTheme("white"), type = "line", log.scale = input$log)
    #'@ if(!is.null(dataInput()[[16]]))
    #'@   chartSeries(dataInput()[[16]], TA=paste0("addTA(", tickers[16], ",on=1)"), 
    #'@               theme = chartTheme("white"), type = "line", log.scale = input$log)
    #'@ if(!is.null(dataInput()[[17]]))
    #'@   chartSeries(dataInput()[[17]], TA=paste0("addTA(", tickers[17], ",on=1)"), 
    #'@               theme = chartTheme("white"), type = "line", log.scale = input$log)
    #'@ if(!is.null(dataInput()[[18]]))
    #'@   chartSeries(dataInput()[[18]], TA=paste0("addTA(", tickers[18], ",on=1)"), 
    #'@               theme = chartTheme("white"), type = "line", log.scale = input$log)
    })
})

