# Define server
shinyServer(function(input, output) {
  library('plyr')
  library('dplyr')
  library('magrittr')
  library('stringr')
  library('stringi')
  library('lubridate')
  
  # Create an environment for storing data
  symbol_env <- new.env()
  
  #*****************************************************************
  # Shared Reactive functions
  # http://rstudio.github.com/shiny/tutorial/#inputs-and-outputs
  #*****************************************************************    	
  # Get stock data
  getData <- reactive(function() {  	
    cat('getData was called\n')
    data <- new.env()
    for(symbol in spl(toupper(input$symbols), delim = "(,\\s+)|,")) {
      #'@ if(is.null(symbol_env[[symbol]]))
      if(identical(symbol_env[[symbol]], emptyenv()))
      ## http://shiny.rstudio.com/articles/validation.html
      #'@ validate(need((symbol_env[[symbol]] == "" | 
      #'@                is.null(symbol_env[[symbol]])) & 
      #'@               !is.vector(symbol_env[[symbol]]), label = "stock"))
        tryCatch({
          symbol_env[[symbol]] = 
            getSymbols(symbol, from = '1970-01-01', src = 'yahoo', 
                       auto.assign = FALSE)
        }, error = function(e) { 
          stop(paste('Problem getting prices for', symbol)) })
      
      #'@ symbol_env[[symbol]] = as.data.frame(
      #'@  tryCatch(suppressAll(getSymbols(symbol, src = 'yahoo', auto.assign = FALSE, from=from, to=to)), 
      #'@           error = function(e) data.frame(Open=NA, High=NA, Low=NA, Close=NA, Volume=NA, Adjusted=NA)))
      #'@ if(nrow(na.omit(symbol_env[[symbol]]))==0){
      #'@  symbol_env[[symbol]] = as.data.frame(
      #'@    tryCatch(suppressAll(getSymbols(symbol, src = 'google', auto.assign = FALSE, 
      #'@                                    from=from, to=to)), error = function(e) 
      #'@                                      data.frame(Open=NA, High=NA, Low=NA, Close=NA, Volume=NA, Adjusted=NA)))
      #'@  symbol_env[[symbol]]$Adjusted <- as.numeric(NA)
      #'@  names(symbol_env[[symbol]])[length(symbol_env[[symbol]])] <- paste0(symbol, '.Adjusted')
      #'@  names(symbol_env[[symbol]]) <- gsub('X888', '888', names(symbol_env[[symbol]]))
      #'@ }
      #'@
      #'@ z = symbol_env[[symbol]] %>% mutate(Date=as.character(rownames(.)), Weekday=factor(weekdays(ymd(Date)))) %>% 
      #'@   select(Date, Weekday)
      #'@ symbol_env[[symbol]] = data.frame(z, symbol_env[[symbol]]) %>% tbl_df
      #'@ rm(z)
      data[[symbol]] = symbol_env[[symbol]]
    }
    bt.prep(data, align='keep.all')
    data
  })
  
  # Determine dates range
  getDateRange <- reactive(function() {  	
    data = getData()
    max(1, nrow(data$prices) - input$dateRange):nrow(data$prices)
  })
  
  # Make table
  makeCorTable <- reactive(function() {  	
    out = getData()
    
    prices = out$prices[getDateRange(), ]
    if( ncol(prices) == 1) return(NULL)
    
    # compute correlation
    ret = prices / mlag(prices) - 1
    100 * cor(coredata(ret), use = 'complete.obs', method = 'pearson')	
  })
  
  #*****************************************************************
  # Not Reactive helper functions
  #*****************************************************************
  # Make stock plot
  makeStockPlot <- function() {  	
    out = getData()
    
    prices = out$prices[getDateRange(), ]
    if(input$plotReturnsFlag) prices = scale.one(prices)
    
    plota.matplot(prices)
    plota.add.copyright()
  }
  
  #*****************************************************************
  # Update plot(s) and table(s)
  #******************************************************************    	
  # Generate a plot
  output$stockPlot <- renderPlot(function() {
    makeStockPlot()
  }, height = 400, width = 600)
  
  # Generate a table
  output$corTable <- reactive(function() {
    temp = makeCorTable()	
    tableColor(as.matrix(temp), digits = 1)		
  })
  
  #*****************************************************************
  # Download
  #******************************************************************    
  # Download pdf report
  output$downloadReport <- downloadHandler(
    filename = 'report.pdf',
    content = function(file) {
      pdf(file = file, width = 8.5, height = 11)
      
      makeStockPlot()
      plot.table(round(makeCorTable(), 1), 'Correlation')
      plota.add.copyright()
      
      dev.off()
    }
  )
  
  # Download csv data
  output$downloadData <- downloadHandler(
    filename = 'data.csv',
    content = function(file) {
      write.csv(makeCorTable(), file)
    }
  )	
  
  #*****************************************************************
  # Update status message 
  #******************************************************************    
  output$status <- renderUI(function() {
    out = tryCatch(getData(), error=function(err) paste(err))
    if(is.character(out))
      HTML(paste("<b>Status</b>: <b><font color='red'>Error:</font></b>", out))
    else
      HTML("<b>Status</b>: <b><font color='green'>Ok</font></b>")
  })
})


