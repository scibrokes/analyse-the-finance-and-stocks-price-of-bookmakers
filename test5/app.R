library(shiny)

spl <- function(s, delim = ',') {
  unlist(strsplit(s, delim))
}

createNonReactiveTextInput <- function(id = id, label = label, 
                                       value = strsplit(value, "(,\\s+)|,")[[1]], button.label = '') {
  value = lapply(substitute(value), function(x) as.character(x))
  
  if(button.label != '')
    list(
      tagList(
        tags$label(label), tags$input(id = id, type = "text", value = value, 
                                      style = "display:none;"), 
        tags$input(id = paste(id, "Temp", sep=''), type = "text", 
                   value = value, style = "display:inline;", 
                   onkeypress = paste("{if (event.keyCode==13) $('#", id, 
                                      "TempChange').click()}", sep = ''))), 
      div(
        tags$button(id = paste(id, "TempChange", sep=''), type = "button", 
                    class = "btn btn-primary", 
                    onclick = paste("$('#", id, "').val($('#", 
                                    id, "Temp').val()).change();", 
                                    sep = ''), button.label)))
  else
    list(
      tagList(
        tags$label(label), tags$input(id = id, type = "text", value = value, 
                                      style = "display:none;"), 
        tags$input(id = paste(id, "Temp", sep = ''), type = "text", 
                   value = value, style = "display:inline;", 
                   onkeypress = paste("{if (event.keyCode==13) $('#", id, 
                                      "').val($('#", id, 
                                      "Temp').val()).change()}", 
                                      sep = ''))))
}

createNonReactiveTextInputCustom <- function(id, label, tag.label = 'input', 
                                             button.label = '', 
                                             enableEnter = TRUE, opts) {
  onkeypress = ''
  if(button.label != '') {
    if(enableEnter)
      onkeypress = paste("{if (event.keyCode==13) $('#", id, 
                         "TempChange').click()}", sep = '')
    
    list(tagList(tags$label(label), 
                 tag(tag.label, c(id = id, style = "display:none;", opts)), 
                 tag(tag.label, c(id = paste(id, "Temp", sep = ''), 
                                  style = "display:inline;", 
                                  onkeypress = onkeypress, opts))),
         div(
           tags$button(id = paste(id, "TempChange", sep = ''), 
                       type = "button", class = "btn btn-primary", onclick = 
                         paste("$('#", id, "').val($('#", id, 
                               "Temp').val()).change();", sep = ''), 
                       button.label)))
  } else {
    if(enableEnter)
      onkeypress = paste("{if (event.keyCode==13) $('#", id, 
                         "').val($('#", id, 
                         "Temp').val()).change()}", sep = '')
    list(
      tagList(
        tags$label(label), tag(tag.label, c(id = id, style = 
                                              "display:none;", opts)), 
        tag(tag.label, c(id = paste(id,"Temp", sep = ''), 
                         style = "display:inline;", 
                         onkeypress = onkeypress, opts))))
  }
}

bt.prep <- function(b, align = c('keep.all', 'remove.na'), dates = NULL, 
                    fill.gaps = FALSE, basic = FALSE) {
  if(!exists('symbolnames', b, inherits = FALSE)) b$symbolnames = ls(b)
  symbolnames = b$symbolnames
  nsymbols = len(symbolnames)
  
  if(nsymbols > 1) {
    out = bt.merge(b, align, dates)
    for(i in 1:nsymbols) {
      temp = coredata(b[[symbolnames[i]]])[out$date.map[, i], , drop = FALSE]
      b[[symbolnames[i]]] = iif(basic, temp, make.xts(temp, out$all.dates))
      map.col = find.names('Close,Volume,Open,High,Low,Adjusted', 
                           b[[symbolnames[i]]])
      
      if(fill.gaps & !is.na(map.col$Close)) {
        close = coredata(b[[symbolnames[i]]][, map.col$Close])
        n = len(close)
        last.n = max(which(!is.na(close)))
        close = ifna.prev(close)
        
        if(last.n + 5 < n) close[last.n:n] = NA
        b[[symbolnames[i]]][, map.col$Close] = close
        index = !is.na(close)
        
        if(!is.na(map.col$Volume)) {
          index1 = is.na(b[[symbolnames[i]]][, map.col$Volume]) & index
          b[[symbolnames[i]]][index1, map.col$Volume] = 0
        }
        
        for(field in spl('Open,High,Low,Adjusted')) {
          j = map.col[[field]]
          if(!is.null(j)) {
            index1 = is.na(b[[symbolnames[i]]][,j]) & index
            b[[symbolnames[i]]][index1, j] = close[index1]
          }
        }
      }
    }
  } else {
    if(!is.null(dates)) b[[symbolnames[1]]] = b[[symbolnames[1]]][dates, ]
    out = list(all.dates = index.xts(b[[symbolnames[1]]]))
    
    if(basic) b[[symbolnames[1]]] = coredata(b[[symbolnames[1]]])
  }
  b$dates = out$all.dates
  dummy.mat = matrix(double(), len(out$all.dates), nsymbols)
  colnames(dummy.mat) = symbolnames
  
  if(!basic) dummy.mat = make.xts(dummy.mat, out$all.dates)
  b$weight = dummy.mat
  b$execution.price = dummy.mat
  for(i in 1:nsymbols) {
    if(has.Cl(b[[symbolnames[i]]])) {
      dummy.mat[, i] = Cl(b[[symbolnames[i]]]);
    }
  }
  b$prices = dummy.mat
}

bt.prep.matrix <- function(b, align = c('keep.all', 'remove.na'), 
                           dates = NULL, basic = FALSE) {
  align = align[1]
  nsymbols = len(b$symbolnames)
  if(!is.null(dates)) {
    temp = make.xts(1:len(b$dates), b$dates)
    temp = temp[dates]
    index = as.vector(temp)
    for(i in b$fields) b[[i]] = b[[i]][index, , drop = FALSE]
    b$dates = b$dates[index]
  }
  
  if(align == 'remove.na') {
    index = which(count(b$Cl, side=1) < nsymbols)
  } else {
    index = which(count(b$Cl, side = 1) < max(1, 0.1 * nsymbols))
  }
  
  if(len(index) > 0) {
    for(i in b$fields) b[[i]] = b[[i]][-index, , drop = FALSE]
    b$dates = b$dates[-index]
  }
  
  dummy.mat = iif(basic, b$Cl, make.xts(b$Cl, b$dates))
  b$weight = NA * dummy.mat
  b$execution.price = NA * dummy.mat
  b$prices = dummy.mat
}

bt.prep.matrix.test <- function() {
  load.packages('quantmod')
  returns = read.xts('Example.csv', date.fn=function(x) paste('1',x), 
                     format='%d %b-%y')
  prices = bt.apply.matrix(1 + returns, cumprod)
  data <- new.env()
  data$symbolnames = colnames(prices)
  data$dates = index(prices)
  data$fields = 'Cl'
  data$Cl = prices
  bt.prep.matrix(data)
  data$weight[] = NA
  data$weight[] = 1
  buy.hold = bt.run.share(data)
  plotbt(buy.hold, plotX = TRUE, log = 'y', LeftMargin = 3)
  mtext('Cumulative Performance', side = 2, line = 1)
}

bt.prep.remove.symbols.min.history <- function(b, min.history = 1000) {
  bt.prep.remove.symbols(b, which(count(b$prices, side=2) < min.history))
}

bt.prep.remove.symbols <- function(b, index) {
  if(len(index) > 0) {
    if(is.character(index)) index = match(index, b$symbolnames)
    b$prices = b$prices[, -index]
    b$weight = b$weight[, -index]
    b$execution.price = b$execution.price[, -index]
    env.rm(b$symbolnames[index], b)
    b$symbolnames = b$symbolnames[-index]
  }
}

bt.prep.trim <- function(b, dates = NULL) {
  if(is.null(dates)) return(b)
  dates.index = dates2index(b$prices, dates)
  data.copy <- new.env()
  
  for(s in b$symbolnames) data.copy[[s]] = b[[s]][dates.index, , drop = F]
  data.copy$symbolnames = b$symbolnames
  data.copy$dates = b$dates[dates.index]
  data.copy$prices = b$prices[dates.index, , drop = F]
  data.copy$weight = b$weight[dates.index, , drop = F]
  data.copy$execution.price = b$execution.price[dates.index, , drop = F]
  return(data.copy)
}

## ===================================================================

# Define UI for application that draws a histogram
ui <- shinyUI(fluidPage(
   
   # Application title
   titlePanel("Old Faithful Geyser Data"),
   
   # Sidebar with a slider input for number of bins 
   sidebarLayout(
      sidebarPanel(
        createNonReactiveTextInput(
          id = "symbols", 
          label = "Yahoo Ticker(s) separated by comma:", 
          value = "AAPL,GOOG", button.label = "Update"),
        br(),
        htmlOutput("status"),
        sliderInput("bins", "Number of bins:",
                    min = 1, max = 50, value = 30)),
      
      # Show a plot of the generated distribution
      mainPanel(
         plotOutput("distPlot")
      )
   )
))

# Define server logic required to draw a histogram
server <- shinyServer(function(input, output) {
  library('plyr')
  library('dplyr')
  library('magrittr')
  library('stringr')
  library('stringi')
  library('lubridate')
  
  # Create an environment for storing data
  symbol_env <- new.env()
  
  getData <- reactive(function() {  	
    cat('getData was called\n')
    data <- new.env()
    for(symbol in spl(toupper(input$symbols), delim = "(,\\s+)|,")) {
      #'@ if(is.null(symbol_env[[symbol]]))
      if(identical(symbol_env[[symbol]], emptyenv()))
      tryCatch({
        symbol_env[[symbol]] = 
          getSymbols(symbol, from = '1970-01-01', src = 'yahoo', 
                     auto.assign = FALSE)
      }, error = function(e) { 
        stop(paste('Problem getting prices for', symbol)) })
      data[[symbol]] = symbol_env[[symbol]]
      #'@ validate(need(input$symbols != "" & !is.null(input$symbols) & is.vector(input$symbols), label = "stock"))
    }
    bt.prep(data, align='keep.all')
    data
  })
  
   output$distPlot <- renderPlot({
      # generate bins based on input$bins from ui.R
      x    <- faithful[, 2] 
      bins <- seq(min(x), max(x), length.out = input$bins + 1)
      
      # draw the histogram with the specified number of bins
      hist(x, breaks = bins, col = 'darkgray', border = 'white')
   })
   
   output$status <- renderUI(function() {
     out = tryCatch(getData(), error=function(err) paste(err))
     if(is.character(out))
       HTML(paste("<b>Status</b>: <b><font color='red'>Error:</font></b>", out))
     else
       HTML("<b>Status</b>: <b><font color='green'>Ok</font></b>")
   })
})

# Run the application 
shinyApp(ui = ui, server = server)

