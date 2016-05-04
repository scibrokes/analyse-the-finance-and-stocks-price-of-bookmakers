# server.R

suppressMessages(library('BBmisc'))
pkgs <- c('shiny', 'shinythemes', 'shinyjs', 'shinyBS', 'shinydashboard', 'shinyAce', 'quantmod', 
          'TTR', 'plyr', 'dplyr', 'magrittr', 'stringr', 'purrr', 'googleCharts', 'lubridate', 
          'googleVis', 'DT', 'sparkline', 'zoo')
suppressAll(source('helpers.R'))
suppressMessages(lib(pkgs))
rm(pkgs)

tickers <- sort(c('BET', 'BOTB', 'BPTY', 'WMH', 'SPO', '888', 'NPT', 'PTEC', 'PCGE', 'TTR', 'GVC', 'WEB', 
                  'BOX', 'RNK', 'LAD', 'STR', 'LVS', 'TNI', 'WYNN', 'MGM', 'MPEL', 'MCRI', 'GXYEY', 'BYD', 
                  'GIGNY', 'CACQ', 'PENN', 'GEBHY', 'WYNMF', 'CWLDY', 'SCHYF', 'SKYZF', 'MCHVY', 'EGT', 
                  'GPIC', 'SGMS', 'NGCRF', 'TACBY', 'IGT', 'GDEN', 'ISLE', 'GLXZ', 'SJMHF', 'PDSSF', 'PNK'))

## =================================================================================================
server <- shinyServer(
  function(input, output, session){
    observe({
      if(input$selectall == 0){
        return(NULL)
      } else if(input$selectall%%2 == 0){
        updateCheckboxGroupInput(session, 'counter', 'Select counter(s)', choices = tickers)
      } else{
        updateCheckboxGroupInput(session, 'counter', 'Select counter(s)', choices = tickers, selected = tickers)
      }
      
      if(input$counter>0 & input$counter<=length(input$counter)){
        isolate({
          values<-renderText({
            input$counter
          })
          dataInput()
        })
      }
    })
    
    ## http://stackoverflow.com/questions/29282524/obtain-the-values-from-a-checkboxgroup-only-once-in-shiny
    dataInput <- reactive({
      data <- llply(input$counter, function(x){
        y = data.frame(Com=x, as.data.frame(
          tryCatch(suppressAll(getSymbols(x, src = 'google', auto.assign = FALSE, 
                                          start = input$dates[1], end = input$dates[2]
          )), error = function(e) 
            data.frame(Open=NA, High=NA, Low=NA, Close=NA, Volume=NA)))); 
        if(nrow(na.omit(y))==0){
          y = data.frame(Com=x, as.data.frame(
            tryCatch(suppressAll(getSymbols(x, src = 'yahoo', auto.assign = FALSE, start = start, end = end)), 
                     error = function(e) data.frame(Open=NA, High=NA, Low=NA, Close=NA, Volume=NA))))
        }else{
          y %<>% mutate(Date=as.character(rownames(.)), Weekday=factor(weekdays(ymd(Date))));
          names(y) <- c('Com', 'Open', 'High', 'Low', 'Close', 'Volume', 'Date', 'Weekday')}; 
        y}, .parallel=TRUE) %>% rbind_all %>% tbl_df %>% 
        select(Date, Weekday, Com, Open, High, Low, Close, Volume) %>% mutate(Com=factor(Com))
      names(data) <- input$counter$value
      return (data)
    })
    
    output$table1 <- renderDataTable({
      data <- dataInput()
      dt <- datatable(data, 
                      caption="Table 1.1.1 : Public Listed Companies Stock Price",
                      extensions = list("ColReorder"=NULL, "ColVis"=NULL, "TableTools"=NULL
                                        #, "FixedColumns"=list(leftColumns=2)
                      ), 
                      options = list(autoWidth=TRUE,
                                     oColReorder=list(realtime=TRUE), #oColVis=list(exclude=c(0, 1), activate='mouseover'),
                                     oTableTools=list(
                                       sSwfPath="//cdnjs.cloudflare.com/ajax/libs/datatables-tabletools/2.1.5/swf/copy_csv_xls.swf",
                                       aButtons=list("copy", "print",
                                                     list(sExtends="collection", sButtonText="Save",
                                                          aButtons=c("csv","xls")))),
                                     dom='CRTrilftp', scrollX=TRUE, scrollCollapse=TRUE,
                                     colVis=list(exclude=c(0), activate='mouseover')))
      return(dt)
    })
    
    output$gvis <- renderGvis({
      data <- dataInput()
      #'@ smp <- subset(data, Com %in% input$counter &
      #'@               Date>=input$dates[1] & Date<=input$dates[2])
      if(input$adjust) data <- adjust(data)   #'@ smp <- adjust(subset(data, Com %in% input$counter &
      #'@ Date>=input$dates[1] & Date<=input$dates[2]))
      if(input$log) #'@ smp <- mutate(subset(data, Com %in% input$counter, 
      #'@                      Date>=input$dates[1] & Date<=input$dates[2]), 
      data <- mutate(Low=log(Low), Open=log(Open), Close=log(Close), High=log(High))
      gvisCandlestickChart(data, xvar='Date', low='Low', open='Open', close='Close', high='High', 
                           options=list(legend='none', gvis.editor="Edit me!"))
    })
    
    datasetInput <- reactive({
      smp <- dataInput()
      switch(smp,
             "table" = smp,
             "chart" = plot,
             "comp" = comp)})
    
    output$table <- renderTable({
      datasetInput()
    })
    
    output$downloadData <- downloadHandler(
      data <- dataInput(),
      filename = function() { 
        paste(input$dataset, '.csv', sep='') 
      },
      content = function(file) {
        write.csv(datasetInput(), file)
      })
    
    output$progressBox <- renderInfoBox({
      infoBox(
        "Progress", paste0(25 + input$count, "%"), icon = icon("list"),
        color = "purple"
      )
    })
    output$approvalBox <- renderInfoBox({
      infoBox(
        "Approval", "80%", icon = icon("thumbs-up", lib = "glyphicon"),
        color = "yellow"
      )
    })
    
    # Same as above, but with fill=TRUE
    output$progressBox2 <- renderInfoBox({
      infoBox(
        "Progress", paste0(25 + input$count, "%"), icon = icon("list"),
        color = "purple", fill = TRUE
      )
    })
    output$approvalBox2 <- renderInfoBox({
      infoBox(
        "Approval", "80%", icon = icon("thumbs-up", lib = "glyphicon"),
        color = "yellow", fill = TRUE
      )
    })
  })
