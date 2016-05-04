# server.R

suppressMessages(library('BBmisc'))
pkgs <- c('shiny', 'shinythemes', 'shinyjs', 'shinyBS', 'shinydashboard', 'shinyAce', 'quantmod', 
          'TTR', 'plyr', 'dplyr', 'magrittr', 'stringr', 'purrr', 'googleCharts', 'lubridate', 
          'googleVis', 'DT', 'sparkline', 'zoo')
suppressMessages(lib(pkgs))
rm(pkgs)

tickers <- c('BET', 'BOTB', 'BPTY', 'WMH', 'SPO', '888', 'NPT', 'PTEC', 'PCGE', 'TTR', 'GVC', 'WEB', 
             'BOX', 'RNK', 'LAD', 'STR', 'LVS', 'TNI', 'WYNN', 'MGM', 'MPEL', 'MCRI', 'GXYEY', 'BYD', 
             'GIGNY', 'CACQ', 'PENN', 'GEBHY', 'WYNMF', 'CWLDY', 'SCHYF', 'SKYZF', 'MCHVY', 'EGT', 
             'GPIC', 'SGMS', 'NGCRF', 'TACBY', 'IGT', 'GDEN', 'ISLE', 'GLXZ', 'SJMHF', 'PDSSF', 'PNK') %>% sort

data <- llply(tickers, function(x) {
  y = data.frame(Com=x, as.data.frame(
    tryCatch(suppressAll(getSymbols(x, src = 'google', auto.assign = FALSE)), error = function(e) 
      data.frame(Open=NA, High=NA, Low=NA, Close=NA, Volume=NA)))); 
  if(nrow(na.omit(y))==0) {
    y = data.frame(Com=x, as.data.frame(
      tryCatch(suppressAll(getSymbols(x, src = 'yahoo', auto.assign = FALSE)), 
               error = function(e) data.frame(Open=NA, High=NA, Low=NA, Close=NA, Volume=NA))))
  }else{
    y %<>% mutate(Date=as.character(rownames(.)), Weekday=factor(weekdays(ymd(Date))));
    names(y) <- c('Com', 'Open', 'High', 'Low', 'Close', 'Volume', 'Date', 'Weekday')}; 
  y}, .parallel=TRUE) %>% rbind_all %>% tbl_df %>% 
  select(Date, Weekday, Com, Open, High, Low, Close, Volume) %>% mutate(Com=factor(Com))

#'@ googleChartsInit(chartTypes = c("ALL", "annotatedtimeline", "area", "bar", "bubble", "calendar", "candlestick", 
#'  "column", "combo", "gauge", "geo", "geomap", "intensitymap", "line", "map", "motion", "org", "pie", "sankey", 
#'  "scatter", "steppedarea", "table", "timeline", "treemap"))
## =================================================================================================

#'@ shiny::runApp(list(
  ui = pageWithSidebar(
    headerPanel("Example"),
    sidebarPanel(
      helpText('Select the stocks to examine. Information will be collected from Google/Yahoo finance.'),
      checkboxGroupInput('counter', 'Select counter(s)', choices = tickers), 
      actionLink('selectall', 'Select/Unselect All'),
      br(),
      br(),
      dateRangeInput('dates', 'Date range', start = '2013-01-01', end = as.character(Sys.Date())),
      br(),
      checkboxInput('log', 'Plot y axis on log scale', value = FALSE),
      checkboxInput('adjust', 'Adjust prices for inflation', value = FALSE)),
    mainPanel(
      tabsetPanel(
        tabPanel('Data', dataTableOutput('table1'),
                 downloadButton('downloadData', 'Download')),
        tabPanel('Plot', htmlOutput('gvis')))
    ))#'@ ,

  server = function(input, output, session) {
    observe({
      if(input$selectall == 0) return(NULL)
      else if(input$selectall%%2 == 0){
        updateCheckboxGroupInput(session, 'counter', 'Select counter(s)', choices = tickers)
      }else{
        updateCheckboxGroupInput(session, 'counter', 'Select counter(s)', choices = tickers, selected = tickers)
      }
    })
    
    output$table1 <- renderDataTable({
      datatable(subset(data, Com %in% input$counter & Date>=input$dates[1] & Date<=input$dates[2]),
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
    })
    
    output$gvis <- renderGvis({
      smp <- subset(data, Com %in% input$counter &
                      Date>=input$dates[1] & Date<=input$dates[2])
      if(input$adjust) smp <- adjust(subset(data, Com %in% input$counter &
                                            Date>=input$dates[1] & Date<=input$dates[2]))
      if(input$log) smp <- mutate(subset(data, Com %in% input$counter, 
                                         Date>=input$dates[1] & Date<=input$dates[2]), 
                                  Low=log(Low), Open=log(Open), Close=log(Close), High=log(High))
      gvisCandlestickChart(smp, xvar='Date', low='Low', open='Open', close='Close', high='High', 
                           options=list(legend='none', gvis.editor="Edit me!"))
    })
  }
#'@ ))

# Run the application 
shinyApp(ui = ui, server = server)
