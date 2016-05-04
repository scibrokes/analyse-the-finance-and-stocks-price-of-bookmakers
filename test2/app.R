## Creating a parallel computing Cluster and support functions.
## Preparing the parallel cluster using the cores
doParallel::registerDoParallel(cores = 16)
pkgs <- c('shiny','shinydashboard','magrittr','plyr','dplyr','googleVis','quantmod','TTR','lubridate','DT','purrr')
suppressMessages(require('BBmisc'))
suppressMessages(lib(pkgs)); rm(pkgs)

## Set the googleVis options first to change the behaviour of plot.gvis, so that only the chart 
##  component of the HTML file is written into the output file.
op <- options(gvis.plot.tag='chart')

tickers <- sort(c('BET', 'BOTB', 'BPTY', 'WMH', 'SPO', '888', 'NPT', 'PTEC', 'PCGE', 'TTR', 'GVC', 'WEB', 
                  'BOX', 'RNK', 'LAD', 'STR', 'LVS', 'TNI', 'WYNN', 'MGM', 'MPEL', 'MCRI', 'GXYEY', 'BYD', 
                  'GIGNY', 'CACQ', 'PENN', 'GEBHY', 'WYNMF', 'CWLDY', 'SCHYF', 'SKYZF', 'MCHVY', 'EGT', 
                  'GPIC', 'SGMS', 'NGCRF', 'TACBY', 'IGT', 'GDEN', 'ISLE', 'GLXZ', 'SJMHF', 'PDSSF', 'PNK'))

ui <- dashboardPage(
  dashboardHeader(title = "Dynamic sidebar"),
  dashboardSidebar(),
  dashboardBody(
    fluidPage(
      titlePanel('Stock Price'),
      sidebarLayout(
        sidebarPanel(
          helpText('Select the stocks to examine. Information will be collected from Google/Yahoo finance.'),
          checkboxGroupInput('counter', 'Select counter(s)', choices = tickers), 
          actionLink('selectall', 'Select/Unselect All'),
          br(),
          br(),
          dateRangeInput('dates', 'Date range', start = '2013-01-01', end = as.character(Sys.Date())),
          checkboxInput('log', 'Plot y axis on log scale', value = FALSE),
          checkboxInput('adjust', 'Adjust prices for inflation', value = FALSE)),
        mainPanel(
          tabsetPanel( #tabBox
            tabPanel('Table', dataTableOutput('table1'), 
                     downloadButton('downloadData', 'Download')),
            tabPanel('Chart', htmlOutput('gvis')))
        )))
  )
)

server <- function(input, output, session) {
  observe({
    if(input$selectall == 0) return(NULL)
    else if(input$selectall%%2 == 0){
      updateCheckboxGroupInput(session, 'counter', 'Select counter(s)', choices = tickers)
    }else{
      updateCheckboxGroupInput(session, 'counter', 'Select counter(s)', choices = tickers, selected = tickers)
    }
  })
  
  getCounter <- function(){
    data <- na.omit(data.frame(Com=NA, Open=NA, High=NA, Low=NA, Close=NA, Volume=NA, Date=NA, Weekday=NA))
    data <- llply(input$counter, function(x) {
      y = data.frame(Com=x, as.data.frame(
        tryCatch(suppressAll(na.omit(getSymbols(x, src = 'google', auto.assign = FALSE, 
                                        start = input$dates[1], end = input$dates[2]))), 
                 error = function(e) data.frame(Open=NA, High=NA, Low=NA, Close=NA, Volume=NA)))); 
      if(nrow(na.omit(y))==0) {
        y = data.frame(Com=x, as.data.frame(
          tryCatch(suppressAll(na.omit(getSymbols(x, src = 'yahoo', auto.assign = FALSE, 
                                          start = input$dates[1], end = input$dates[2]))), 
                   error = function(e) data.frame(Open=NA, High=NA, Low=NA, Close=NA, Volume=NA))))
      }else{
        y %<>% mutate(Date=as.character(rownames(.)), Weekday=factor(weekdays(ymd(Date))));
        names(y) <- c('Com', 'Open', 'High', 'Low', 'Close', 'Volume', 'Date', 'Weekday')}; 
      y}, .parallel=TRUE) %>% rbind_all %>% tbl_df %>% 
      select(Date, Weekday, Com, Open, High, Low, Close, Volume) %>% mutate(Com=factor(Com)) %>% na.omit
    return(data)
  }
  
  output$table1 <- renderDataTable({
    data <- getCounter()
    data <- datatable(data, 
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
    data <- getCounter()
    if(input$adjust) data <- adjust(data)
    if(input$log) data <- mutate(data, Low=log(Low), Open=log(Open), Close=log(Close), High=log(High))
    gvisCandlestickChart(data, xvar='Date', low='Low', open='Open', close='Close', high='High', 
                         options=list(legend='none', gvis.editor="Edit me!"))
  })
}

#'@ runApp(shinyApp(ui, server), display.mode='showcase')
shinyApp(ui, server)
