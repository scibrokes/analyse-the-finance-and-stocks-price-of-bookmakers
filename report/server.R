## Loading multiple packages at once
pkgs <- c('plyr', 'dplyr', 'magrittr', 'shiny', 'DT', 'rCharts')
suppressMessages(l_ply(pkgs, require, character.only = TRUE, quietly = TRUE))
rm(pkgs)

options(RCHART_TEMPLATE = 'Rickshaw.html', RCHART_LIB = 'morris') #RCHART_LIB = 'polycharts'


server <- shinyServer(
  function(input, output, session){
    observe({
      if(input$selectall == 0) return(NULL)
      else if(input$selectall %% 2 == 0){
        updateCheckboxGroupInput(session, 'counter', 'Select counter(s)', choices = tickers)
      }else{
        updateCheckboxGroupInput(session, 'counter', 'Select counter(s)', choices = tickers, selected = tickers)
      }
    })
    
    ## Define a reactive expression for update data
    terms <- reactive({
      ## Change when the "Submit" button is pressed...
      input$submit
      
      ## ...but not for anything else
      isolate({
        withProgress({
          setProgress(message = 'Processing data...')
          filterData(input$counter)
          })
        })
      })
    
    output$table1 <- renderDataTable({
      dataTable <- terms()$dataTable
      dataTable <- filter(dataTable, Com %in% input$counter & Date>=input$dates[1] & Date<=input$dates[2])
      DT::datatable(dataTable,
                    caption='Table 1.1.1 : Public Listed Companies Stock Price',
                    extensions = list('ColReorder'=NULL, 'ColVis'=NULL, 'TableTools'=NULL
                                      #, 'FixedColumns'=list(leftColumns=2)
                    ), 
                    options = list(autoWidth=TRUE, 
                                   oColReorder=list(realtime=TRUE), #oColVis=list(exclude=c(0, 1), activate='mouseover'),
                                   oTableTools=list(
                                     sSwfPath='//cdnjs.cloudflare.com/ajax/libs/datatables-tabletools/2.1.5/swf/copy_csv_xls.swf',
                                     aButtons=list('copy', 'print',
                                                   list(sExtends='collection', sButtonText='Save', aButtons=c('csv','xls')))),
                                   dom='CRTrilftp', scrollX=TRUE, scrollCollapse=TRUE,
                                   colVis=list(exclude=c(0), activate='mouseover')))
    })
    
    ## http://blog.revolutionanalytics.com/2014/01/quantitative-finance-applications-in-r-plotting-xts-time-series.html
    output$plot <- renderPlot({
      dataPlot <- terms()$dataPlot
      
      llply(dataPlot, function(x) chartSeries(x, TA=c(addVo(),addBBands()), theme = chartTheme("white"),
                                              type = "line", log.scale = input$log))
      
    })
    
    ## comparison between rCharts and googleVis
    ## https://rpubs.com/miguelpatricio/r_charts_googleVis
    output$rChart <- renderPlot({
      
    })
    
    output$gvis1 <- renderGvis({
      dataPlot <- terms()$dataPlot
      dataPlot %<>% filter(Date>=input$dates[1] & Date<=input$dates[2])
      gvisCandlestickChart(smp, xvar='Date', low='Low', open='Open', close='Close', high='High', 
                           options=list(legend='none', gvis.editor='Edit me!'))
    })
    
    output$gvis2 <- renderGvis({
      ggPlot <- terms()$ggPlot
      ggPlot %<>% filter(Date>=input$dates[1] & Date<=input$dates[2])
      if(input$adjust) ggPlot <- adjust(ggPlot)
      if(input$log) ggPlot <- mutate(subset(ggPlot, Com %in% input$counter, 
                                            Date>=input$dates[1] & Date<=input$dates[2]), Close=log(Close))
      gvis.options <- list(title='Stock Closing Price Trend', series='[{targetAxisIndex:0},{targetAxisIndex:1}]', 
                           hAxis="{title:'Date'}', vAxis='{title:'Close'}", 
                           width='automatic', height='automatic',
                           titleTextStyle="{color:'red', fontName:'Courier', fontSize:16}",
                           curveType='function',  pointSize=9,
                           series="[{targetAxisIndex:0, color:'red'}, {targetAxisIndex:1, color:'blue'}]",
                           vAxes="[{title:'Percent', format:'#,###%',
                             titleTextStyle: {color: 'blue'},
                             textStyle:{color: 'blue'},
                             textPosition: 'out'}, 
                            {title:'Millions', format:'#,###',
                             titleTextStyle: {color: 'red'},  
                             textStyle:{color: 'red'},
                             textPosition: 'out'}]",
                           hAxes="[{title:'Date', textPosition: 'out'}]",
                           gvis.editor='Edit me!')
      gvisLineChart(xvar='Date', yvar=c('Price'), data=ggPlot, options=gvis.options)
    })
    
    ## http://rcharts.io/gallery/
    ## http://timelyportfolio.github.io/rCharts_rickshaw_gettingstarted/
    output$rChart1 <- renderChart2({
      ggPlot <- terms()$ggPlot
      ggPlot %<>% filter(Com==grep(input$counter, Com, value=TRUE), Date>=input$dates[1] & Date<=input$dates[2])
      r1 <- Rickshaw$new()
      r1$layer(y ~ x, groups='Com', data = ggPlot, type = "line", height = 240, width = 540)
      # add a helpful slider this easily; other features TRUE as a default
      r1$set(dom = 'plot', hoverDetail = TRUE, shelving = TRUE, legend = TRUE, slider = TRUE, highlight = TRUE)
      return(r1)
    })
})