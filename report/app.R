#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

suppressMessages(library('BBmisc'))
pkgs <- c('shiny', 'shinythemes', 'shinyjs', 'shinyBS', 'shinydashboard', 'shinyAce', 'quantmod', 
          'TTR', 'plyr', 'dplyr', 'stringr', 'purrr', 'googleCharts')
suppressMessages(lib(pkgs))
suppressMessages(source('stocks.R'))
rm(pkgs)

## =================================================================================
messageData <- data_frame(from=c('Andy', 'Terry', 'Jimmy'), 
                          message=c('Your song: 世界第一等 & 欢聚歌', 'Gaming introduction',
                                    'World Vision Donation'))

dbHeader <- dashboardHeader(title = 'Reporting Dashboard',
                dropdownMenuOutput('messageMenu'),
                dropdownMenu(type = 'notifications',
                             notificationItem(text = '5 new users today', icon('users')),
                             notificationItem(text = '12 items delivered', 
                                              icon('truck'), status = 'success'),
                             notificationItem(text = 'Server load at 86%', 
                                              icon = icon('exclamation-triangle'), 
                                              status = 'warning')),
                dropdownMenu(type = 'tasks',
                             badgeStatus = 'success',
                             taskItem(value = 90, color = 'green', 'Documentation'),
                             taskItem(value = 17, color = 'aqua', 'Project X'),
                             taskItem(value = 75, color = 'yellow', 'Server deployment'),
                             taskItem(value = 80, color = 'red', 'Overall project')))
dbHeader$children$children <- HTML("<a href='http://www.scibrokes.com' target='_blank'>
                        <img align='right' alt='Logo' src='www/oda-army.jpg'/></a>")

## ---------------------------------------------------------------------------------
ui <- dashboardPage(
  dbHeader,
  
  ## Sidebar content
  dashboardSidebar(
    sidebarMenu(
      menuItem('Dashboard', tabName = 'dashboard', icon = icon('dashboard')),
      menuItem('Widgets', tabName = 'widgets', icon = icon('th')),
      menuItem('Stocks Prices', tabName = 'line-charts', icon = icon('line-chart')),
      menuItem('Financial Statements', tabName = 'bar-charts', icon = icon('bar-chart'),
      menuSubItem('Annual Summary Charts', tabName='bar-charts', icon = icon('bar-chart'))),
      menuItem('Calendar', tabName = 'calendar', icon = icon('calendar')),
      menuItem('References', tabName = 'file-pdf', icon = icon('file-pdf-o')),
      menuItem('Source code', icon = icon('file-code-o'), 
               href = 'https://github.com/scibrokes/analyse-the-finance-and-stocks-price-of-bookmakers')
    )),
  ## Body content
  dashboardBody(
    tabItems(
      # First tab content
      tabItem(tabName = 'dashboard',
              h3(HTML("<a href='http://success.wp.shu.edu.tw/1050323%EF%BC%9A%E8%B2%A1%E7%B6%93%E9%96%8B%E6%94%BE%E6%95%B8%E6%93%9A%E5%B9%B3%E5%8F%B0stock-ai%E5%92%8Cr%E5%88%86%E6%9E%90%E4%BB%8B%E9%9D%A2/'>財經開放數據平台Stock AI和R分析介面</a>")),
              fluidRow(
                box(plotOutput('plot1', height = 250)),
                box(
                  title = 'Controls',
                  sliderInput('slider', 'Number of observations:', 1, 100, 50)
                ))),
      # 2nd tab content
      tabItem(tabName = 'line-charts',
              h2('Stocks Prices Line Charts'),
              dateRangeInput("dates", "Date range",
                             start = "2013-01-01", end = as.character(Sys.Date())),
              br(),
              fluidRow(column(1, checkboxInput('888', '888 Holdings PLC', value = FALSE)),
                       column(2, checkboxInput('BET', 'Betfair Group Ltd', value = FALSE)),
                       column(3, checkboxInput('BOTB', 'Best of the Best PLC', value = FALSE))),
              fluidRow(column(1, checkboxInput('BOX', 'Boxhill Technologies PLC', value = FALSE)),
                       column(2, checkboxInput('BPTY', 'Bwin.party Digital Entertainment PLC', value = FALSE)),
                       column(3, checkboxInput('GVC', 'GVC Holdings PLC', value = FALSE))),
              fluidRow(column(1, checkboxInput('LAD', 'Ladbrokes PLC', value = FALSE)),
                       column(2, checkboxInput('NPT', 'Netplay TV PLC', value = FALSE)),
                       column(3, checkboxInput('PCGE', 'Betfair PLC', value = FALSE))),
              fluidRow(column(1, checkboxInput('PTEC', 'Playtech PLC', value = FALSE)),
                       column(2, checkboxInput('RNK', 'Rank Group PLC', value = FALSE)),
                       column(3, checkboxInput('SBT', 'Sportingbet Ltd', value = FALSE))),
              fluidRow(column(1, checkboxInput('SPO', 'Sportech PLC', value = FALSE)),
                       column(2, checkboxInput('STR', 'Stride Gaming PLC', value = FALSE)),
                       column(3, checkboxInput('TNI', 'Trinity Mirror PLC', value = FALSE))),
              fluidRow(column(1, checkboxInput('TTR', '32Red PLC', value = FALSE)),
                       column(2, checkboxInput('WEB', 'Webis Holdings PLC', value = FALSE)),
                       column(3, checkboxInput('WMH', 'William Hill PLC', value = FALSE))),
              mainPanel(
                
              )
              ),
      # 3rd tab content
      tabItem(tabName = 'bar-charts',
              h2('Financial Statements Bar Charts'),
              fluidPage(
                title = 'Examples of DataTables',
                sidebarLayout(
                  sidebarPanel(
                    conditionalPanel(
                      'input.dataset === "diamonds"',
                      checkboxGroupInput('show_vars', 'Columns in diamonds to show:',
                                         names(diamonds), selected = names(diamonds))
                    ),
                    conditionalPanel(
                      'input.dataset === "mtcars"',
                      helpText('Click the column header to sort a column.')
                    ),
                    conditionalPanel(
                      'input.dataset === "iris"',
                      helpText('Display 5 records by default.')
                    ),
                    titlePanel('Downloading Data'),
                    sidebarLayout(
                      sidebarPanel(
                        selectInput("dataset", "Choose a dataset:", 
                                    choices = c("rock", "pressure", "cars")),
                        downloadButton('downloadData', 'Download')
                      ),
                      mainPanel(
                        tableOutput('table')
                      )
                    )),
                  mainPanel(
                    tabsetPanel(
                      id = 'dataset',
                      tabPanel('diamonds', DT::dataTableOutput('mytable1')),
                      tabPanel('mtcars', DT::dataTableOutput('mytable2')),
                      tabPanel('iris', DT::dataTableOutput('mytable3'))
                    )
                  )
                )
              )),
      # 4th tab content
      tabItem(tabName = 'widgets',
              h2('Widgets tab content')),
      # 5th tab content
      tabItem(tabName = 'file-pdf',
              h2('Reference tab content'),
              tabPanel("Help",
                HTML('<iframe src=\"https://englianhu.files.wordpress.com/2016/03/financial-statements-a-step-by-step-guide-to-understanding-and-creating-financial-reports.pdf" 
                         width=\"900\" height=\"600\"></iframe>')
              ),
              imageOutput("imp_pdf",width="500px",height="500px")),
      # 6th tab content
      tabItem(tabName = 'calendar',
              fluidRow(
                column(3, dateRangeInput('dates', label=h3('Date range')))
              )))))

server <- function(input, output) {
  set.seed(122)
  histdata <- rnorm(500)
  
  output$plot1 <- renderPlot({
    data <- histdata[seq_len(input$slider)]
    hist(data)
  })
  
  output$messageMenu <- renderMenu({
    # Code to generate each of the messageItems here, in a list. This assumes
    # that messageData is a data frame with two columns, 'from' and 'message'.
    msgs <- apply(messageData, 1, function(row) {
      messageItem(from = row[['from']], message = row[['message']])
    })
    
    # This is equivalent to calling:
    #   dropdownMenu(type='messages', msgs[[1]], msgs[[2]], ...)
    dropdownMenu(type = 'messages', .list = msgs)
  })
  
  datasetInput <- reactive({
    switch(input$dataset,
           "rock" = rock,
           "pressure" = pressure,
           "cars" = cars)})
  
  output$table <- renderTable({
    datasetInput()
  })
  
  output$downloadData <- downloadHandler(
    filename = function() { 
      paste(input$dataset, '.csv', sep='') 
    },
    content = function(file) {
      write.csv(datasetInput(), file)
    })
  
  # Provide explicit colors for regions, so they don't get recoded when the
  # different series happen to be ordered differently from year to year.
  # http://andrewgelman.com/2014/09/11/mysterious-shiny-things/
  defaultColors <- c("#3366cc", "#dc3912", "#ff9900", "#109618", "#990099", "#0099c6", "#dd4477")
  series <- structure(
    lapply(defaultColors, function(color) { list(color=color) }),
    names = levels(iris$Species)
  )
  
  yearData <- reactive({
    # Filter to the desired year, and put the columns
    # in the order that Google's Bubble Chart expects
    # them (name, x, y, color, size). Also sort by region
    # so that Google Charts orders and colors the regions
    # consistently.
    df <- data %.% arrange(Sepal.Length)
  })
  
  output$chart <- reactive({
    # Return the data and options
    list(
      data = googleDataTable(yearData()),
      options = list(
        title = sprintf(
          "Health expenditure vs. life expectancy, %s",
          input$Species),
        series = series
      )
    )
  })

}

shinyApp(ui, server)


