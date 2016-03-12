#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

suppressMessages(library('BBmisc'))
pkgs <- c('shiny', 'shinythemes', 'shinyjs', 'shinydashboard', 'shinyAce', 'quantmod', 
          'TTR', 'plyr', 'dplyr', 'stringr', 'purrr')
suppressMessages(lib(pkgs))
rm(pkgs)

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
                        <img align='right' alt='Logo' src='./oda-army.jpg'/></a>")

ui <- dashboardPage(
  dbHeader,
  
  ## Sidebar content
  dashboardSidebar(
    sidebarMenu(
      menuItem('Dashboard', tabName = 'dashboard', icon = icon('dashboard')),
      menuItem('Widgets', tabName = 'widgets', icon = icon('th')),
      menuItem('Stocks Prices', tabName = 'line-charts', icon = icon('line-chart')),
      menuItem('Financial Statements', tabName = 'bar-charts', icon = icon('bar-chart')),
      menuItem('Calendar', tabName = 'calendar', icon = icon('calendar')),
      menuItem('References', tabName = 'file-pdf', icon = icon('file-pdf')),
      menuItem('Source code', icon = icon('file-code-o'), 
               href = 'https://github.com/scibrokes/analyse-the-finance-and-stocks-price-of-bookmakers')
    )),
  ## Body content
  dashboardBody(
    tabItems(
      # First tab content
      tabItem(tabName = 'dashboard',
              fluidRow(
                box(plotOutput('plot1', height = 250)),
                box(
                  title = 'Controls',
                  sliderInput('slider', 'Number of observations:', 1, 100, 50)
                ))),
      # 2nd tab content
      tabItem(tabName = 'line-charts',
              h2('Stocks Prices Line Charts')),
      # 3rd tab content
      tabItem(tabName = 'bar-charts',
              h2('Financial Statements Bar Charts')),
      # 4th tab content
      tabItem(tabName = 'widgets',
              h2('Widgets tab content')),
      # 5th tab content
      tabItem(tabName = 'file-pdf',
              h2('Reference tab content')),
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
}

shinyApp(ui, server)
