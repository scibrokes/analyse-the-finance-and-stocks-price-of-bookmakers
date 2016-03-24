suppressMessages(library('BBmisc'))
pkgs <- c('shiny', 'shinythemes', 'shinyjs', 'shinyBS', 'shinydashboard', 'shinyAce', 'quantmod', 
          'TTR', 'plyr', 'dplyr', 'stringr', 'purrr', 'googleCharts', 'lubridate', 'googleVis', 'zoo')
rm(pkgs)

tickers <- c('BET', 'BPTY', 'WMH', 'SPO', '888', 'PTEC', 'PCGE', 'TTR', 'GVC', 'WEB', 
             'BOX', 'RNK', 'LAD', 'STR') %>% sort

## =================================================================================
messageData <- data_frame(from=c('Sayaka', '鬼塚先生', '小松拓也'), 
                          message=c('Your song: Ever Since', 'Gaming introduction',
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
                                         taskItem(value = 80, color = 'red', 'Overall project')),
                            tags$li(class = 'dropdown',
                                    tags$a(href='http://www.scibrokes.com', target='_blank', 
                                           tags$img(height = '20px', alt='Scibrokes Logo', align='right', 
                                                    src='https://avatars0.githubusercontent.com/u/13562894?v=3&s=200')
                                    )
                            ))

## ---------------------------------------------------------------------------------
ui <- dashboardPage(
  dbHeader,
  ## Sidebar content
  dashboardSidebar(
    sidebarMenu(
      menuItem('Dashboard', tabName = 'dashboard', icon = icon('dashboard')),
      menuItem('Widgets', tabName = 'widgets', icon = icon('th')),
      menuItem('The Economics of Soccer', tabName = 'soccer', icon = icon('bar-chart'),
               menuSubItem('Summary Charts', tabName='report1', icon = icon('bar-chart')),
               menuSubItem('Econometrics in Soccer Betting', tabName='econometrics', icon = icon('bar-chart'))),
      menuItem('Company Profile', tabName = 'profile', icon = icon('bar-chart'),
               menuSubItem('Annual Summary Charts', tabName='report2', icon = icon('bar-chart')),
               menuSubItem('Stock Market', tabName='stocks', icon = icon('line-chart')),
               menuSubItem('Financial Report', tabName='finance', icon = icon('bar-chart'))),
      menuItem('References', tabName = 'reference', icon = icon('file-pdf-o')),
      menuItem('Source code', icon = icon('file-code-o'), 
               href = 'https://github.com/scibrokes/analyse-the-finance-and-stocks-price-of-bookmakers')
    )),
  ## Body content
  dashboardBody(
    tabItems(
      # First tab content
      tabItem(tabName = 'dashboard',
              h3(HTML("<a href='http://success.wp.shu.edu.tw/1050323%EF%BC%9A%E8%B2%A1%E7%B6%93%E9%96%8B
                      %E6%94%BE%E6%95%B8%E6%93%9A%E5%B9%B3%E5%8F%B0stock-ai%E5%92%8Cr%E5%88%86%E6%9E%90%
                      E4%BB%8B%E9%9D%A2/'>財經開放數據平台Stock AI和R分析介面</a>")),
              fluidRow(
                box(plotOutput('plot1', height = 250)),
                box(
                  title = 'Controls',
                  sliderInput('slider', 'Number of observations:', 1, 100, 50)
                ))),
      # 2nd tab content
      tabItem(tabName = 'widgets',
              h2('Widget content')),
      # 3rd tab content
      tabItem(tabName = 'soccer',
              h2('The Economics of Soccer')),
      # 4th tab content
      tabItem(tabName = 'report1',
              h2('Summary Charts')),
      # 5th tab content
      tabItem(tabName = 'econometrics',
              h2('Econometrics in Soccer Betting')),
      # 6th tab content
      tabItem(tabName = 'profile',
              h2('Company Profile')),
      # 7th tab content
      tabItem(tabName = 'report2',
              h2('Annual Summary Chart')),
      # 8th tab content
      tabItem(tabName = 'stocks',
              h2('Stock Market'),
              fluidPage(
                titlePanel('Stock Price'),
                sidebarLayout(
                  sidebarPanel(
                    helpText('Select the stocks to examine. Information will be collected from google finance.'),
                    checkboxGroupInput('counter', 'Select counter(s)', choices = tickers), 
                    #'@ checkboxInput('symb01', '888 : 01) 888 Holdings PLC', value = FALSE),
                    #'@ checkboxInput('symb02', 'BET : 02) Betfair Group Ltd', value = FALSE),
                    #'@ checkboxInput('symb03', 'BOTB : 03) Best of the Best PLC', value = FALSE),
                    #'@ checkboxInput('symb04', 'BOX : 04) Boxhill Technologies PLC', value = FALSE),
                    #'@ checkboxInput('symb05', 'BPTY : 05) Bwin.party Digital Entertainment PLC', value = FALSE),
                    #'@ checkboxInput('symb06', 'GVC : 06) GVC Holdings PLC', value = FALSE),
                    #'@ checkboxInput('symb07', 'LAD : 07) Ladbrokes PLC', value = FALSE),
                    #'@ checkboxInput('symb08', 'NPT : 08) Netplay TV PLC', value = FALSE),
                    #'@ checkboxInput('symb09', 'PCGE : 09) Betfair PLC', value = FALSE),
                    #'@ checkboxInput('symb10', 'PTEC : 10) Playtech PLC', value = FALSE),
                    #'@ checkboxInput('symb11', 'RNK : 11) Rank Group PLC', value = FALSE),
                    #'@ checkboxInput('symb12', 'SBT : 12) Sportingbet Ltd', value = FALSE),
                    #'@ checkboxInput('symb13', 'SPO : 13) Sportech PLC', value = FALSE),
                    #'@ checkboxInput('symb14', 'STR : 14) Stride Gaming PLC', value = FALSE),
                    #'@ checkboxInput('symb15', 'TNI : 15) Trinity Mirror PLC', value = FALSE),
                    #'@ checkboxInput('symb16', 'TTR : 16) 32Red PLC', value = FALSE),
                    #'@ checkboxInput('symb17', 'WEB : 17) Webis Holdings PLC', value = FALSE),
                    #'@ checkboxInput('symb18', 'WMH : 18) William Hill PLC', value = FALSE),
                    actionLink('selectall', 'Select All'),
                    br(),
                    br(),
                    dateRangeInput('dates', 'Date range',
                                   start = '2013-01-01', 
                                   end = as.character(Sys.Date())),
                    br(),
                    checkboxInput('log', 'Plot y axis on log scale', value = FALSE),
                    checkboxInput('adjust', 'Adjust prices for inflation', value = FALSE)
                    ),
                  mainPanel(plotOutput("plot"))))),
      # 9th tab content
      tabItem(tabName = 'finance',
              h2('Financial Report')),
      # 10th tab content
      tabItem(tabName = 'reference',
              h2('Reference:'),
              tabPanel("Help",
                       HTML('<iframe src=\"https://englianhu.files.wordpress.com/2016/03/financial-statements-a-step-by-step-guide-to-understanding-and-creating-financial-reports.pdf" width=\"900\" height=\"600\"></iframe>')),
              imageOutput('imp_pdf', width='500px', height='500px')))))


