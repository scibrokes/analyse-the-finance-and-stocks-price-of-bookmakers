## Stocks analysis
## the stocks price analysis page for shiny apps
## 
## Please refer to shiny apps inside 'report/app.R'.
## 
## Setup Options, Loading Required Libraries and Preparing Environment
## Loading the packages and setting adjustment
suppressMessages(library('BBmisc'))
pkgs <- c('shiny', 'shinythemes', 'shinyjs', 'shinyBS', 'shinydashboard', 'shinyAce', 'quantmod', 
          'TTR', 'plyr', 'dplyr', 'stringr', 'purrr', 'googleCharts')
suppressMessages(lib(pkgs))
suppressMessages(source('libs.R'))
rm(pkgs)

## Due to bwin-party and Betfair has M&A recently and close market for stocks trading, 
##  here I unable to get the data.
what_metrics <- yahooQF(c('Price/Sales', 'P/E Ratio', 'Price/EPS Estimate Next Year',
                          'PEG Ratio', 'Dividend Yield', 'Market Capitalization'))

tickers <- c('BET', 'BOTB', 'BPTY', 'WMH', 'SPO', '888', 'NPT', 'PTEC', 'PCGE', 'TTR', 'GVC', 'WEB', 
             'BOX', 'RNK', 'LAD', 'STR', 'LVS', 'TNI', 'WYNN', 'MGM', 'MPEL', 'MCRI', 'GXYEY', 'BYD', 
             'GIGNY', 'CACQ', 'PENN', 'GEBHY', 'WYNMF', 'CWLDY', 'SCHYF', 'SKYZF', 'MCHVY', 'EGT', 
             'GPIC', 'SGMS', 'NGCRF', 'TACBY', 'IGT', 'GDEN', 'ISLE', 'GLXZ', 'SJMHF', 'PDSSF', 'PNK') %>% sort

data <- llply(tickers, function(x){
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

# Not all the metrics are returned by Yahoo.
metrics <- getQuote(paste(tickers, sep='', collapse=';'), what=what_metrics)

#Add tickers as the first column and remove the first column which had date stamps
metrics <- data.frame(Symbol=tickers, metrics[, 2:length(metrics)]) 

#Change colnames
colnames(metrics) <- c('Symbol', 'Revenue Multiple', 'Earnings Multiple', 
                       'Earnings Multiple (Forward)', 'Price-to-Earnings-Growth', 
                       'Div Yield', 'Market Cap')

rm(what_metrics)
#'@ getFinancials(tickers)


if (!exists(".inflation")) {
  .inflation <- getSymbols('CPIAUCNS', src = 'FRED', 
                           auto.assign = FALSE)
}  

# adjusts yahoo finance data with the monthly consumer price index 
# values provided by the Federal Reserve of St. Louis
# historical prices are returned in present values 
adjust <- function(data) {
  
  latestcpi <- last(.inflation)[[1]]
  inf.latest <- time(last(.inflation))
  months <- split(data)               
  
  adjust_month <- function(month) {               
    date <- substr(min(time(month[1]), inf.latest), 1, 7)
    coredata(month) * latestcpi / .inflation[date][[1]]
  }
  
  adjs <- lapply(months, adjust_month)
  adj <- do.call("rbind", adjs)
  axts <- xts(adj, order.by = time(data))
  axts[ , 5] <- Vo(data)
  axts
}


