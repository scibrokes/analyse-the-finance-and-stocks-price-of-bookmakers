## Stocks analysis
## the stocks price analysis page for shiny apps
## 
## Please refer to shiny apps inside 'report/app.R'.
## 
## Setup Options, Loading Required Libraries and Preparing Environment
## Loading the packages and setting adjustment
suppressMessages(source('libs.R'))

## Due to bwin-party and Betfair has M&A recently and close market for stocks trading, 
##  here I unable to get the data.
what_metrics <- yahooQF(c('Price/Sales', 'P/E Ratio', 'Price/EPS Estimate Next Year',
                          'PEG Ratio', 'Dividend Yield', 'Market Capitalization'))

tickers <- c('BET', 'BPTY', 'WMH', 'SPO', '888', 'PTEC', 'PCGE', 'TTR', 'GVC', 'WEB', 
             'BOX', 'RNK', 'LAD', 'STR')

comp <- llply(tickers, function(x){
  y <- data.frame(Com=x, as.data.frame(getSymbols(x, auto.assign=FALSE, src='google'))) %>% 
    mutate(Date=ymd(rownames(.)), Weekday=factor(weekdays(Date))); 
  names(y) <- c('Com', 'Open', 'High', 'Low', 'Close', 'Volume', 'Date', 'Weekday'); 
  y}) %>% rbind_all %>% tbl_df %>% 
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





