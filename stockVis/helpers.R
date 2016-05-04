if (!exists(".inflation")) {
  .inflation <- getSymbols('CPIAUCNS', src = 'FRED', 
     auto.assign = FALSE)
}

# adjusts Google/Yahoo finance data with the monthly consumer price index 
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

# 
getCounter <- function(ticked){
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
}

