getData <- function(ticker, from='2007-01-01', to=as.character(Sys.Date())){
  
  if(!require('BBmisc')) suppressMessages(install.packages('BBmisc'))
  suppressMessages(BBmisc::lib(c('plyr', 'dplyr', 'tidyr', 'lubridate', 'TTR', 'quantmod', 'doParallel')))
  
  ## Load below packages globally (system path) in case of BBmisc::lib() only install locally (user path) 
  ##   and doesn't work during runApp() shinyapps. 
  #'@ suppressMessages(library('plyr'))
  #'@ suppressMessages(library('dplyr'))
  #'@ suppressMessages(library('tidyr'))
  #'@ suppressMessages(library('quantmod'))
  #'@ suppressMessages(library('TTR'))
  #'@ suppressMessages(library('lubridate'))
  #'@ suppressMessages(library('doParallel'))
  
  ## Creating a parallel computing Cluster and support functions.
  ## Preparing the parallel cluster using the cores
  doParallel::registerDoParallel(cores = 16)
  
  data <- llply(ticker, function(x){
    y = as.data.frame(
      tryCatch(suppressAll(getSymbols(x, src = 'yahoo', auto.assign = FALSE, from=from, to=to)), 
               error = function(e) data.frame(Open=NA, High=NA, Low=NA, Close=NA, Volume=NA, Adjusted=NA)))
    if(nrow(na.omit(y))==0){
      y = as.data.frame(
        tryCatch(suppressAll(getSymbols(x, src = 'google', auto.assign = FALSE, 
                                        from=from, to=to)), error = function(e) 
                                          data.frame(Open=NA, High=NA, Low=NA, Close=NA, Volume=NA, Adjusted=NA)))
      y$Adjusted <- as.numeric(NA)
      names(y)[length(y)] <- paste0(x, '.Adjusted')
      names(y) <- gsub('X888', '888', names(y))
    }
    
    z = y %>% mutate(Date=as.character(rownames(.)), Weekday=factor(weekdays(ymd(Date)))) %>% 
      select(Date, Weekday)
    dt = data.frame(z, y) %>% tbl_df
  }, .parallel=TRUE) %>% Reduce(function(x, y) merge(x, y, all=T, by=c('Date', 'Weekday')), ., accumulate=F) %>% 
    tbl_df %>% filter(nchar(Date)==10 & substr(Date,5,5)=='-' & substr(Date,8,8)=='-')
  names(data) <- gsub(substr(grep('X888', names(data), value=TRUE),1,4),'888', names(data))
  
  cate <- c('Open', 'High', 'Low', 'Close', 'Volume', 'Adjusted')
  dataPlotLine <- llply(cate, function(x) data[c('Date', 'Weekday', grep(x, names(data), value=TRUE))])
  names(dataPlotLine) <- cate
  
  dataComb <- llply(seq(cate), function(i) {
    dataPlotLine[[i]] %>% gather(Com, Price, -c(Date:Weekday))
  }) %>% rbind_all %>% tbl_df # %>% arrange(Date, Com)
  dataComb2 <- dataComb %>% separate(Com, c('Com', 'Status')) %>% na.omit %>% spread(Status, Price) %>% 
    arrange(Date, Com) %>% select(Date, Weekday, Com, Open, High, Low, Close, Volume, Adjusted)
  #'@ dataComb %>% spread(Com, Price) #similar with the above dataframe object 'data'
  
  return(list(dataTable=dataComb2, ggPlot=dataComb, listPlot=dataPlotLine, dataPlot=data))
}

