if(!require('memoise')) suppressMessages(devtools::install_github('hadley/memoise'))
suppressMessages(require('memoise', quietly = TRUE))

## The list of valid counters
## Using <<- but not <- for globally use but not locally use, the objects will be load inside environment and 
##   be used anytime across all functions.
tickers <<- c('BET', 'BOTB', 'BPTY', 'WMH', 'SPO', '888', 'NPT', 'PTEC', 'PCGE', 'TTR', 'GVC', 'WEB', 
             'BOX', 'RNK', 'LAD', 'STR', 'LVS', 'TNI', 'WYNN', 'MGM', 'MPEL', 'MCRI', 'GXYEY', 'BYD', 
             'GIGNY', 'CACQ', 'PENN', 'GEBHY', 'WYNMF', 'CWLDY', 'SCHYF', 'SKYZF', 'MCHVY', 'EGT', 
             'GPIC', 'SGMS', 'NGCRF', 'TACBY', 'IGT', 'GDEN', 'ISLE', 'GLXZ', 'SJMHF', 'PDSSF', 'PNK') %>% sort

## input <<- list(counter=tickers[sample(c(0,1), length(tickers), replace=TRUE)==1], 
##                dates=c('2015-10-10', as.character(Sys.Date())))

## Using "memoise" to automatically cache the results
## http://www.r-bloggers.com/fibonacci-sequence-in-r-with-memoization/
filterData <- memoise(function(ticker) {
  ## Careful not to let just any name slip in here; a malicious user could manipulate this value.
  
  if(any(!(ticker %in% tickers))){
    stop("Unknown ticker! Please select tickers among : '888', 'BET', 'BOTB', 'BOX', 'BPTY', 'BYD', 'CACQ', 
          'CWLDY', 'EGT', 'GDEN', 'GEBHY', 'GIGNY', 'GLXZ', 'GPIC', 'GVC', 'GXYEY', 'IGT', 'NGCRF', 'NPT', 
          'ISLE', 'LAD', 'LVS', 'MCHVY', 'MCRI', 'MGM', 'MPEL', 'PCGE', 'PDSSF', 'PENN', 'PNK', 'PTEC', 
          'RNK', 'SCHYF', 'SGMS', 'SJMHF', 'SKYZF', 'SPO', 'STR', 'TACBY', 'TNI', 'TTR', 'WEB', 'WMH', 
          'WYNMF', 'WYNN'")
  }else{
    getData(ticker=ticker)
  }
})

