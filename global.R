if(!require('memoise')) suppressMessages(devtools::install_github('hadley/memoise'))
require('memoise', quietly=TRUE)
suppressMessages(source('function/getData.R'))

## The list of valid counters
## Using <<- but not <- for globally use but not locally use, the objects will be load inside environment and 
##   be used anytime across all functions.
tickers <<- list('888', 'BET', 'BOTB', 'BOX', 'BPTY', 'BYD', 'CACQ', 'CWLDY', 'EGT', 'GDEN', 'GEBHY', 'GIGNY', 
              'GLXZ', 'GPIC', 'GVC', 'GXYEY', 'IGT', 'ISLE', 'LAD', 'LVS', 'MCHVY', 'MCRI', 'MGM', 'MPEL', 
              'NGCRF', 'NPT', 'PCGE', 'PDSSF', 'PENN', 'PNK', 'PTEC', 'RNK', 'SCHYF', 'SGMS', 'SJMHF', 'SKYZF', 
              'SPO', 'STR', 'TACBY', 'TNI', 'TTR', 'WEB', 'WMH', 'WYNMF', 'WYNN')

## Using "memoise" to automatically cache the results
filterData <- memoise(function(ticker=ticker, from=from, to=to) {
  ## Careful not to let just any name slip in here; a malicious user could manipulate this value.
  if(!(ticker %in% tickers)){
    stop("Unknown ticker! Please select tickers among : '888', 'BET', 'BOTB', 'BOX', 'BPTY', 'BYD', 'CACQ', 
          'CWLDY', 'EGT', 'GDEN', 'GEBHY', 'GIGNY', 'GLXZ', 'GPIC', 'GVC', 'GXYEY', 'IGT', 'NGCRF', 'NPT', 
          'ISLE', 'LAD', 'LVS', 'MCHVY', 'MCRI', 'MGM', 'MPEL', 'PCGE', 'PDSSF', 'PENN', 'PNK', 'PTEC', 
          'RNK', 'SCHYF', 'SGMS', 'SJMHF', 'SKYZF', 'SPO', 'STR', 'TACBY', 'TNI', 'TTR', 'WEB', 'WMH', 
          'WYNMF', 'WYNN'")
  }else{
    getData(ticker=ticker, from=from, to=to)
  }
})

