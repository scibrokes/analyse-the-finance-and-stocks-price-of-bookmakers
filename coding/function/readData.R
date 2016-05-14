readData <- function(pth, sht) {
  ## Only annual report of SB1 or SB2 without looping
  ##  as I didn't make a nested looping in read_excel() inside getData(), 
  ##  and below df_to_blocks()
  if(!is.vector(pth) && (length(pth) < 1))
    stop('`pth` is the path of workbook(s) which must be a vector and length > 0.')
  if(!is.vector(sht) && (length(sht) != 1))
    stop('`sht` is the worksheet name and it must be a vector and length == 1.')
  
  ## get Data
  source('./coding/function/getData.R')
  
  ## Load packages
  suppressMessages(library('plyr'))
  suppressMessages(library('dplyr'))
  suppressMessages(library('magrittr'))
  suppressMessages(library('turner'))
  suppressMessages(library('zoo'))
  
  vdata <- llply(seq(pth), function(i) {
    getData(pth = pth[i], sht = sht)
    })
  
  names(vdata) <- c('SB 1', 'SB 2')
  rm(pth, sht)
  
  com1 <- NULL
  for(i in seq(length(vdata))) {
    ## df_to_blocks only allow one data frame to be splited to multiple data frames 
    ##   within a list
    com1[[i]] <- llply(vdata[[i]], function(x) {
      df_to_blocks(data.frame(x), blocks = c(1, rep(4, (ncol(x) - 1)/4)), byrow = FALSE)
      })
  }
  rm(vdata, i)
  
  dat <- suppressWarnings(llply(seq(length(com1)), function(i) {
    llply(seq(length(com1[[i]])), function(j) {
      dat = llply(seq(length(com1[[i]][[j]])), function(k) {
        y = data.frame(Month = com1[[i]][[j]][[k]][[1]][1], com1[[i]][[j]][[1]], 
                       com1[[i]][[j]][[k]]) %>% .[c(2:16, 47:60), ]
        if(ncol(y) == 6)
          names(y) = c('Month', 'Bet_Type', 'Turnover', 'Ticket', 'Comp_Turn', 'PL')
        y = y[-1, ]
        data.frame(Pay_Type = c(rep('Credit', nrow(y)/2), rep('Cash', nrow(y)/2)), y)
        })
      dat[[1]] = NULL
      dat %>% rbind_all
      }) %>% rbind_all %>% 
      mutate(Month = strsplit(Month, ' ')[[1]] %>% paste(substr(.[1], 1, 3), .[2]) %>% 
               .[2] %>% substring(6) %>% as.yearmon, Turnover = as.numeric(Turnover), 
             Ticket = as.numeric(Ticket), Comp_Turn = as.numeric(Comp_Turn), 
             PL = as.numeric(PL))
    }))
  rm(com1)
  names(dat) <- c('SB 1', 'SB 2')
  
  return(dat)
  }

