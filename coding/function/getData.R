## getData
getData <- function(pth, sht){
  suppressMessages(library('plyr'))
  suppressMessages(library('dplyr'))
  suppressMessages(library('magrittr'))
  suppressMessages(library('readxl'))
  
  mdata = llply(seq(dir(pth)), function(i) {
    y = read_excel(paste0(pth, '/', dir(pth)[[i]]), sheet = sht) %>% 
      data.frame %>% tbl_df
    nm = as.character(y[1, 1])
    y %<>% .[-c(1:4), ]#-c(2:21)]
    return(list(data = y, names = nm))
  })
  
  vdata = llply(seq(length(mdata)), function(i){
    names(mdata[[i]]) = mdata[[i]][[2]]; mdata[[i]][[1]]
  })
  
  names(vdata) = sapply(mdata, `[[`, 2)
  rm(pth, mdata)
  return(vdata)
}

