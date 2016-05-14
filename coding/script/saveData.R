## setwd()
pth1 <- 1:2; sheet = c('2006-07 Source')#, '2007-08 Source')
#pth <- paste0(<mypath>, pth1, '/League Listing')
rm(pth1)

source('./coding/function/readData.R')
SB0607 <- readData(pth = pth, sht = sheet)
saveRDS(SB0607, './data/report.rds')
rm(pth, sheet, getData, readData)

