## Setup Options, Loading Required Libraries and Preparing Environment

## Setting to omit all warnings
options(warn=-1)

## Loading the package 'BBmisc'
if(suppressMessages(!require('BBmisc'))){
  install.packages('BBmisc')}
suppressMessages(library('BBmisc'))

## http://www.r-bloggers.com/new-package-dplyrr-utilities-for-comfortable-use-of-dplyr-with-databases
## direct connect to database (if any)
#'@ if(!'dplyrr' %in% installed.packages()){
#'@   devtools::install_github("hoxo-m/dplyrr")}
#'@ install.packages('nycflights13')
#'@ library(c('dplyrr', 'nycflights13'))

## Loading multiple packages at once
pkgs <- c('devtools', 'zoo', 'chron', 'stringr', 'reshape', 'reshape2', 'data.table', 'sparkline', 'stringdist',
          'DT', 'plyr', 'dplyr', 'purrr', 'magrittr', 'foreach', 'manipulate', 'ggplot2', 'PerformanceAnalytics',
          'highlightHTML', 'knitr', 'rmarkdown', 'scales', 'lubridate', 'tidyr', 'whisker', 'gtable', 'grid',
          'doParallel', 'gridExtra', 'googleVis', 'quantmod', 'TTR')
suppressAll(lib(pkgs)); rm(pkgs)

## Creating a parallel computing Cluster and support functions.
## Preparing the parallel cluster using the cores
doParallel::registerDoParallel(cores = 16)
#'@ BiocParallel::register(MulticoreParam(workers=8))

## Set the googleVis options first to change the behaviour of plot.gvis, so that only the chart 
##  component of the HTML file is written into the output file.
op <- options(gvis.plot.tag='chart')


