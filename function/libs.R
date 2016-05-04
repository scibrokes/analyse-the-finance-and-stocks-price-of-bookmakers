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
pkgs <- c('devtools', 'stringr', 'reshape', 'reshape2', 'sparkline', 'DT', 'plyr', 'dplyr', 'purrr', 
          'magrittr', 'foreach', 'knitr', 'rmarkdown', 'scales', 'lubridate', 'tidyr', 'memoise', 
          'doParallel', 'rCharts', 'quantmod', 'TTR', 'farecast', 'shiny')
suppressAll(lib(pkgs)); rm(pkgs)

## Creating a parallel computing Cluster and support functions.
## Preparing the parallel cluster using the cores
doParallel::registerDoParallel(cores = 16)

## Set the rCharts options first to change the behaviour of plot, so that only the chart 
##  component of the HTML file is written into the output file.


## Load the functions
funs <- c('adjust.R', 'global.R')
l_ply(funs, function(x) source(paste0('function/', x))); rm(funs)

