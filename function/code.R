###############################################################################
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
###############################################################################
# Systematic Investor Toolbox (SIT)
# Copyright (C) 2012  Michael Kapler
#
# For more information please visit my blog at 
# www.SystematicInvestor.wordpress.com
# or drop me a line at TheSystematicInvestor at gmail
###############################################################################
# Systematic Investor Toolbox (SIT)
#
# Systematic Investor Toolbox is a collection of tools that I use
# in my investment research. I will demonstrate and document 
# various uses of toolbox in the Systematic Investor blog at
#	www.SystematicInvestor.wordpress.com
#
#
###############################################################################
# Example Usage:
###############################################################################
#
#
###############################################################################
# Load Systematic Investor Toolbox (SIT)
# http://systematicinvestor.wordpress.com/systematic-investor-toolbox/
###############################################################################
#con = gzcon(url('http://www.systematicportfolio.com/sit.gz', 'rb'))
#source(con)
#close(con)
#
#
###############################################################################
# Load Systematic Investor Toolbox (SIT): Windows only
############################################################################### 
# Load Systematic Investor Toolbox (SIT)
#setInternet2(TRUE)
#con = gzcon(url('https://github.com/systematicinvestor/SIT/raw/master/sit.gz',
#                 'rb'))
#source(con)
#close(con)
#
#
###############################################################################
# Load Systematic Investor Toolbox (SIT): from file, if for example you saved 
# sit.gz to c:/temp/sit.gz
############################################################################### 
#con = gzcon(file('c:/temp/sit.gz', 'rb'))
#source(con)
#close(con)
#
#
###############################################################################
# Load Systematic Investor Toolbox (SIT): Requires RCurl package
############################################################################### 
#require(RCurl)
#sit = getURLContent('https://github.com/systematicinvestor/SIT/raw/master/
#                     sit.gz', binary=TRUE, followlocation = TRUE, 
#                     ssl.verifypeer = FALSE)
#con = gzcon(rawConnection(sit, 'rb'))
#source(con)
#close(con)
#
#
###############################################################################
# Example Usage:
############################################################################### 
# Run plota test
#plota.test()
#
#
#
#
#
#
#More to come,
#
#Michael Kapler
#TheSystematicInvestor at gmail
#

bl.compute.risk.aversion <- function(bench, risk.free = 0) {
  lambda = mean(coredata(bench) - coredata(risk.free)) / var(coredata(bench))
  return( as.double(lambda))
  }

bl.compute.eqret <- function(risk.aversion, cov, cap.weight, risk.free = 0) {
  return(risk.aversion * cov %*% as.vector(cap.weight) + as.double(risk.free))
  }

bl.compute.posterior <- function(mu, cov, pmat = NULL, qmat = NULL, 
                                 tau = 0.025, confidences = NULL){
  out = list()
  if(!is.null(pmat)) {
    if(is.null(confidences)) {
      omega = diag(c(1, diag(tau * pmat %*% cov %*% t(pmat))))[-1, -1]
      } else {
        omega = diag(c(1,confidences))[-1,-1]
        }
    temp = solve(solve(tau * cov) + t(pmat) %*% solve(omega) %*% pmat)
    out$cov = cov + temp
    out$expected.return = temp %*% (solve(tau * cov) %*% mu + t(pmat) %*% 
                                      solve(omega) %*% qmat)
    } else {
      temp = tau * cov
      out$cov = cov + temp
      out$expected.return = temp %*% (solve(tau * cov) %*% mu)
      }
  return(out)
  }

bl.compute.optimal <- function(risk.aversion, mu, cov){
  return( (1/risk.aversion) * solve(cov) %*% mu )
  }

aa.black.litterman.examples <- function(){
  data = '1,0.4880,0.4780,0.5150,0.4390,0.5120,0.4910
0.4880,1,0.6640,0.6550,0.3100,0.6080,0.7790
0.4780,0.6640,1,0.8610,0.3550,0.7830,0.6680
0.5150,0.6550,0.8610,1,0.3540,0.7770,0.6530
0.4390,0.3100,0.3550,0.3540,1,0.4050,0.3060
0.5120,0.6080,0.7830,0.7770,0.4050,1,0.6520
0.4910,0.7790,0.6680,0.6530,0.3060,0.6520,1'
  
  Corrmat = matrix(as.double(spl(gsub('\n', ',', data), ',')),
                   nrow = len(spl(data, '\n')), byrow = TRUE)
  
  RiskAversion = 2.5
  stdevs = c(16.0, 20.3, 24.8, 27.1, 21.0,  20.0, 18.7)/100
  MktWeight = c(1.6, 2.2, 5.2, 5.5, 11.6, 12.4, 61.5)/100
  tau = 0.05
  Covmat = Corrmat * (stdevs %*% t(stdevs))
  EqRiskPrem = RiskAversion * Covmat %*% MktWeight
  EqRiskPrem = bl.compute.eqret(RiskAversion, Covmat, MktWeight)
  AssetNames = c('Australia', 'Canada', 'France', 'Germany', 'Japan', 'UK', 
                 'USA')
  Table2 = cbind(AssetNames, 
                 round(cbind(stdevs, MktWeight, EqRiskPrem) * 100, 1))
  colnames(Table2) = c('Assets','Std Dev','Weq','PI')
  Table2
  
  P = matrix(c(0, 0, -29.5, 100, 0, -70.5, 0)/100, nrow=1)
  Q = 5/100
  Omega = diag(c(1,diag(tau * P %*% Covmat %*% t(P))))[-1,-1]
  PostCov = solve(solve(tau*Covmat) + (t(P) %*% solve(Omega) %*% P))
  SigmaP = Covmat + PostCov
  ExpRet = PostCov %*% (solve(tau * Covmat) %*% EqRiskPrem + t(P) %*% 
                          solve(Omega) %*% Q)
  post = bl.compute.posterior(EqRiskPrem, Covmat, P, Q, tau)
  ExpRet = post$expected.return
  
  SigmaP = post$cov
  OptimalWeights = (1/RiskAversion) * solve(SigmaP) %*% ExpRet
  OptimalWeights = bl.compute.optimal(RiskAversion, ExpRet, SigmaP)
  Tab4Col4 = OptimalWeights - (MktWeight)/(1+tau)
  Table4 = cbind(AssetNames, 
                 round(cbind(t(P), ExpRet, OptimalWeights, 
                             round(Tab4Col4 * 1000)/1000) * 100, 1))
  colnames(Table4) = c('Assets', 'P', 'MU', 'W','W - Weq/1+tau')
  Table4
  
  x = c( 0.001005,  0.001328, -0.000579, -0.000675,  0.000121,  0.000128, 
        -0.000445, -0.000437,  0.001328,  0.007277, -0.001307, -0.000610, 
        -0.002237, -0.000989,  0.001442, -0.001535, -0.000579, -0.001307, 
         0.059852,  0.027588,  0.063497,  0.023036,  0.032967,  0.048039, 
        -0.000675, -0.000610,  0.027588,  0.029609,  0.026572,  0.021465, 
         0.020697,  0.029854,  0.000121, -0.002237,  0.063497,  0.026572, 
         0.102488,  0.042744,  0.039943,  0.065994,  0.000128, -0.000989, 
         0.023036,  0.021465,  0.042744,  0.032056,  0.019881,  0.032235, 
        -0.000445,  0.001442,  0.032967,  0.020697,  0.039943,  0.019881, 
         0.028355,  0.035064, -0.000437, -0.001535,  0.048039,  0.029854, 
         0.065994,  0.032235,  0.035064,  0.079958)
  
  varCov <- matrix(x, ncol = 8, nrow = 8)
  mu <- c(0.08, 0.67,6.41, 4.08, 7.43, 3.70, 4.80, 6.60) / 100
  pick <- matrix(0, ncol = 8, nrow = 3, dimnames = list(NULL, letters[1:8]))
  pick[1,7] <- 1
  pick[2,1] <- -1; pick[2,2] <- 1
  pick[3, 3:6] <- c(0.9, -0.9, .1, -.1)
  post = bl.compute.posterior(mu, varCov, pick, c(0.0525, 0.0025, 0.02), 
                              tau = 0.025)
  
  library(BLCOP)
  confidences <- 1 / c(0.000709, 0.000141, 0.000866)
  myViews <- BLViews(pick, c(0.0525, 0.0025, 0.02), confidences, letters[1:8])
  myPosterior <- posteriorEst(myViews, tau = 0.025, mu, varCov )
  myPosterior
  myPosterior@posteriorMean - post$expected.return
  myPosterior@posteriorCovar - post$cov
  }

Rglpk.read.model <- 
  function(file, type = c("MPS_fixed", "MPS_free", "CPLEX_LP", "MathProg"), 
           ignore_first_row = FALSE, verbose = FALSE) {
    
    if(!file.exists(file)) stop(paste("There is no file called", file, "!"))
    type_db <- c("MPS_fixed" = 1L, "MPS_free" = 2L, "CPLEX_LP" = 3L, 
                 "MathProg" = 4L)
  obj <- list(file = tools::file_path_as_absolute(file), 
              type = type_db[match.arg(type)])
  meta_data <- Rglpk:::glp_get_meta_data_from_file(obj, verbose)
  milp_data <- Rglpk:::glp_retrieve_MP_from_file(
    meta_data, ignore_first_row, verbose)
  MP_data <- Rglpk:::glp_merge_MP_data(meta_data, milp_data)
  dir_db <- c("free" = 1L, ">=" = 2L, "<=" = 3L, "DB" = 4L, "==" = 5L)
  MP_data$direction_of_constraints <- 
    names(dir_db[MP_data$direction_of_constraints])
  types <- rep("C", length.out = MP_data$n_objective_vars)
  
  if(any(MP_data$objective_var_is_integer)) 
    types[MP_data$objective_var_is_integer] <- "I"
  
  if(any(MP_data$objective_var_is_binary)) 
    types[MP_data$objective_var_is_binary] <- "B"
  
  MP_data$types = types
  index = which(MP_data$direction_of_constraints == 'free')
  
  if(length(index) > 0 ) {
    MP_data$constraint_matrix = as.matrix(MP_data$constraint_matrix)[-index, ]
    MP_data$direction_of_constraints = 
      MP_data$direction_of_constraints[-index]
    MP_data$right_hand_side = MP_data$right_hand_side[-index]
    }
  MP_data
  }

Rglpk.create.constraints <- function(prob) {
  n = prob$n_objective_vars
  lb = rep(NA, n)
  lb[prob$bounds$lower$ind] = prob$bounds$lower$val
  ub = rep(NA, n)
  ub[prob$bounds$upper$ind] = prob$bounds$upper$val
  
  constraints = new.constraints(n, lb = lb, ub = ub)
  constraints$binary.index = which(prob$objective_var_is_binary == 1)
  
  if(len(constraints$binary.index) == 0) constraints$binary.index = 0
  if(is.null(dim(prob$constraint_matrix))) {
    prob$constraint_matrix = matrix(prob$constraint_matrix)
    } else {
      prob$constraint_matrix = t(prob$constraint_matrix)
    }
  
  index = which(prob$direction_of_constraints == '==')
  if(len(index) > 0) constraints = add.constraints(
    prob$constraint_matrix[, index], type = '=', 
    b = prob$right_hand_side[index], constraints)
  
  index = which(prob$direction_of_constraints == '<=')
  
  if(len(index) > 0) constraints = add.constraints(
    prob$constraint_matrix[, index], type = '<=', 
    b = prob$right_hand_side[index], constraints)
  index = which(prob$direction_of_constraints == '>=')
  
  if(len(index) > 0) constraints = add.constraints(
    prob$constraint_matrix[,index], type = '>=', 
    b = prob$right_hand_side[index], constraints)
  
  f.obj = prob$objective_coefficients
  dir = ifelse(prob$maximize, 'max', 'min')
  prob$names = prob$objective_vars_names
  prob$tickers = prob$objective_vars_names
  
  if(len(grep('\\[',prob$objective_vars_names)) > 0) {
    temp = matrix(spl(gsub(']', '', prob$objective_vars_names), '\\['), 
                  nr = 2)
    prob$names = temp[1, ]
    prob$tickers = temp[2, ]
    }
  return(list(constraints = constraints, f.obj = f.obj, dir = dir, 
              prob = prob))
  }

parse.views <- function(symbolnames, views) {
  load.packages('Rglpk')
  views = parse.expr(views)
  
  if (len(views)==0) 
    return(list(A = matrix(0, nr = 0, nc = len(symbolnames)), b = c(), 
                meq = 0))
  
  model.file = tempfile('temp.model')
  on.exit(unlink(model.file))
  
  cat("", join(paste0('var ', symbolnames, '>=0;'), '\n'), 
      "minimize objective : ", join(symbolnames, '+'), ";", 
      join(paste0('V', 1:len(views), ':', views, ';'), '\n'), "", 
      file = model.file, append = FALSE)
  
  model = Rglpk.read.model(model.file,type = 'MathProg')
  temp = Rglpk.create.constraints(model)$constraints
  A = t(as.matrix(temp$A))
  colnames(A) = symbolnames
  list(A = A, b = temp$b, meq = temp$meq)
  }

min.var.portfolio.gmpl <- function(ia, constraints) {
  load.packages('quadprog, corpcor')
  cov.temp = ia$cov
  n0 = ia$n
  n = nrow(constraints$A)
  
  if(n != nrow(cov.temp)) {
    temp =  matrix(0, n, n)
    temp[1:n0, 1:n0] = cov.temp[1:n0, 1:n0]
    cov.temp = temp
  }
  
  if(!is.positive.definite(cov.temp)) {
    cov.temp <- make.positive.definite(cov.temp, 0.000000001)
  }
  
  binary.vec = 0
  if(!is.null(constraints$binary.index)) 
    binary.vec = constraints$binary.index
  sol = solve.QP.bounds(Dmat = cov.temp, dvec = rep(0, nrow(cov.temp)) , 
                        Amat = constraints$A, bvec = constraints$b, 
                        constraints$meq, lb = constraints$lb, 
                        ub = constraints$ub, binary.vec = binary.vec)
  
  if(binary.vec[1] != 0) cat(sol$counter,'QP calls made to solve problem 
                             with', len(binary.vec), 'binary variables using 
                             Branch&Bound', '\n')
  x = sol$solution[1:ia$n]
  names(x) = ia$symbols
  return(x)
  }

portopt.mathprog.test <- function() {
  load.packages('quantmod,quadprog,corpcor')
  tickers = dow.jones.components()
  ia = aa.test.create.ia.custom(tickers, dates = '2000::2010')
  n = ia$n
  constraints = new.constraints(n, lb = 0, ub = 1)
  constraints = add.constraints(rep(1, n), 1, type = '=', constraints)
  x = min.var.portfolio.gmpl(ia, constraints)
  png(filename = 'plot1.png', width = 600, height = 400, units = 'px', 
      pointsize = 12, bg = 'white')
  barplot(100 * x, las = 2, main = 'Minimum Variance Portfolio')
  dev.off()
  load.packages('Rglpk')
  model.file = 'model1.mod'
  
  cat('
set SYMBOLS ;
var weight{i in SYMBOLS} >= 0, <= 1 ;
minimize alpha : sum{i in SYMBOLS} weight[i] ;
subject to fully_invested : sum{i in SYMBOLS} weight[i] = 1 ;
data;
set SYMBOLS := ', ia$symbols, ';', file = model.file, append = FALSE)

model = Rglpk.read.model(model.file, type = 'MathProg')
constraints = Rglpk.create.constraints(model)$constraints
x = min.var.portfolio.gmpl(ia, constraints)
png(filename = 'plot2.png', width = 600, height = 400, units = 'px', 
    pointsize = 12, bg = 'white')
barplot(100 * x, las = 2, main='Minimum Variance Portfolio using GNU MathProg 
        model')
dev.off()
model.file = 'model2.mod'
cat('
set SYMBOLS ;
var weight{i in SYMBOLS} >= 0, <= 1 ;
var held{SYMBOLS} binary;
minimize alpha : sum{i in SYMBOLS} weight[i] ;
subject to fully_invested : sum{i in SYMBOLS} weight[i] = 1 ;
subject to MinWgt {i in SYMBOLS} : weight[i] >= 0.025 * held[i];
subject to MaxWgt {i in SYMBOLS} : weight[i] <= .20 * held[i] ;
subject to MaxAssetsLB : 0 <= sum {i in SYMBOLS} held[i] ;
subject to MaxAssetsUB : sum {i in SYMBOLS} held[i] <= 6 ;
data;
set SYMBOLS := ', ia$symbols, ';
', file = model.file, append = FALSE)
model = Rglpk.read.model(model.file,type = 'MathProg')
constraints = Rglpk.create.constraints(model)$constraints
x = min.var.portfolio.gmpl(ia, constraints)
png(filename = 'plot3.png', width = 600, height = 400, units = 'px', 
    pointsize = 12, bg = 'white')
barplot(100 * x, las = 2, main = 'Minimum Variance Portfolio using GNU 
        MathProg model \n with Minimum Investment and Number of Assets 
        Constraints')
dev.off()
model.file = 'model3.mod'
cat('
set SYMBOLS ;
var long {i in SYMBOLS} >= 0, <= 0.8 ;
var short{i in SYMBOLS} >= 0, <= 0.5 ;
var islong{SYMBOLS} binary;
minimize alpha : sum{i in SYMBOLS} long[i] ;
subject to fully_invested : sum{i in SYMBOLS} (long[i] - short[i]) = 1 ;
subject to leverage : sum{i in SYMBOLS} (long[i] + short[i]) = 1.6 ;
subject to long_flag  {i in SYMBOLS} : long[i]  <= islong[i] ;
subject to short_flag {i in SYMBOLS} : short[i] <= (1 - islong[i]) ;
data;
set SYMBOLS := ', ia$symbols, ';', file = model.file, append = FALSE)
model = Rglpk.read.model(model.file,type = 'MathProg')
constraints = Rglpk.create.constraints(model)$constraints
x = min.var.portfolio.gmpl(aa.test.ia.add.short(ia), constraints)
x = x[1:ia$n] - x[-c(1:ia$n)]
png(filename = 'plot4.png', width = 600, height = 400, units = 'px', 
    pointsize = 12, bg = 'white')
barplot(100 * x, las = 2, main = 'Minimum Variance Portfolio using GNU 
        MathProg model \n with 130:30 Constraints')
dev.off()
ia = aa.test.create.ia.custom(tickers[1:15], dates = '2000::2010')
model.file = 'model4.mod'
param = ia$cov[, 1, drop=F]
colnames(param) = 'CurWgt'
param[, 'CurWgt'] = 1/ia$n
cat('
set SYMBOLS ;
param CurWgt{SYMBOLS} ;
var weight{i in SYMBOLS} >= 0, <= 1 ;
var TradePos{i in SYMBOLS} >= 0 ;
var TradeNeg{i in SYMBOLS} >= 0 ;
var TradeFlag{SYMBOLS} binary;
var trade{SYMBOLS} binary;
minimize alpha : sum{i in SYMBOLS} weight[i] ;
subject to fully_invested : sum{i in SYMBOLS} weight[i] = 1 ;
subject to TradeRange {i in SYMBOLS} : (CurWgt[i] - weight[i]) = (TradePos[i] - TradeNeg[i]) ;
subject to TradeFlagPos {i in SYMBOLS} : TradePos[i] <= 100 * TradeFlag[i];
subject to TradeFlagNeg {i in SYMBOLS} : TradeNeg[i] <= 100 * (1 - TradeFlag[i]);
subject to MinTradeSize {i in SYMBOLS} : (TradePos[i] + TradeNeg[i]) >= 0.01 * trade[i];
subject to MaxTradeSize {i in SYMBOLS} : (TradePos[i] + TradeNeg[i]) <= .90 * trade[i] ;
subject to MaxTrade : sum {i in SYMBOLS} trade[i] <= 48 ;
data;
set SYMBOLS := ', ia$symbols, ';
param : CurWgt:=', file = model.file, append = FALSE)
write.table(param, sep = '\t', quote = F, col.names = F, file = model.file, 
            append = TRUE)
cat(';', file = model.file, append = TRUE)
model = Rglpk.read.model(model.file, type = 'MathProg')
constraints = Rglpk.create.constraints(model)$constraints
x = min.var.portfolio.gmpl(ia, constraints)
sqrt(x %*% ia$cov %*% x)
png(filename = 'plot5.png', width = 600, height = 400, units = 'px', 
    pointsize = 12, bg = 'white')
barplot(100 * x, las = 2, main = 'Minimum Variance Portfolio using GNU 
        MathProg model \n with Turnover Constraints')
dev.off()
ia = aa.test.create.ia.custom(tickers[1:10], dates = '2000::2010')
model.file = 'model4.mod'
param = ia$cov[, 1, drop = F]
colnames(param) = 'CurWgt'
param[, 'CurWgt'] = 1/ia$n
cat('
set SYMBOLS ;
param CurWgt{SYMBOLS} ;
var weight{i in SYMBOLS} >= 0, <= 1 ;
var TradePos{i in SYMBOLS} >= 0 ;
var TradeNeg{i in SYMBOLS} >= 0 ;
var TradeFlag{SYMBOLS} binary;
var trade{SYMBOLS} binary;
minimize alpha : sum{i in SYMBOLS} weight[i] ;
subject to fully_invested : sum{i in SYMBOLS} weight[i] = 1 ;
subject to TradeRange {i in SYMBOLS} : (CurWgt[i] - weight[i]) = (TradePos[i] - TradeNeg[i]) ;
subject to TradeFlagPos {i in SYMBOLS} : TradePos[i] <= 100 * TradeFlag[i];
subject to TradeFlagNeg {i in SYMBOLS} : TradeNeg[i] <= 100 * (1 - TradeFlag[i]);
subject to MinTradeSize {i in SYMBOLS} : (TradePos[i] + TradeNeg[i]) >= 0.05 * trade[i];
subject to MaxTradeSize {i in SYMBOLS} : (TradePos[i] + TradeNeg[i]) <= .20 * trade[i] ;
subject to MaxTrade : sum {i in SYMBOLS} trade[i] <= 8 ;
data;
set SYMBOLS := ', ia$symbols, ';
param : CurWgt:=', file = model.file, append = FALSE)

write.table(param, sep = '\t', quote = F, col.names = F, file = model.file, 
            append = TRUE)
cat(';', file = model.file, append = TRUE)
model = Rglpk.read.model(model.file,type = 'MathProg')
constraints = Rglpk.create.constraints(model)$constraints
x = min.var.portfolio.gmpl(ia, constraints)
sqrt(x %*% ia$cov %*% x)
png(filename = 'plot6.png', width = 600, height = 400, units = 'px', 
    pointsize = 12, bg = 'white')
barplot(100 * x, las = 2, main = 'Minimum Variance Portfolio using GNU 
        MathProg model \n with Turnover Constraints')
dev.off()
}

add.constraint.omega <- function(ia, value, type = c('=', '>=', '<='), 
                                 constraints) {
  if(is.null(ia$parameters.omega)) omega = 0 else 
    omega = ia$parameters.omega
  
  n0 = ncol(ia$hist.returns)
  n = nrow(constraints$A)
  nt = nrow(ia$hist.returns)
  constraints = add.variables(2 * nt + 1, constraints, lb = c(rep(0, 2 * nt), 
                                                              -Inf))
  constraints$A[n + 2 * nt + 1, ] = -constraints$b
  constraints$b[] = 0
  index = which(constraints$ub[1:n] < +Inf)
  if(len(index) > 0) {
    a = rbind(diag(n), matrix(0, 2 * nt, n), -constraints$ub[1:n])
    constraints = add.constraints(a[, index], rep(0, len(index)), '<=', 
                                  constraints)
    }
  index = which(constraints$lb[1:n] > -Inf)
  
  if( len(index) > 0 ) {
    a = rbind(diag(n), matrix(0, 2 * nt, n), -constraints$lb[1:n])
    constraints = add.constraints(a[, index], rep(0, len(index)), '>=', 
                                  constraints)
    }
  
  constraints$lb[1:n] = -Inf
  constraints$ub[1:n] = Inf
  a = rbind(matrix(0, n, nt), -diag(nt), diag(nt), -omega)
  a[1:n0, ] = t(ia$hist.returns)
  constraints = add.constraints(a, rep(0, nt), '=', constraints)
  constraints = add.constraints(c( rep(0, n), rep(0, nt), (1/nt) * 
                                     rep(1, nt), 0), 1, '=', constraints)
  constraints = add.constraints(c(rep(0, n), (1/nt) * rep(1, nt), 
                                  rep(0, nt), 0), value, type[1], constraints)
  return(constraints)
  }

portfolio.omega <- function(weight, ia) {
  weight = weight[, 1:ia$n]
  
  if(is.null(ia$parameters.omega)) omega = 0 else 
    omega = ia$parameters.omega
  
  portfolio.returns = weight %*% t(ia$hist.returns)
  return(apply(portfolio.returns, 1, function(x) 
    mean(pmax(x - omega,0)) / mean(pmax(omega - x,0))))
  }

max.omega.portfolio <- function(ia, constraints, 
                                type = c('mixed', 'lp', 'nlp')) {
  n = nrow(constraints$A)
  nt = nrow(ia$hist.returns)
  type = type[1]
  
  if(type == 'mixed'	|| type == 'lp') {
    sol = optimize.portfolio(ia, constraints, add.constraint.omega, 
                             portfolio.omega, 'max', T)
    x = rep(NA, n)
    
    if(!inherits(sol, 'try-error') && sol$status ==0) {
      x0 = sol$solution[1:n]
      u = sol$solution[(1 + n):(n + nt)]
      d = sol$solution[(n + nt + 1):(n + 2 * nt)]
      t = sol$solution[(n + 2 * nt + 1):(n + 2 * nt + 1)]
      x = x0/t
      }
    }
  
  if((type == 'mixed' && (sol$status != 0 || any(u * d != 0))) || 
     type == 'nlp') {
    if(is.null(ia$parameters.omega)) omega = 0 else 
      omega = ia$parameters.omega
    
    fn <- function(x) {
      portfolio.returns = x %*% t(ia$hist.returns)
      mean(pmax(portfolio.returns - omega, 0)) / 
        mean(pmax(omega - portfolio.returns, 0))
      }
    
    x = optimize.portfolio.nlp(ia, constraints, fn, direction = 'max')
    }
  
  return(x)
  }

portopt.omega <- function(ia, constraints = NULL, nportfolios = 50, 
                          name = 'Omega') {
  out = list(weight = matrix(NA, nportfolios, nrow(constraints$A)))
  colnames(out$weight) = rep('', ncol(out$weight))
  colnames(out$weight)[1:ia$n] = ia$symbols
  ef.risk = portopt(ia, constraints, 2)
  out$weight[nportfolios, ] = ef.risk$weight[2, ]
  out$weight[1, ] = ef.risk$weight[1, ]
  constraints$x0 = out$weight[1, ]
  out$return = portfolio.return(out$weight, ia)
  target = seq(out$return[1], out$return[nportfolios], 
               length.out = nportfolios)
  constraints = add.constraints(c(ia$expected.return, 
                                  rep(0, nrow(constraints$A) - ia$n)), 
                                target[1], type = '<=', constraints)
  
  for(i in 1:nportfolios) {
    cat('i =', i, '\n')
    constraints$b[len(constraints$b)] = -target[i]
    out$weight[i, ] = max.omega.portfolio(ia, constraints)
    constraints$x0 = out$weight[i, ]
    }
  
  out$return = portfolio.return(out$weight, ia)
  out$risk = portfolio.risk(out$weight, ia)
  out$name = name
  return(out)
  }

plot.omega <- function(weight, ia) {
  omegafn = function(x,L) { mean(pmax(x - L,0)) / mean(pmax(L - x,0)) }
  
  if(is.null(ia$parameters.omega)) omega = 0 else 
    omega = ia$parameters.omega
  
  weight = weight[, 1:ia$n, drop = F]
  portfolio.returns = weight %*% t(ia$hist.returns)
  threshhold = quantile(portfolio.returns, probs = c(0.05, 0.95))
  threshhold = seq(threshhold[1], threshhold[2], length.out = 100)
  par(mar = c(4 ,4, 2, 1), cex = 0.8)
  
  for(i in 1:nrow(weight)) {
    data = sapply(threshhold, function(L) omegafn(portfolio.returns[i, ], L))
    if(i==1) plot(threshhold, log(data), type = 'l', col = i, 
                  xlab = 'Threshhold', ylab='Log(Omega)', 
                  main='Portfolio Omega')
    lines(threshhold, log(data), col = i)
    }
  
  abline(v = omega, col = 'orange')
  grid()
  plota.legend(rownames(weight), 1:nrow(weight), x = 'bottomleft')
  }

new.constraints <- function(n, A = NULL, b = NULL, type = c('=', '>=', '<='), 
                            lb = NA, ub = NA) {
  meq = 0
  if(is.null(A) || is.na(A) || is.null(b) || is.na(b)) {
    A = matrix(0, n, 0)
    b = c()
  } else {
    if(is.null(dim(A))) dim(A) = c(len(A), 1)
    if(type[1] == '=') meq = len(b)
    if(type[1] == '<=') {
      A = -A
      b = -b
    }
  }
  
  if(is.null(lb) || is.na(lb)) lb = rep(NA, n)
  if(len(lb) != n) lb = rep(lb[1], n)
  if(is.null(ub) || is.na(ub)) ub = rep(NA, n)
  if(len(ub) != n) ub = rep(ub[1], n)
  return(list(n = n, A = A, b = b, meq = meq, lb = lb, ub = ub))
  }

add.constraints <- function(A, b, type = c('=', '>=', '<='), constraints) {
  if(is.null(constraints)) constraints = new.constraints(n = nrow(A))
  if(is.null(dim(A))) A = matrix(A)
  if(len(b) == 1) b = rep(b, ncol(A))
  
  if(type[1] == '=') {
    constraints$A = cbind(A, constraints$A)
    constraints$b = c(b, constraints$b)
    constraints$meq = constraints$meq + len(b)
    }
  
  if(type[1] == '>=') {
    constraints$A = cbind(constraints$A, A)
    constraints$b = c(constraints$b, b)
    }
  
  if(type[1] == '<=') {
    constraints$A = cbind(constraints$A, -A)
    constraints$b = c(constraints$b, -b)
    }
  
  return(constraints)
  }

add.variables <- function(n, constraints, lb = NA, ub = NA) {
  constraints$A = rbind(constraints$A, matrix(0, n, len(constraints$b)))
  if(is.null(lb) || is.na(lb)) lb = rep(NA, n)
  if(len(lb) != n) lb = rep(lb[1], n)
  if(is.null(ub) || is.na(ub)) ub = rep(NA, n)
  if(len(ub) != n) ub = rep(ub[1], n)
  constraints$lb = c(constraints$lb, lb)
  constraints$ub = c(constraints$ub, ub)
  constraints$n = constraints$n + n
  return(constraints)
  }

delete.constraints <- function(delete.index, constraints) {
  constraints$A = constraints$A[, -delete.index, drop=F]
  constraints$b = constraints$b[-delete.index]
  constraints$meq = constraints$meq - len(intersect((1:constraints$meq), 
                                                    delete.index))
  return(constraints)
  }

type.constraints <- function(constraints) {
  c(rep('=', constraints$meq), rep('>=', len(constraints$b) - 
                                     constraints$meq))
  }

create.basic.constraints <- function(n, const.lb = 0, const.ub = 1, 
                                     const.sum = 1) {
  if(len(const.lb) == 1) const.lb = rep(const.lb, n)
  if(len(const.ub) == 1) const.ub = rep(const.ub, n)
  
  constraints = new.constraints(n, lb = const.lb, ub = const.ub)
  constraints = add.constraints(diag(n), type='>=', b = const.lb, constraints)
  constraints = add.constraints(diag(n), type='<=', b=const.ub, constraints)
  
  if(!is.na(const.sum))
    constraints = add.constraints(rep(1, n), type = '=', b = const.sum, 
                                  constraints)
  return(constraints)
  }

merge.constraints <- function(constraints1, constraints2) {
  if(constraints1$n != constraints2$n) 
    stop('merge.constraints: both constraints must be based on same number of 
         assets')
  
  if(constraints2$meq > 0) {
    constraints1$A = cbind( constraints2$A[, 1:constraints2$meq], 
                            constraints1$A, 
                            constraints2$A[, -c(1:constraints2$meq)])
    constraints1$b = c(constraints2$b[1:constraints2$meq], constraints1$b, 
                       constraints2$b[-c(1:constraints2$meq)])
    constraints1$meq = constraints1$meq + constraints2$meq
  } else {
    constraints1$A = cbind( constraints1$A, constraints2$A)
    constraints1$b = c(constraints1$b, constraints2$b)
    }
  
  constraints1$lb =  pmax(constraints1$lb, constraints2$lb, na.rm = T)
  constraints1$ub =  pmin(constraints1$ub, constraints2$ub, na.rm = T)
  constraints1
  }

min.portfolio <- function(ia, constraints, add.constraint.fn, min.risk.fn) {
  optimize.portfolio(ia, constraints, add.constraint.fn, min.risk.fn)
  }

optimize.portfolio <- function(ia, constraints, add.constraint.fn, 
                               min.risk.fn, direction = 'min', 
                               full.solution = F) {
  load.packages('quadprog,corpcor,lpSolve,kernlab')
  n = nrow(constraints$A)
  nt = nrow(ia$hist.returns)
  constraints = match.fun(add.constraint.fn)(ia, 0, '>=', constraints)
  f.obj = constraints$A[, ncol(constraints$A)]
  constraints = delete.constraints(ncol(constraints$A), constraints)
  f.con = constraints$A
  f.dir = c(rep('=', constraints$meq), rep('>=', len(constraints$b) 
                                           - constraints$meq))
  f.rhs = constraints$b
  x = NA
  binary.vec = 0
  
  if(!is.null(constraints$binary.index)) 
    binary.vec = constraints$binary.index
  sol = try(solve.LP.bounds(direction, f.obj, t(f.con), f.dir, f.rhs, 
                            lb = constraints$lb, ub = constraints$ub, 
                            binary.vec = binary.vec, default.lb = -100), TRUE)
  
  if(!inherits(sol, 'try-error')) {
    x = sol$solution[1:n]
    if(F) {
      f.obj %*% sol$solution  - match.fun(min.risk.fn)(t(x), ia)
      }
    }
  
  if(full.solution) x = sol
  return(x)
  }

optimize.portfolio.nlp <- function(ia, constraints, fn, nl.constraints = NULL, 
                                   direction = 'min', full.solution = F) {
  load.packages('Rdonlp2', repos ='http://R-Forge.R-project.org')
  if(direction == 'min') fnscale = 1 else fnscale = -1
  cntl = donlp2Control()
  cntl$silent = T
  cntl$fnscale = fnscale
  cntl$iterma =10000
  cntl$nstep = 100
  cntl$epsx = 1e-10
  par.l = constraints$lb
  par.u = constraints$ub
  p = rep(1, nrow(constraints$A))
  
  if(!is.null(constraints$x0)) {
    if(sum(is.na(constraints$x0)) == 0) p = constraints$x0
    }
  
  A = t(constraints$A)
  lin.l = constraints$b
  lin.u = constraints$b
  lin.u[-c(1:constraints$meq)] = +Inf
  x = NA
  if(!is.null(nl.constraints)) {
    sol = donlp2(p, fn, par.lower = par.l, par.upper = par.u, A = A, 
                 lin.u = lin.u, lin.l = lin.l, control = cntl, 
                 nlin = nl.constraints$constraints, 
                 nlin.upper = nl.constraints$upper, 
                 nlin.lower = nl.constraints$lower)
  } else {
    sol = donlp2(p, fn, par.lower = par.l, par.upper = par.u, A = A, 
                 lin.u = lin.u, lin.l = lin.l, control = cntl)
  }
  
  if(!inherits(sol, 'try-error')) {
    x = sol$par
    }
  if(full.solution) x = sol
  return(x)
  }

add.constraint.maxloss <- function(ia, value, type = c('=', '>=', '<='), 
                                   constraints) {
  n0 = ncol(ia$hist.returns)
  n = nrow(constraints$A)
  nt = nrow(ia$hist.returns)
  constraints = add.variables(1, constraints)
  a = rbind(matrix(0, n, nt), 1)
  a[1:n0, ] = t(ia$hist.returns)
  constraints = add.constraints(a, rep(0, nt), '>=', constraints)
  constraints = add.constraints(c(rep(0, n), 1), value, type[1], constraints)
  return(constraints)
  }

portfolio.maxloss <- function(weight, ia) {
  weight = weight[, 1:ia$n, drop = F]
  portfolio.returns = weight %*% t(ia$hist.returns)
  return(-apply(portfolio.returns, 1, min))
  }

min.maxloss.portfolio <- function(ia, constraints) {
  min.portfolio(ia, constraints, add.constraint.maxloss, portfolio.maxloss)
  }

add.constraint.mad <- function(ia, value, type = c('=', '>=', '<='), 
                               constraints) {
  n0 = ncol(ia$hist.returns)
  n = nrow(constraints$A)
  nt = nrow(ia$hist.returns)
  constraints = add.variables(2 * nt, constraints, lb = 0)
  a = rbind(matrix(0, n, nt), -diag(nt), diag(nt))
  a[1:n0, ] = t(ia$hist.returns) - repmat(colMeans(ia$hist.returns), 1, nt)
  constraints = add.constraints(a, rep(0, nt), '=', constraints)
  constraints = add.constraints(c(rep(0, n), (1/nt) * rep(1, 2 * nt)), 
                                value, type[1], constraints)
  return(constraints)
  }

portfolio.mad <- function(weight, ia) {
  if(is.null(dim(weight))) dim(weight) = c(1, len(weight))
  weight = weight[, 1:ia$n, drop = F]
  portfolio.returns = weight %*% t(ia$hist.returns)
  return(apply(portfolio.returns, 1, function(x) mean(abs(x - mean(x)))))
  }

min.mad.portfolio <- function(ia, constraints) {
  min.portfolio(ia, constraints, add.constraint.mad, portfolio.mad)
  }

add.constraint.cvar <- function(ia, value, type = c('=', '>=', '<='), 
                                constraints) {
  if(is.null(ia$parameters.alpha)) alpha = 0.95 else 
    alpha = ia$parameters.alpha
  
  n0 = ncol(ia$hist.returns)
  n = nrow(constraints$A)
  nt = nrow(ia$hist.returns)
  constraints = add.variables(nt + 1, constraints, lb = c(rep(0, nt), -Inf))
  a = rbind( matrix(0, n, nt), diag(nt), 1)
  a[1:n0, ] = t(ia$hist.returns)
  constraints = add.constraints(a, rep(0, nt), '>=', constraints)
  constraints = add.constraints(c(rep(0, n), (1/(1 - alpha)) * (1/nt) * 
                                    rep(1, nt), 1), value, type[1], 
                                constraints)
  return(constraints)
  }

portfolio.cvar <- function(weight, ia) {
  weight = weight[, 1:ia$n, drop = F]
  if(is.null(ia$parameters.alpha)) alpha = 0.95 else 
    alpha = ia$parameters.alpha
  
  portfolio.returns = weight %*% t(ia$hist.returns)
  return(apply(portfolio.returns, 1, function(x) -compute.cvar(x, 1 - alpha)))
  }

min.cvar.portfolio <- function(ia, constraints) {
  min.portfolio(ia, constraints, add.constraint.cvar, portfolio.cvar)
  }

portfolio.var <- function(weight, ia) {
  weight = weight[, 1:ia$n, drop = F]
  if(is.null(ia$parameters.alpha)) alpha = 0.95 else 
    alpha = ia$parameters.alpha
  
  portfolio.returns = weight %*% t(ia$hist.returns)
  return(apply(portfolio.returns, 1, function(x) -compute.var(x, 1 - alpha)))
  }

add.constraint.cdar <- function(ia, value, type = c('=', '>=', '<='), 
                                constraints) {
  if(is.null(ia$parameters.alpha)) alpha = 0.95 else 
    alpha = ia$parameters.alpha
  
  n0 = ncol(ia$hist.returns)
  n = nrow(constraints$A)
  nt = nrow(ia$hist.returns)
  constraints = add.variables(2 * nt + 1, constraints, 
                              lb = c(rep(0, nt), rep(-Inf, nt + 1)))
  a = rbind(matrix(0, n, nt), diag(nt), 1, -diag(nt))
  a[1:n0, ] = t(apply(t(ia$hist.returns), 1, cumsum))
  constraints = add.constraints(a, rep(0, nt), '>=', constraints)
  a = rbind( matrix(0, n, nt), 0 * diag(nt), 0, diag(nt))
  a[1:n0, ] = -t(apply(t(ia$hist.returns), 1, cumsum))
  constraints = add.constraints(a, rep(0, nt), '>=', constraints)
  temp = diag(nt);
  temp[-nt, -1] = -diag((nt-1))
  diag(temp) = 1
  a = rbind(matrix(0, n, nt), 0 * diag(nt), 0, temp)
  a = a[, -1]
  constraints = add.constraints(a, rep(0, (nt - 1)), '>=', constraints)
  constraints = add.constraints(c(rep(0, n), (1/(1 - alpha)) * (1/nt) * 
                                    rep(1, nt), 1, rep(0, nt)), value, 
                                type[1], constraints)
  return(constraints)
  }

portfolio.cdar <- function(weight, ia) {
  weight = weight[, 1:ia$n, drop = F]
  if(is.null(ia$parameters.alpha)) alpha = 0.95 else 
    alpha = ia$parameters.alpha
  
  portfolio.returns = weight %*% t(ia$hist.returns)
  return(apply(portfolio.returns, 1, function(x) {
    x = cumsum(x)
    x = x - cummax(x) -compute.cvar(x, 1-alpha)
    }))
  }

min.cdar.portfolio <- function(ia, constraints) {
  min.portfolio(ia, constraints, add.constraint.cdar, portfolio.cdar)
  }

portfolio.cdar.real <- function(weight, ia) {
  weight = weight[, 1:ia$n, drop = F]
  if(is.null(ia$parameters.alpha)) alpha = 0.95 else 
    alpha = ia$parameters.alpha
  
  portfolio.returns = weight %*% t(ia$hist.returns)
  out = rep(0, nrow(weight))
  
  for(i in 1:nrow(weight) ) {
    portfolio.equity = cumprod(1 + portfolio.returns[i, ])
    x = compute.drawdowns(portfolio.equity)
    out[i] = compute.cvar(x, alpha)
    }
  return(out)
  }

compute.drawdowns <- function(portfolio.equity, make.plot = FALSE) {
  temp = portfolio.equity / cummax(portfolio.equity) - 1
  temp = c(temp, 0)
  drawdown.start = which(temp == 0 & mlag(temp, -1) != 0)
  drawdown.end = which(temp == 0 & mlag(temp, 1) != 0)
  if(make.plot) {
    plot((1:len(temp)), temp, type = 'l')
    points((1:len(temp))[drawdown.start] , temp[drawdown.start], col = 'red')
    points((1:len(temp))[drawdown.end] , temp[drawdown.end], col = 'blue')
    }
  
  return(apply(cbind(drawdown.start, drawdown.end), 1, function(x){ 
    min(temp[x[1]:x[2]], na.rm = T)})
    )
  }

min.avgcor.portfolio <- function(ia, constraints) {
  cov = ia$cov[1:ia$n, 1:ia$n]
  s = sqrt(diag(cov))
  fn <- function(x){
    sd_x = sqrt( t(x) %*% cov %*% x )
    mean((x %*% cov) / (s * sd_x))
    }
  
  x = optimize.portfolio.nlp(ia, constraints, fn)
  return(x)
  }

portfolio.avgcor <- function(weight, ia) {
  weight = weight[, 1:ia$n, drop = F]
  cov = ia$cov[1:ia$n, 1:ia$n]
  s = sqrt(diag(cov))
  return(apply(weight, 1, function(x) {
    sd_x = sqrt(t(x) %*% cov %*% x)
    mean((x %*% cov) / (s * sd_x))
    }))
  }

min.cor.insteadof.cov.portfolio <- function(ia, constraints) {
  if(is.null(ia$cov.temp)) ia$cov.temp = ia$cov
  sol = solve.QP.bounds(Dmat = ia$correlation, 
                        dvec = rep(0, nrow(ia$cov.temp)) , 
                        Amat = constraints$A, bvec = constraints$b, 
                        constraints$meq, lb = constraints$lb, 
                        ub = constraints$ub)
  return(sol$solution)
  }

portfolio.avgcor.real <- function(weight, ia) {
  weight = weight[, 1:ia$n, drop = F]
  portfolio.returns = weight %*% t(ia$hist.returns)
  return(apply(portfolio.returns, 1, function(x) 
    mean(cor(ia$hist.returns, x))))
  }

add.constraint.mad.downside <- function(ia, value, type = c('=', '>=', '<='), 
                                        constraints) {
  n0 = ncol(ia$hist.returns)
  n = nrow(constraints$A)
  nt = nrow(ia$hist.returns)
  constraints = add.variables(nt, constraints, lb = 0)
  a = rbind(matrix(0, n, nt), diag(nt))
  
  if(is.null(ia$parameters.mar) || is.na(ia$parameters.mar)) {
    a[1:n0, ] = t(ia$hist.returns) - repmat(colMeans(ia$hist.returns), 1, nt)
    constraints = add.constraints(a, rep(0, nt), '>=', constraints)
  } else {
    a[1:n0, ] = t(ia$hist.returns)
    constraints = add.constraints(a, rep(ia$parameters.mar, nt), '>=', 
                                  constraints)
  }
  
  constraints = add.constraints(c(rep(0, n), (1/nt) * rep(1, nt)), value, 
                                type[1], constraints)
  return(constraints)
  }

portfolio.mad.downside <- function(weight, ia) {
  weight = weight[, 1:ia$n, drop = F]
  portfolio.returns = weight %*% t(ia$hist.returns)
  
  if(is.null(ia$parameters.mar) || is.na(ia$parameters.mar)) {
    return(apply(portfolio.returns, 1, function(x) 
      mean(pmax(mean(x) - x, 0))))
  } else {
    return(apply(portfolio.returns, 1, function(x) 
      mean(pmax(ia$parameters.mar - x, 0))))
    }
  }

min.mad.downside.portfolio <- function(ia, constraints) {
  min.portfolio(ia, constraints, add.constraint.mad.downside, 
                portfolio.mad.downside)
  }

portfolio.risk.downside <- function(weight, ia) {
  weight = weight[, 1:ia$n, drop = F]
  portfolio.returns = weight %*% t(ia$hist.returns)
  
  if(is.null(ia$parameters.mar) || is.na(ia$parameters.mar)) {
    return(apply(portfolio.returns, 1, function(x) 
      sqrt(mean(pmax(mean(x) - x, 0)^2))))
  } else {
    return(apply(portfolio.returns, 1, function(x) 
      sqrt(mean(pmax(ia$parameters.mar - x, 0)^2))))
    }
  }

min.risk.downside.portfolio <- function(ia, constraints) {
  n = nrow(constraints$A)
  nt = nrow(ia$hist.returns)
  
  constraints = add.constraint.mad.downside(ia, 0, '>=', constraints)
  f.obj = constraints$A[, ncol(constraints$A)]
  constraints = delete.constraints(ncol(constraints$A), constraints)
  Dmat = diag(len(f.obj))
  diag(Dmat) = f.obj
  
  if(!is.positive.definite(Dmat)) {
    Dmat <- make.positive.definite(Dmat)
    }
  x = NA
  binary.vec = 0
  
  if(!is.null(constraints$binary.index)) 
    binary.vec = constraints$binary.index
  
  sol = try(solve.QP.bounds(Dmat = Dmat, dvec = rep(0, nrow(Dmat)), 
                            Amat = constraints$A, bvec = constraints$b, 
                            constraints$meq, lb = constraints$lb, 
                            ub = constraints$ub, binary.vec = binary.vec), 
            TRUE)
  
  if(!inherits(sol, 'try-error')) {
    x = sol$solution[1:n]
    if(F) {
      sol$solution %*% Dmat %*% (sol$solution) - 
        portfolio.risk.downside(t(x), ia)^2
      }
    }
  return(x)
  }

add.constraint.gini <- function(ia, value, type = c('=', '>=', '<='), 
                                constraints) {
  n0 = ncol(ia$hist.returns)
  n = nrow(constraints$A)
  nt = nrow(ia$hist.returns)
  constraints = add.variables(nt * (nt - 1), constraints, lb = 0)
  a = matrix(0, n0 + nt * (nt - 1), nt * (nt - 1)/2)
  diag(a[(n0 + 1):(n0 + nt * (nt - 1)/2), ]) = -1
  diag(a[(n0 + 1 + nt * (nt - 1)/2) : (n0 + nt * (nt - 1)), ]) = 1
  hist.returns = as.matrix(ia$hist.returns)
  i.start = 0
  
  for(t in 1:(nt-1)) {
    index = (i.start+1):(i.start + nt - t)
    for(i in 1:n0) {
      a[i, index] = (hist.returns[t, i] - hist.returns[, i])[(t + 1):nt]
      }
    i.start = i.start + nt - t
    }
  
  constraints = add.constraints(a, 0, '=', constraints)
  constraints = add.constraints(c(rep(0, n), rep(1, nt * (nt - 1))), value, 
                                type[1], constraints)
  return(constraints)
  }

min.gini.portfolio <- function(ia, constraints) {
  min.portfolio(ia, constraints, add.constraint.gini, 
                portfolio.gini.coefficient)
  }

portfolio.gini.coefficient <- function(weight, ia) {
  if(is.null(dim(weight))) dim(weight) = c(1, len(weight))
  weight = weight[, 1:ia$n, drop = F]
  portfolio.returns = weight %*% t(ia$hist.returns)
  n = ncol(portfolio.returns)
  one.to.n = 1:n
  out = weight[, 1] * NA
  out[] = apply(portfolio.returns, 1, function(x) {
    temp = sort(x, decreasing = F)
    sum((2 * one.to.n - n - 1) * temp)
    })
  
  out = 2 * out /(n * (n - 1))
  return(out)
  }

lp.obj.portfolio <- function(ia, constraints, 
                             f.obj = c(ia$expected.return, 
                                       rep(0, nrow(constraints$A) - ia$n)), 
                             direction = 'min') {
  x = NA
  binary.vec = 0
  
  if(!is.null(constraints$binary.index)) 
    binary.vec = constraints$binary.index
  
  sol = try(solve.LP.bounds(direction, f.obj, t(constraints$A), 
                            c(rep('=', constraints$meq), 
                              rep('>=', len(constraints$b) - 
                                    constraints$meq)), constraints$b, 
                            lb = constraints$lb, ub = constraints$ub, 
                            binary.vec = binary.vec), TRUE)
  
  if(!inherits(sol, 'try-error')) x = sol$solution
  
  return(x)
  }

max.return.portfolio <- function(ia, constraints) {
  lp.obj.portfolio(ia, constraints, direction = 'max')
  }

portfolio.return <- function(weight, ia) {
  if(is.null(dim(weight))) dim(weight) = c(1, len(weight))
  weight = weight[, 1:ia$n, drop = F]
  portfolio.return = weight %*% ia$expected.return
  return(portfolio.return)
  }

portfolio.geometric.return <- function(weight, ia) {
  weight = weight[, 1:ia$n, drop = F]
  portfolio.returns = weight %*% t(ia$hist.returns)
  return( apply(portfolio.returns, 1, function(x) 
    (prod(1 + x)^(1 / len(x)))^ia$annual.factor - 1 ))
  }

max.geometric.return.portfolio <- function(ia, constraints, min.risk, 
                                           max.risk) {
  fn <- function(x){
    portfolio.returns = x %*% t(ia$hist.returns)
    prod(1 + portfolio.returns)
    }
  
  nlcon1 <- function(x){
    sqrt(t(x) %*% ia$cov %*% x)
    }
  
  nl.constraints = list()
  nl.constraints$constraints = list(nlcon1)
  nl.constraints$upper = c(max.risk)
  nl.constraints$lower = c(min.risk)
  x = optimize.portfolio.nlp(ia, constraints, fn, nl.constraints, 
                             direction = 'max')
  return(x)
  }

portfolio.unrebalanced.return <- function(weight, ia) {
  weight = weight[, 1:ia$n, drop = F]
  total.return = apply(1 + ia$hist.returns, 2, prod)
  total.portfolio.return = weight %*% total.return / rowSums(weight)
  total.portfolio.return = 
    (total.portfolio.return^(1 / nrow(ia$hist.returns)))^
    ia$annual.factor - 1
  return(total.portfolio.return)
  }

geom2aritm <- function(G, V, a, b) {
  (2 * G + a * V^2) / (1 - b * G + sqrt((1 + b * G)^2 + 2 * a * b * V^2))
  }

aritm2geom <- function(R, V, a, b) {
  R - a * V^2 / (2 * (1 + b * R))
  }

geom2aritm4 <- function(G, V) {
  (1 + G) * sqrt(1/2 + 1/2 * sqrt(1 + 4 * V^2/(1 + G)^2)) - 1
  }

aritm2geom4 <- function(R, V) {
  (1 + R)/(sqrt(1 + V^2/(1 + R)^2)) - 1
  }

target.return.portfolio.helper <- function(ia, constraints, target.return) {
  constraints.target = add.constraints(ia$expected.return, type = '>=', 
                                       b = target.return, constraints)
  sol = try(min.var.portfolio(ia, constraints.target), silent = TRUE)
  if(inherits(sol, 'try-error')) sol = max.return.portfolio(ia, constraints)
  sol
  }

target.return.portfolio <- function(target.return, annual.factor = 252) {
  target.return = as.double(target.return[1])
  if(target.return > 1) target.return = target.return / 100
  target.return = target.return / annual.factor
  function(ia, constraints) {
    target.return.portfolio.helper(ia, constraints, target.return)
    }
  }

target.risk.portfolio.helper <- function(ia, constraints, target.risk, 
                                         silent = TRUE, min.w = NA, 
                                         max.w = NA) {
  if(is.na(max.w)) max.w = max.return.portfolio(ia, constraints)
  if(is.na(min.w)) min.w = min.var.portfolio(ia, constraints)
  max.r = portfolio.return(max.w, ia)
  min.r = portfolio.return(min.w, ia)
  max.s = portfolio.risk(max.w, ia)
  min.s = portfolio.risk(min.w, ia)
  
if(target.risk >= min.s & target.risk <= max.s) {
  f <- function(x, ia, constraints, target.risk) {
    portfolio.risk(target.return.portfolio.helper(ia, constraints, x), ia) - 
      target.risk
    }
  
  f.lower = min.s - target.risk
  f.upper = max.s - target.risk
  sol = uniroot(f, c(min.r, max.r), f.lower = f.lower, f.upper = f.upper, 
                tol = 0.0001, ia = ia, constraints = constraints, 
                target.risk = target.risk)
  
  if(!silent) cat('Found solution in', sol$iter, 'itterations', '\n')
  return(target.return.portfolio.helper(ia, constraints, sol$root))
  } else if(target.risk < min.s) {
    return(min.w)
  } else {
    return(max.w)
  }
  
  stop(paste('target.risk =', target.risk, 'is not possible, max risk =', 
             max.s, ', min risk =', min.s))
  }

target.risk.portfolio <- function(target.risk, annual.factor = 252) {
  target.risk = as.double(target.risk[1])
  if(target.risk > 1) target.risk = target.risk / 100
  target.risk = target.risk / sqrt(annual.factor)
  function(ia, constraints) {
    target.risk.portfolio.helper(ia, constraints, target.risk)
    }
  }

min.risk.portfolio <- function(ia, constraints) {
  x = NA
  binary.vec = 0
  if(!is.null(constraints$binary.index)) 
    binary.vec = constraints$binary.index
  
  if(is.null(ia$cov.temp)) ia$cov.temp = ia$cov
  sol = try(solve.QP.bounds(Dmat = ia$cov.temp, 
                            dvec = rep(0, nrow(ia$cov.temp)), 
                            Amat = constraints$A, bvec = constraints$b, 
                            constraints$meq, lb = constraints$lb, 
                            ub = constraints$ub, binary.vec = binary.vec), 
            TRUE)
  
  if(!inherits(sol, 'try-error')) {
    if(binary.vec[1] != 0) cat(sol$counter,'QP calls made to solve problem 
                               with', len(constraints$binary.index), 
                               'binary variables using Branch&Bound', '\n')
    x = sol$solution;
    }
  return(x)
  }

portfolio.risk <- function(weight, ia) {
  if(is.null(dim(weight))) dim(weight) = c(1, len(weight))
  weight = weight[, 1:ia$n, drop = FALSE]
  cov = ia$cov[1:ia$n, 1:ia$n]
  return(apply(weight, 1, function(x) sqrt(t(x) %*% cov %*% x)))
  }

find.erc.portfolio <- function(ia, constraints) {
  cov = ia$cov[1:ia$n, 1:ia$n]
  fn <- function(x){
    risk.contribution = (x * (cov %*% x))
    sum(abs(risk.contribution - mean(risk.contribution)))
    }
  
  x = optimize.portfolio.nlp(ia, constraints, fn)
  return(x)
  }

find.erc.portfolio.simple <- function(ia, constraints) {
  cov = ia$cov[1:ia$n, 1:ia$n]
  fn <- function(x){
    if (sum(x) == 0) x = x + 1e-2
    x  = x / sum(x)
    risk.contribution = (x * (cov %*% x))
    var(as.double(risk.contribution))
    }
  
  x0 = 1/sqrt(diag(cov))
  x0 = x0 / sum(x0)
  x = nlminb(start = x0, objective = fn, lower = constraints$lb, 
             upper = constraints$ub)
  x$par = x$par / sum(x$par)
  return(x$par)
  }

find.erc.portfolio.test <- function() {
  ia = aa.test.create.ia.rebal()
  n = ia$n
  x0 = 1/sqrt(diag(ia$cov))
  temp = x0 / sum(x0)
  rc.temp = portfolio.risk.contribution(temp, ia)
  rc.temp = abs(as.vector(rc.temp))
  plot(rc.temp,ylim=c(0, 0.4))
  diff(range(rc.temp))
  sd(rc.temp)
  
  constraints = new.constraints(n, lb = 0, ub = 1)
  constraints = add.constraints(rep(1, n), 1, type = '=', constraints)
  temp = find.erc.portfolio(ia, constraints)
  rc.temp = portfolio.risk.contribution(temp, ia)
  rc.temp = abs(as.vector(rc.temp))
  plot(rc.temp,ylim=c(0, 0.4))
  diff(range(rc.temp))
  sd(rc.temp)
  
  temp = find.erc.portfolio.simple(ia, constraints)
  temp = temp / sum(temp)
  rc.temp = portfolio.risk.contribution(temp, ia)
  rc.temp = abs(as.vector(rc.temp))
  plot(rc.temp,ylim=c(0, 0.4))
  diff(range(rc.temp))
  sd(rc.temp)
  
  temp = equal.risk.contribution.portfolio(ia, constraints)
  temp = temp / sum(temp)
  rc.temp = portfolio.risk.contribution(temp, ia)
  rc.temp = abs(as.vector(rc.temp))
  plot(rc.temp, ylim = c(0, 0.4))
  diff(range(rc.temp))
  sd(rc.temp)
  }

portfolio.risk.contribution <- function(weight, ia) {
  if(is.null(dim(weight))) dim(weight) = c(1, len(weight))
  weight = weight[, 1:ia$n, drop = FALSE]
  cov = ia$cov[1:ia$n, 1:ia$n]
  out = weight
  out[] = t(apply(weight, 1, function(x) 
    (x * (cov %*% x)) / (t(x) %*% cov %*% x)[1]))
  return(out)
  }

portopt <- function(ia, constraints = NULL, nportfolios = 50, name = 'Risk', 
                    min.risk.fn = min.risk.portfolio, 
                    equally.spaced.risk = FALSE) {
  load.packages('quadprog,corpcor,lpSolve,kernlab')
  if(is.null(constraints)) {
    constraints = new.constraints(rep(0, ia$n), 0, type = '>=')
    }
  ia$risk = iif(ia$risk == 0, 0.000001, ia$risk)
  
  if(is.null(ia$cov)) ia$cov = ia$correlation * (ia$risk %*% t(ia$risk))
  ia$cov.temp = ia$cov
  n0 = ia$n
  n = nrow(constraints$A)
  if(n != nrow(ia$cov.temp)) {
    temp =  matrix(0, n, n)
    temp[1:n0, 1:n0] = ia$cov.temp[1:n0, 1:n0]
    ia$cov.temp = temp
    }
  
  if(!is.positive.definite(ia$cov.temp, method = 'chol')) {
    ia$cov.temp <- make.positive.definite(ia$cov.temp, 0.000000001)
    }
  
  if(nportfolios<2) nportfolios = 2
  out = list(weight = matrix(NA, nportfolios, nrow(constraints$A)))
  colnames(out$weight) = rep('', ncol(out$weight))
  colnames(out$weight)[1:ia$n] = ia$symbols
  out$weight[nportfolios, ] = max.return.portfolio(ia, constraints)
  out$weight[1, ] = match.fun(min.risk.fn)(ia, constraints)
  constraints$x0 = out$weight[1, ]
  
  if(nportfolios > 2) {
    out$return = portfolio.return(out$weight, ia)
    target = seq(out$return[1], out$return[nportfolios], 
                 length.out = nportfolios)
    
    constraints = add.constraints(c(ia$expected.return, 
                                    rep(0, nrow(constraints$A) - ia$n)), 
                                  target[1], type = '>=', constraints)
    
    for(i in 2:(nportfolios - 1)) {
      constraints$b[len(constraints$b)] = target[i]
      out$weight[i, ] = match.fun(min.risk.fn)(ia, constraints)
      constraints$x0 = out$weight[i, ]
      }
    
    if(equally.spaced.risk) {
      out$risk = portfolio.risk(out$weight, ia)
      temp = diff(out$risk)
      index = which(temp >= median(temp) + mad(temp))
      
      if( len(index) > 0 ) {
        index = min(index)
        proper.spacing = ceiling((out$risk[nportfolios] - 
                                    out$risk[index]) / temp[(index - 1)]) - 1
        nportfolios1 = proper.spacing + 2
        
        if(nportfolios1 > 2) {
          out$return = portfolio.return(out$weight, ia)
          out$risk = portfolio.risk(out$weight, ia)
          temp = spline(out$risk, out$return, n = nportfolios, 
                        method = 'natural')
          target = temp$y[which(temp$y > out$return[index] & 
                                  temp$y < out$return[nportfolios] & 
                                  temp$x > out$risk[index] & 
                                  temp$x < out$risk[nportfolios])]
          
          target = c(out$return[index], target, out$return[nportfolios])
          nportfolios1 = len(target)
          out1 = list(weight = matrix(NA, nportfolios1, nrow(constraints$A)))
          out1$weight[1, ] = out$weight[index, ]
          out1$weight[nportfolios1, ] = out$weight[nportfolios, ]
          constraints$x0 = out1$weight[1, ]
          
          for(i in 2:(nportfolios1 - 1)) {
            constraints$b[len(constraints$b)] = target[i]
            out1$weight[i, ] = match.fun(min.risk.fn)(ia, constraints)
            constraints$x0 = out1$weight[i, ]
            }
          
          out$weight = rbind(out$weight[-c(index:nportfolios), ], 
                             out1$weight)
          }
        }
      }
    }
  rm.index = is.na(rowSums(out$weight))
  
  if(any(rm.index)) out$weight = out$weight[!rm.index, ]
  out$return = portfolio.return(out$weight, ia)
  out$risk = portfolio.risk(out$weight, ia)
  out$name = name
  return(out)
  }

plot.ia <- function(ia, layout = NULL) {
  if(is.null(layout)) layout(matrix(1:2, nr = 1))
  temp = cbind(ia$expected.return, ia$risk)
  temp[] = plota.format(100 * temp[], 1, '', '%')
  colnames(temp) = spl('Return,Risk')
  plot.table(temp, 'Symbol')
  temp = ia$correlation
  temp[lower.tri(temp, TRUE)] = NA
  temp = temp[-ia$n, -1]
  temp[] = plota.format(100 * temp[], 1, '', '%')
  plot.table(temp, highlight = TRUE, colorbar = TRUE)
  }

plot.ef <- function(ia, efs, portfolio.risk.fn = portfolio.risk, 
                    transition.map = TRUE, layout = NULL) {
  risk.label = as.character(substitute(portfolio.risk.fn))
  n = ia$n
  x = match.fun(portfolio.risk.fn)(diag(n), ia)
  y = ia$expected.return
  xlim = range(c(0, x, max(sapply(efs, function(x) 
    max(match.fun(portfolio.risk.fn)(x$weight,ia))))), na.rm = TRUE)
  ylim = range(c(0, y, min(sapply(efs, function(x) 
    min(portfolio.return(x$weight, ia)))), 
    max(sapply(efs, function(x) max(portfolio.return(x$weight, ia))))), 
    na.rm = TRUE)
  
  x = 100 * x
  y = 100 * y
  xlim = 100 * xlim
  ylim = 100 * ylim
  
  if(!transition.map) layout = TRUE
  if(is.null(layout)) layout(1:2)
  par(mar = c(4, 3, 2 , 1), cex = 0.8)
  plot(x, y, xlim = xlim, ylim = ylim, xlab='', ylab='', 
       main=paste(risk.label, 'vs Return'), col='black')
  mtext('Return', side = 2, line = 2, cex = par('cex'))
  mtext(risk.label, side = 1, line = 2, cex = par('cex'))
  grid();
  text(x, y, ia$symbols,	col = 'blue', adj = c(1, 1), cex = 0.8)
  
  for(i in len(efs):1) {
    ef = efs[[i]]
    x = 100 * match.fun(portfolio.risk.fn)(ef$weight, ia)
    y = 100 * ef$return
    lines(x, y, col=i)
    }
  
  plota.legend(sapply(efs, function(x) x$name), 1:len(efs))
  if(transition.map) {
    plot.transition.map(efs[[i]]$weight, x, risk.label, efs[[i]]$name)
    }
  }

plot.add.portfolios = function(ia, portfolio.risk.fn = portfolio.risk, ...) {
  portfolios = lst(...)
  col = plota.colors(portfolios)
  
  for(i in 1:len(portfolios)) 
    points(100 * portfolio.risk.fn(portfolios[[i]], ia), 100 * 
             portfolio.return(portfolios[[i]], ia), pch = 15, col = col[i])
  plota.legend(names(portfolios), col, x = 'bottomright')
  }

plot.transitopn.map <- function(x, y, xlab = 'Risk', name = '', 
                                type = c('s','l')) {
  plot.transition.map(x, y , xlab, name, type)
  }

plot.transition.map <- function(y, x, xlab = 'Risk', name = '', 
                                type = c('s', 'l'), col = NA) {
  if(is.list(y)) {
    name = y$name
    x = 100 * y$risk
    y = y$weight
    }
  
  y[is.na(y)] = 0
  par(mar = c(4, 3, 2, 1), cex = 0.8)
  plota.stacked(x, y, xlab = xlab, main = paste('Transition Map for', name), 
                type = type[1], col = ifna(col, plota.colors(ncol(y))))
  }

portfolio.turnover <- function(weight) {
  if(is.null(dim(weight))) dim(weight) = c(1, len(weight))
  out = weight[, 1] * NA
  out[] = rowSums(abs(weight - mlag(weight))) / 2
  return(out)
  }

portfolio.concentration.herfindahl.index <- function(weight) {
  if(is.null(dim(weight))) dim(weight) = c(1, len(weight))
  one.over.n = 1 / rowSums(!is.na(weight))
  out = weight[, 1] * NA
  out[] = (rowSums(weight^2, na.rm = TRUE) - one.over.n) / (1 - one.over.n)
  return(out)
  
  one.over.n = 1 / ncol(weight)
  out = weight[, 1] * NA
  out[] = (rowSums(weight^2) - one.over.n) / (1 - one.over.n)
  return(out)
  }

portfolio.concentration.gini.coefficient <- function(weight) {
  if(is.null(dim(weight))) dim(weight) = c(1, len(weight))
  n = ncol(weight)
  one.to.n = 1:n
  out = weight[, 1] * NA
  
  for(i in 1:nrow(weight)) {
    x = coredata(weight[i, ])
    index = !is.na(x)
    n1 = sum(index)
    if(n1 > 0) {
      temp = sort(x[index], decreasing = TRUE)
      if(temp[n1] < 0) temp = temp - temp[n1]
      out[i] = (n1+1) / (n1-1) - 2 * sum(temp * one.to.n[1:n1]) / 
        (n1 * (n1 - 1) * sum(temp) / n1)
      }
    }
  return(out)
  
  out[] = apply(weight, 1, function(x) {
    temp = sort(x, decreasing = TRUE)
    sum(temp * one.to.n)
    })
  
  out = (n + 1) / (n - 1) - 2 * out /(n * (n - 1) * apply(weight, 1, mean))
  return(out)
  }

binary_branch_bound <- function(index_binvar, bbb_data, bbb_solve, 
                                control = bbb_control()) {
  fbest = Inf
  xbest = 0 * bbb_data$x0
  counter = 0
  nbinvar = length(index_binvar)
  flag = 7
  
  stack = new.env()
  stack$data = list()
  stack$cost = c()
  stack$pointer = c()
  stack$data[[1]] = list(lb = bbb_data$lb, ub = bbb_data$ub, var = 1:nbinvar, 
                         path = rep(0, nbinvar), level = 0, fval = Inf)
  stack$cost = 0
  stack$pointer = 1
  control$proborder.selected = control$proborder
  
  if(F) {
    lb = bbb_data$lb
    ub = bbb_data$ub
    for(i in 0:1) {
      lb[] = i
      ub[] = i
      sol = match.fun(bbb_solve)(bbb_data, lb, ub)
      
      if(sol$ok) {
        x = sol$x
        fval = sol$fval
        xi = x[index_binvar]
        
        if (max(abs(round(xi, 0) - xi)) < control$bineps) {
          fbest = fval
          xbest = x
          flag = 1
          
          if(!control$silent) cat('FOUND SOLUTION =', fbest, '\n');
          }
        }
      }
    }
  
  while(length(stack$data) > 0) {
    subprob = bbb_pop(stack)
    
    if(!control$silent) {
      cat('-----------------------------------------------------', '\n')
      
      if(max(subprob$path) > 0) {
        temp.index = order(-subprob$path)[1:sum(subprob$path > 0)]
        cat('\t', paste('b', temp.index, ' = ', subprob$lb[temp.index], 
                        sep = ''), '\n')
      } else {
        cat(counter, '\t', 'FIRST NODE', '\n')
      }
      
      cat(counter, '\t', subprob$lb, '\t', subprob$var, '\t', subprob$fval, 
          '\t', fbest, '\n')
      cat('\t', subprob$ub, '\n')
      cat('stack size =', len(stack$pointer), '\n')
      }
    
    if(is.finite(subprob$fval) & is.finite(fbest) & fbest <= subprob$fval) {
      
      if(!control$silent) 
        cat('SKIP this problem because a solution with lower FVAL already found\n')
      } else {
        counter = counter + 1
        sol = match.fun(bbb_solve)(bbb_data, subprob$lb, subprob$ub)
        
        if(!sol$ok) {
          
          if(!control$silent) cat('NO SOLUTION EXISTS\n\n');
          } else {
            x = sol$x
            fval = sol$fval
            
            if(!control$silent) {
              cat('SOLUTION OK', '\t', sol$fval, '\n')
              cat('\t', round(x[index_binvar[subprob$var]], 3), '\n\n')
              }
            
            if(flag != 1) flag = 5
            
            if(fval <= fbest) {
              if(length(subprob$var) == 0) {
                fbest = fval
                xbest = x
                flag = 1
                
                if(!control$silent) cat('FOUND SOLUTION =', fbest, '\n');
                } else {
                  xi = x[index_binvar[subprob$var]]
                  
                  if(max(abs(round(xi, 0) - xi)) < control$bineps) {
                    fbest = fval
                    xbest = x
                    flag = 1
                    
                    if(!control$silent) cat('FOUND SOLUTION =', fbest, '\n');
                    } else {
                      branchvar = bbb_decision(xi, control)
                      probs = bbb_separate(subprob, branchvar, fval)
                      p0 = probs$p0
                      p1 = probs$p1
                      if(!control$silent) cat('Branch on =', 
                                              subprob$var[branchvar], '\n');
                      if(control$searchdir == 0) {
                        cost = 1 / (subprob$level + 1)
                        } else if(control$searchdir == 1) {
                          cost = subprob$level + 1
                          } else if(control$searchdir == 2) {
                            cost = fval
                            } else if(control$searchdir == 3) {
                              cost = fval/(subprob$level + 1)
                              }
                      if(control$proborder == 2) {
                        control$proborder.selected = round(xi[branchvar], 0)
                        }
                      if(control$proborder.selected == 0) {
                        bbb_push(stack, p1, p0, cost)
                        } else {
                          bbb_push(stack, p0, p1, cost)
                        }
                    }
                }
            }
          }
        if(F) {
          cat('counter =', counter, '\n')
          cat('fbest     =', fbest, '\n')
          cat('stack$pointer =', stack$pointer, '\n')
          cat('\n')
        }
      }
    }
  rm(list=ls(stack,all=TRUE), envir=stack)
  
  return(list(xmin = xbest, fmin = fbest, counter = counter, flag = flag))
  }

bbb_decision <- function(xi, control) {
  if(control$branchvar == 0) {
    branchvar = 1
    } else if(control$branchvar == 1) {
      branchvar = which.max(abs(xi - round(xi, 0)))
      } else if(control$branchvar == 2) {
        branchvar = which.min(abs(xi - round(xi, 0)))
        } else {
          branchvar = 1
        }
  return(branchvar)
  }

bbb_pop <- function(stack) {
  i = stack$pointer[length(stack$data)]
  subprob = stack$data[[i]]
  stack$pointer[stack$pointer > i] = stack$pointer[stack$pointer > i] - 1
  stack$data[[i]] = NULL
  length(stack$cost) = length(stack$data)
  length(stack$pointer) = length(stack$data)
  return(subprob)
  }

bbb_push <- function(stack, element1, element2, cost) {
  n = length(stack$data)
  i = match(TRUE, stack$cost <= cost)
  
  if(is.na(i)) i = n else i = i - 1
  stack$data[[(n + 1)]] = element1
  stack$data[[(n + 2)]] = element2
  
  if(i == 0) {
    stack$pointer = c((n + 1), (n + 2), stack$pointer)
    stack$cost = c(cost, cost, stack$cost)
    } else {
      stack$pointer = c(stack$pointer[1:i], (n + 1), (n + 2), 
                        stack$pointer[-c(1:i)])
      stack$cost = c(stack$cost[1:i], cost, cost, stack$cost[-c(1:i)])
      }
    }

bbb_separate <- function(prob, branchvar, fval) {
  if(length(prob$var) >= 1) {
    p0 = prob
    p0$fval = fval
    p0$level = prob$level + 1
    p0$var = prob$var[-branchvar]
    p0$path[prob$var[branchvar]] = 1 + max(p0$path)
    p1 = p0
    p0$lb[prob$var[branchvar]] = 0
    p0$ub[prob$var[branchvar]] = 0
    p1$lb[prob$var[branchvar]] = 1
    p1$ub[prob$var[branchvar]] = 1
    } else {
      stop('no more integer variables to branch on')
      }
  return(list(p0 = p0, p1 = p1))
  }

bt.merge <- function(b, align = c('keep.all', 'remove.na'), dates = NULL) {
  align = align[1]
  symbolnames = b$symbolnames
  nsymbols = len(symbolnames)
  ncount = sapply(symbolnames, function(i) nrow(b[[i]]))
  all.dates = double(sum(ncount))
  itemp = 1
  
  for(i in 1:nsymbols) {
    all.dates[itemp:(itemp + ncount[i] -1)] = attr(b[[symbolnames[i]]], 
                                                   'index')
    itemp = itemp + ncount[i]
    }
  
  temp = sort(all.dates)
  unique.dates = c(temp[1], temp[-1][diff(temp) != 0])
  if(!is.null(dates)) {
    class(unique.dates) = c('POSIXct', 'POSIXt')
    temp = make.xts(integer(len(unique.dates)), unique.dates)
    unique.dates = attr(temp[dates], 'index')
    }
  
  date.map = matrix(NA, nr = len(unique.dates), nsymbols)
  itemp = 1
  for(i in 1:nsymbols) {
    index = match(all.dates[itemp:(itemp + ncount[i] - 1)], unique.dates)
    sub.index = which(!is.na(index))
    date.map[index[sub.index], i] = sub.index
    itemp = itemp + ncount[i]
    }
  
  index = c()
  if(align == 'remove.na') {
    index = which(count(date.map, side = 1) < nsymbols)
    }
  
  if(len(index) > 0) {
    date.map = date.map[-index, , drop = FALSE]
    unique.dates = unique.dates[-index]
    }
  
  class(unique.dates) = c('POSIXct', 'POSIXt')
  return( list(all.dates = unique.dates, date.map = date.map))
  }

bt.prep <- function(b, align = c('keep.all', 'remove.na'), dates = NULL, 
                    fill.gaps = FALSE, basic = FALSE) {
  if(!exists('symbolnames', b, inherits = FALSE)) b$symbolnames = ls(b)
  symbolnames = b$symbolnames
  nsymbols = len(symbolnames)
  
  if(nsymbols > 1) {
    out = bt.merge(b, align, dates)
    for(i in 1:nsymbols) {
      temp = coredata(b[[symbolnames[i]]])[out$date.map[, i], , drop = FALSE]
      b[[symbolnames[i]]] = iif(basic, temp, make.xts(temp, out$all.dates))
      map.col = find.names('Close,Volume,Open,High,Low,Adjusted', 
                           b[[symbolnames[i]]])
      
      if(fill.gaps & !is.na(map.col$Close)) {
        close = coredata(b[[symbolnames[i]]][, map.col$Close])
        n = len(close)
        last.n = max(which(!is.na(close)))
        close = ifna.prev(close)
        
        if(last.n + 5 < n) close[last.n:n] = NA
        b[[symbolnames[i]]][, map.col$Close] = close
        index = !is.na(close)
        
        if(!is.na(map.col$Volume)) {
          index1 = is.na(b[[symbolnames[i]]][, map.col$Volume]) & index
          b[[symbolnames[i]]][index1, map.col$Volume] = 0
          }
        
        for(field in spl('Open,High,Low,Adjusted')) {
          j = map.col[[field]]
          if(!is.null(j)) {
            index1 = is.na(b[[symbolnames[i]]][,j]) & index
            b[[symbolnames[i]]][index1, j] = close[index1]
            }
          }
        }
      }
    } else {
      if(!is.null(dates)) b[[symbolnames[1]]] = b[[symbolnames[1]]][dates, ]
      out = list(all.dates = index.xts(b[[symbolnames[1]]]))
      
      if(basic) b[[symbolnames[1]]] = coredata(b[[symbolnames[1]]])
      }
  b$dates = out$all.dates
  dummy.mat = matrix(double(), len(out$all.dates), nsymbols)
  colnames(dummy.mat) = symbolnames
  
  if(!basic) dummy.mat = make.xts(dummy.mat, out$all.dates)
  b$weight = dummy.mat
  b$execution.price = dummy.mat
  for(i in 1:nsymbols) {
    if(has.Cl(b[[symbolnames[i]]])) {
      dummy.mat[, i] = Cl(b[[symbolnames[i]]]);
      }
    }
  b$prices = dummy.mat
  }

bt.prep.matrix <- function(b, align = c('keep.all', 'remove.na'), 
                           dates = NULL, basic = FALSE) {
  align = align[1]
  nsymbols = len(b$symbolnames)
  if(!is.null(dates)) {
    temp = make.xts(1:len(b$dates), b$dates)
    temp = temp[dates]
    index = as.vector(temp)
    for(i in b$fields) b[[i]] = b[[i]][index, , drop = FALSE]
    b$dates = b$dates[index]
    }
  
  if(align == 'remove.na') {
    index = which(count(b$Cl, side=1) < nsymbols)
    } else {
      index = which(count(b$Cl, side = 1) < max(1, 0.1 * nsymbols))
      }
  
  if(len(index) > 0) {
    for(i in b$fields) b[[i]] = b[[i]][-index, , drop = FALSE]
    b$dates = b$dates[-index]
    }
  
  dummy.mat = iif(basic, b$Cl, make.xts(b$Cl, b$dates))
  b$weight = NA * dummy.mat
  b$execution.price = NA * dummy.mat
  b$prices = dummy.mat
  }

bt.prep.matrix.test <- function() {
  load.packages('quantmod')
  returns = read.xts('Example.csv', date.fn=function(x) paste('1',x), 
                     format='%d %b-%y')
  prices = bt.apply.matrix(1 + returns, cumprod)
  data <- new.env()
  data$symbolnames = colnames(prices)
  data$dates = index(prices)
  data$fields = 'Cl'
  data$Cl = prices
  bt.prep.matrix(data)
  data$weight[] = NA
  data$weight[] = 1
  buy.hold = bt.run.share(data)
  plotbt(buy.hold, plotX = TRUE, log = 'y', LeftMargin = 3)
  mtext('Cumulative Performance', side = 2, line = 1)
  }

bt.prep.remove.symbols.min.history <- function(b, min.history = 1000) {
  bt.prep.remove.symbols(b, which(count(b$prices, side=2) < min.history))
  }

bt.prep.remove.symbols <- function(b, index) {
  if(len(index) > 0) {
    if(is.character(index)) index = match(index, b$symbolnames)
    b$prices = b$prices[, -index]
    b$weight = b$weight[, -index]
    b$execution.price = b$execution.price[, -index]
    env.rm(b$symbolnames[index], b)
    b$symbolnames = b$symbolnames[-index]
    }
  }

bt.prep.trim <- function(b, dates = NULL) {
  if(is.null(dates)) return(b)
  dates.index = dates2index(b$prices, dates)
  data.copy <- new.env()
  
  for(s in b$symbolnames) data.copy[[s]] = b[[s]][dates.index, , drop = F]
  data.copy$symbolnames = b$symbolnames
  data.copy$dates = b$dates[dates.index]
  data.copy$prices = b$prices[dates.index, , drop = F]
  data.copy$weight = b$weight[dates.index, , drop = F]
  data.copy$execution.price = b$execution.price[dates.index, , drop = F]
  return(data.copy)
  }

bt.run.share <- function(b, prices = b$prices, clean.signal = T, 
                         trade.summary = F, do.lag = 1, 
                         do.CarryLastObservationForwardIfNA = T, silent = F, 
                         capital = 100000, commission = 0, weight = b$weight, 
                         dates = 1:nrow(b$prices)) {
  prices[] = bt.apply.matrix(coredata(prices), ifna.prev)
  weight = mlag(weight, do.lag - 1)
  do.lag = 1
  
  if(clean.signal) weight[] = bt.exrem(weight)
  weight = (capital / prices) * weight
  bt.run(b, trade.summary = trade.summary, do.lag = do.lag, 
         do.CarryLastObservationForwardIfNA = 
           do.CarryLastObservationForwardIfNA, type = 'share', 
         silent = silent, capital = capital, commission = commission, 
         weight = weight, dates = dates)
  }

bt.run <- function(b, trade.summary = F, do.lag = 1, 
                   do.CarryLastObservationForwardIfNA = T, 
                   type = c('weight', 'share'), silent = F, capital = 100000, 
                   commission = 0, weight = b$weight, 
                   dates = 1:nrow(b$prices)) {
  dates.index = dates2index(b$prices, dates)
  type = type[1]
  weight[] = ifna(weight, NA)
  
  if(do.lag > 0) weight = mlag(weight, do.lag)
  if(do.CarryLastObservationForwardIfNA)
    weight[] = apply(coredata(weight), 2, ifna.prev)
  weight[is.na(weight)] = 0
  weight1 = mlag(weight, -1)
  tstart = weight != weight1 & weight1 != 0
  tend = weight != 0 & weight != weight1
  trade = ifna(tstart | tend, FALSE)
  prices = b$prices
  
  if(sum(trade) > 0) {
    execution.price = coredata(b$execution.price)
    prices1 = coredata(b$prices)
    prices1[trade] = iif(is.na(execution.price[trade]), prices1[trade], 
                        execution.price[trade])
    prices[] = prices1
    }
  
  if(type == 'weight') {
    ret = prices / mlag(prices) - 1
    ret[] = ifna(ret, NA)
    ret[is.na(ret)] = 0
    } else {
      ret = prices
      }
  
  temp = b$weight
  temp[] = weight
  weight = temp
  bt = bt.summary(weight, ret, type, b$prices, capital, commission)
  bt$dates.index = dates.index
  bt = bt.run.trim.helper(bt, dates.index)
  
  if(trade.summary) bt$trade.summary = bt.trade.summary(b, bt)
  if(!silent) {
    cat('Latest weights :\n')
    print(round(100 * last(bt$weight), 2))
    cat('\n')
    cat('Performance summary :\n')
    cat('', spl('CAGR,Best,Worst'), '\n', sep = '\t')
    cat('', sapply(cbind(bt$cagr, bt$best, bt$worst), function(x) 
      round(100 * x, 1)), '\n', sep = '\t')
    cat('\n')
    }
  return(bt)
  }

bt.run.trim.helper = function(bt, dates.index) {
  n.dates = len(dates.index)
  for(n in ls(bt)) {
    if(!is.null(dim(bt[[n]]))) {
      if(nrow(bt[[n]]) > n.dates)
        bt[[n]] = bt[[n]][dates.index, , drop = F]
      } else if(len(bt[[n]]) > n.dates)
        bt[[n]] = bt[[n]][dates.index]
      }
  
  bt$equity = bt$equity / as.double(bt$equity[1])
  bt$best = max(bt$ret)
  bt$worst = min(bt$ret)
  bt$cagr = compute.cagr(bt$equity)
  bt
  }

bt.summary <- function(weight, ret, type = c('weight', 'share'), 
                       close.prices, capital = 100000, commission = 0) {
  if(!is.list(commission)) {
    if(type == 'weight')
      commission = list(cps = 0.0, fixed = 0.0, percentage = commission)
    else
      commission = list(cps = commission, fixed = 0.0, percentage = 0.0)
    }
  
  type = type[1]
  n = nrow(ret)
  bt = list()
  bt$weight = weight
  bt$type = type
  com.weight = mlag(weight, -1)
  if(type == 'weight') {
    temp = ret[, 1]
    temp[] = rowSums(ret * weight) - 
      rowSums(abs(com.weight - mlag(com.weight)) * commission$percentage, 
              na.rm = T) - rowSums(sign(abs(com.weight - mlag(com.weight))
                                        ) * commission$fixed, na.rm = T)
    bt$ret = temp
    } else {
      bt$share = weight
      bt$capital = capital
      prices = ret
      prices[] = bt.apply.matrix(coredata(prices), ifna.prev)
      close.prices[] = bt.apply.matrix(coredata(close.prices), ifna.prev)
      cash = capital - rowSums(bt$share * mlag(close.prices), na.rm = T)
      share.nextday = mlag(bt$share, -1)
      tstart = bt$share != share.nextday & share.nextday != 0
      tend = bt$share != 0 & bt$share != share.nextday
      trade = ifna(tstart | tend, FALSE)
      tstart = trade
      
      index = mlag(apply(tstart, 1, any))
      index = ifna(index, FALSE)
      index[1] = T
      totalcash = NA * cash
      totalcash[index] = cash[index]
      totalcash = ifna.prev(totalcash)
      totalcash = ifna(totalcash, 0)
      portfolio.ret = (totalcash + rowSums(bt$share * prices, na.rm = T) - 
                         rowSums(abs(com.weight - mlag(com.weight)) * 
                                   commission$cps, na.rm = T) - 
                         rowSums(sign(abs(com.weight - mlag(com.weight))) * 
                                   commission$fixed, na.rm = T) - 
                         rowSums(prices * abs(com.weight - mlag(com.weight)
                                              ) * commission$percentage, 
                                 na.rm = T)) / (totalcash + 
                                                  rowSums(bt$share * 
                                                            mlag(prices), 
                                                          na.rm = T) ) - 1
      bt$weight = bt$share * mlag(prices) / 
        (totalcash + rowSums(bt$share * mlag(prices), na.rm = T))
      bt$weight[is.na(bt$weight)] = 0
      temp = ret[, 1]
      temp[] = ifna(portfolio.ret, 0)
      temp[1] = 0
      bt$ret = temp
      }
  
  bt$best = max(bt$ret)
  bt$worst = min(bt$ret)
  bankrupt = which(bt$ret <= -1)
  
  if(len(bankrupt) > 0) bt$ret[bankrupt[1]:n] = -1
  bt$equity = cumprod(1 + bt$ret)
  bt$cagr = compute.cagr(bt$equity)
  return(bt)
  }

bt.summary.test <- function() {
  load.packages('quantmod')
  data <- new.env()
  getSymbols('EEM', src = 'yahoo', from = '1980-01-01', env = data, 
             auto.assign = T)
  bt.prep(data, align='keep.all', dates='2013:08::2013:08:10')
  buy.date = '2013:08:05'
  sell.date = '2013:08:06'
  coredata(data$prices) <-c(10, 20, 40, 60, 20, 160, 60)
  prices = data$prices
  data$weight[] = NA
  data$weight[buy.date] = 1
  data$weight[sell.date] = 0
  commission = list(cps = 0.0, fixed = 0.0, percentage = 1/100)
  model3 = bt.run(data, commission = commission, silent = T)
  model3$ret
  data$weight[] = NA
  data$weight[buy.date] = 1
  data$weight[sell.date] = 0
  commission = list(cps = 0.0, fixed = 0.0, percentage = 1/100)
  model3 = bt.run.share(data, commission = commission, capital = 100000, 
                        silent = T)
  model3$ret
  }

bt.trim <- function(..., dates = '::') {
  models = variable.number.arguments( ... )
  for(i in 1:len(models)) {
    bt = models[[i]]
    n = len(bt$equity)
    first = which.max(!is.na(bt$equity) & bt$equity != 1)
    
    if(first > 1 && !is.na(bt$equity[(first-1)]))
      first = first - 1
    if(first < n) {
      index = first:n
      dates.range = range(dates2index(bt$equity[index], dates))
      index = index[dates.range[1]] : index[dates.range[2]]
      bt$dates.index = bt$dates.index[index]
      bt$equity = bt$equity[index]
      bt$equity = bt$equity / as.double(bt$equity[1])
      bt$ret = bt$ret[index]
      bt$weight = bt$weight[index, , drop = F]
      if (!is.null(bt$share)) bt$share = bt$share[index, , drop = F]
      bt$best = max(bt$ret)
      bt$worst = min(bt$ret)
      bt$cagr = compute.cagr(bt$equity)
      }
    models[[i]] = bt
    }
  return (models)
  }

bt.trim.test <- function() {
  load.packages('quantmod')
  data <- new.env()
  getSymbols(spl('SPY,GLD'), src = 'yahoo', from = '1980-01-01', env = data, 
             auto.assign = T)
  
  for(i in ls(data)) data[[i]] = adjustOHLC(data[[i]], use.Adjusted = T)
  bt.prep(data, align='keep.all')
  models = list()
  data$weight[] = NA
  data$weight$SPY[] = 1
  models$SPY = bt.run.share(data, clean.signal = F)
  data$weight[] = NA
  data$weight$GLD[] = 1
  models$GLD = bt.run.share(data, clean.signal = F)
  strategy.performance.snapshoot(bt.trim(models), T)
  }

bt.run.weight.fast <- function(b, do.lag = 1, 
                               do.CarryLastObservationForwardIfNA = TRUE) {
  weight = ifna(coredata(b$weight), NA)
  if(do.lag > 0) weight = mlag(weight, do.lag)
  if(do.CarryLastObservationForwardIfNA) weight[] = apply(coredata(weight), 
                                                          2, ifna.prev)
  weight[is.na(weight)] = 0
  prices = coredata(b$prices)
  ret = prices / mlag(prices) - 1
  ret[] = ifna(ret, 0)
  ret = rowSums(ret * weight)
  list(weight = weight, ret = ret, equity = cumprod(1 + ret))
  }

compute.turnover <- function(bt, b) {
  year.ends = unique(c(endpoints(bt$weight, 'years'), nrow(bt$weight)))
  year.ends = year.ends[year.ends > 0]
  nr = len(year.ends)
  period.index = c(1, year.ends)
  
  if(bt$type == 'weight') {
    portfolio.value = rowSums(abs(bt$weight), na.rm = T)
    portfolio.turnover = rowSums(abs(bt$weight - mlag(bt$weight)), na.rm = T)
    portfolio.turnover[rowSums(!is.na(bt$weight) & !is.na(mlag(bt$weight))
                               ) == 0] = NA
    } else {
      prices = mlag(b$prices[bt$dates.index, , drop = F])
      
      if(is.null(bt$cash)) {
        cash = bt$capital - rowSums(bt$share * prices, na.rm = T)
        share.nextday = mlag(bt$share, -1)
        tstart = bt$share != share.nextday & share.nextday != 0
        index = mlag(apply(tstart, 1, any))
        index = ifna(index, FALSE)
        totalcash = NA * cash
        totalcash[index] = cash[index]
        totalcash = ifna.prev(totalcash)
        } else
          totalcash = bt$cash
        portfolio.value = totalcash + rowSums(bt$share * prices, na.rm = T)
        portfolio.turnover = rowSums(prices * abs(bt$share - mlag(bt$share)), 
                                     na.rm = T)
        portfolio.turnover[rowSums(!is.na(bt$share) & 
                                     !is.na(mlag(bt$share)) & 
                                     !is.na(prices)) == 0] = NA
        }
  
  portfolio.turnover[1:2] = 0
  temp = NA * period.index
  
  for(iyear in 2:len(period.index)) {
    temp[iyear] = sum(portfolio.turnover[
      period.index[(iyear - 1)] : period.index[iyear]], na.rm = T) /
      mean(portfolio.value[period.index[(iyear - 1)] : period.index[iyear]], 
           na.rm = T)
    }
  return(ifna(mean(temp, na.rm=T),0))
  }

compute.max.deviation <- function(bt, target.allocation) {
  weight = bt$weight[-1, ]
  max(abs(weight - repmat(target.allocation, nrow(weight), 1)))
  }

bt.trade.summary <- function(b, bt) {
  
  if(bt$type == 'weight') weight = bt$weight else weight = bt$share
  out = NULL
  weight1 = mlag(weight, -1)
  tstart = weight != weight1 & weight1 != 0
  tend = weight != 0 & weight != weight1
  tstart[1, weight[1,] != 0] = T
  n = nrow(weight)
  tend[n, weight[n,] != 0] = T
  trade = ifna(tstart | tend, FALSE)
  prices = b$prices[bt$dates.index, , drop = F]
  
  if(sum(trade) > 0) {
    execution.price = coredata(b$execution.price[bt$dates.index, , drop = F])
    prices1 = coredata(b$prices[bt$dates.index, , drop = F])
    prices1[trade] = iif(is.na(execution.price[trade]), prices1[trade], 
                         execution.price[trade])
    prices1[is.na(prices1)] = ifna(mlag(prices1), NA)[is.na(prices1)]
    prices[] = prices1
    weight = bt$weight
    
    symbolnames = b$symbolnames
    nsymbols = len(symbolnames)
    ntrades = max(sum(tstart, na.rm = T), sum(tend, na.rm = T))
    trades = matrix(NA, nr = ntrades, nc = 7)
    colnames(trades) = 
      spl('date,symbol,weight,entry.date,exit.date,entry.price,exit.price')
    itrade = 1
    
    for(i in 1:nsymbols) {
      tstarti = which(tstart[, i])
      tendi = which(tend[, i])
      
      if(len(tstarti) > 0) {
        if(len(tendi) > len(tstarti)) tstarti = c(1, tstarti)
        ntrade = len(tstarti)
        ntrade.index = itrade:(itrade + ntrade - 1)
        trades[ntrade.index, ] = cbind((tstarti + 1), i, 
                                      coredata(weight[(tstarti + 1), i]), 
                                      tstarti, tendi, coredata(
                                        prices[tstarti, i]), 
                                      coredata(prices[tendi, i]))
        itrade = itrade + ntrade
        }
      }
    
    out = list()
    out$stats = cbind(
      bt.trade.summary.helper(trades), bt.trade.summary.helper(
        trades[trades[, 'weight'] >= 0, ]), bt.trade.summary.helper(
          trades[trades[, 'weight'] <0, ]))
    
    colnames(out$stats) = spl('All,Long,Short')
    dates = index(weight)
    dates0 = format(dates, '%Y-%m-%d')
    index = order(dates[trades[, 'entry.date']])
    temp = matrix('', nr = nrow(trades), nc = 8)
    colnames(temp) = 
      spl('symbol,weight,entry.date,exit.date,nhold,entry.price,exit.price,return')
    temp[, 'symbol'] = symbolnames[trades[index, 'symbol']]
    temp[, 'weight'] = round(100 * trades[index, 'weight'], 1)
    temp[, 'entry.date'] = dates0[trades[index, 'entry.date']]
    temp[, 'exit.date'] = dates0[trades[index, 'exit.date']]
    temp[, 'nhold'] = as.numeric(dates[trades[index, 'exit.date']] - 
                                   dates[trades[index, 'entry.date']])
    temp[, 'entry.price'] = round(trades[index, 'entry.price'], 2)
    temp[, 'exit.price'] = round(trades[index, 'exit.price'], 2)
    temp[, 'return'] = round(100 * trades[index, 'weight'] * 
                               (trades[index, 'exit.price'] / 
                                  trades[index,'entry.price'] - 1), 2)
    out$trades = temp
    }
  return(out)
  }

bt.trade.summary.old <- function(b, bt) {
  if(bt$type == 'weight') weight = bt$weight else weight = bt$share
  out = NULL
  weight1 = mlag(weight, -1)
  tstart = weight != weight1 & weight1 != 0
  tend = weight != 0 & weight != weight1
  tstart[1, weight[1, ] != 0] = T
  n = nrow(weight)
  tend[n, weight[n, ] != 0] = T
  trade = ifna(tstart | tend, FALSE)
  prices = b$prices[bt$dates.index, , drop = F]
  
  if(sum(trade) > 0) {
    execution.price = coredata(b$execution.price[bt$dates.index, , drop = F])
    prices1 = coredata(b$prices[bt$dates.index, , drop = F])
    prices1[trade] = iif(is.na(execution.price[trade]), prices1[trade], 
                         execution.price[trade])
    prices1[is.na(prices1)] = ifna(mlag(prices1), NA)[is.na(prices1)]
    prices[] = prices1
    weight = bt$weight
    symbolnames = b$symbolnames
    nsymbols = len(symbolnames)
    trades = c()
    
    for(i in 1:nsymbols) {
      tstarti = which(tstart[, i])
      tendi = which(tend[, i])
      if(len(tstarti) > 0) {
        if(len(tendi) > len(tstarti)) tstarti = c(1, tstarti)
        trades = rbind(trades, cbind(i, weight[(tstarti + 1), i], tstarti, 
                                     tendi, tendi - tstarti, 
                                     coredata(prices[tstarti, i]), 
                                     coredata(prices[tendi, i])))
        }
      }
    
    colnames(trades) = 
      spl('symbol,weight,entry.date,exit.date,nhold,entry.price,exit.price')
    out = list()
    out$stats = cbind(
      bt.trade.summary.helper(trades), 
      bt.trade.summary.helper(trades[trades[, 'weight'] >= 0, ]), 
      bt.trade.summary.helper(trades[trades[, 'weight'] <  0, ]))
    
    colnames(out$stats) = spl('All,Long,Short')
    temp.x = index.xts(weight)
    trades = data.frame(coredata(trades))
    trades$symbol = symbolnames[trades$symbol]
    trades$nhold = as.numeric(temp.x[trades$exit.date] - 
                                temp.x[trades$entry.date])
    trades$entry.date = temp.x[trades$entry.date]
    trades$exit.date = temp.x[trades$exit.date]
    trades$return = round(100 * (trades$weight) * 
                            (trades$exit.price / trades$entry.price - 1), 2)
    trades$entry.price = round(trades$entry.price, 2)
    trades$exit.price = round(trades$exit.price, 2)
    trades$weight = round(100 * (trades$weight), 1)
    out$trades = as.matrix(trades)
    }
  return(out)
  }

bt.trade.summary.test <- function() {
  test = list(weight1 = matrix(c(0, 0, 0, 1), nc = 1), 
              weight2 = matrix(c(0, 1, 0, 0), nc = 1), 
              weight3 = matrix(c(1, 1, 1, 1), nc = 1), 
              weight4 = matrix(c(1, 0, 0, 0), nc = 1), 
              weight5 = matrix(c(1, 2, 0, 1, 2), nc = 1))
  
  for(i in 1:len(test)) {
    weight = test[[i]]
    weight1 = mlag(weight, -1)
    tstart = weight != weight1 & weight1 != 0
    tend = weight != 0 & weight != weight1
    n = nrow(weight)
    tend[n, weight[n, ] != 0] = T
    tend[1, ] = NA
    trade = ifna(tstart | tend, FALSE)
    tstarti = which(tstart)
    tendi = which(tend)
    
    if(len(tendi) > len(tstarti)) tstarti = c(1, tstarti)
    cat(len(tstarti), len(tendi), '\n')
    }
  }

bt.trade.summary.helper <- function(trades) {
  if(nrow(trades) <= 0) return(NA)
  out = list()
  tpnl = trades[, 'weight'] * (trades[, 'exit.price'] / 
                                 trades[, 'entry.price'] - 1)
  tlen = trades[, 'exit.date'] - trades[, 'entry.date']
  
  out$ntrades = nrow(trades)
  out$avg.pnl = mean(tpnl)
  out$len = mean(tlen)
  out$win.prob = len(which(tpnl > 0)) / out$ntrades
  out$win.avg.pnl = mean(tpnl[tpnl > 0])
  out$win.len = mean(tlen[tpnl > 0])
  out$loss.prob = 1 - out$win.prob
  out$loss.avg.pnl = mean(tpnl[tpnl < 0])
  out$loss.len = mean(tlen[tpnl < 0])
  out$expectancy = (out$win.prob * out$win.avg.pnl + out$loss.prob * 
                      out$loss.avg.pnl) / 100
  out$profitfactor = -(out$win.prob * out$win.avg.pnl) / 
    (out$loss.prob * out$loss.avg.pnl)
  return(as.matrix(unlist(out)))
  }

bt.change.periodicity <- function(b, periodicity = 'months', 
                                  period.ends = NULL, date.map.fn = NULL) {
  require(xts)
  b1 = env()
  for(n in ls(b))
    if(is.xts(b[[n]])) {
      
      if(!is.null(periodicity))
        period.ends = endpoints(b[[n]], periodicity)
      temp = b[[n]][period.ends, ]
      
      if(!is.null(date.map.fn))
        index(temp) = date.map.fn(index(temp))
      colnames(temp) = colnames(b[[n]])
      b1[[n]] = temp
      } else
        b1[[n]] = b[[n]]
  if(!is.null(b$dates))
    b1$dates = index(b1$prices)
  b1
  }

bt.apply <- function(b, xfun = Cl, ...) {
  out = b$weight
  out[] = NA
  symbolnames = b$symbolnames
  nsymbols = length(symbolnames)
  xfun = match.fun(xfun)
  
  for(i in 1:nsymbols) {
    msg = try(xfun(coredata(b[[symbolnames[i]]]), ...), silent = T)
    if(class(msg)[1] == 'try-error')
      warning(i, msg, '\n')
    else
      out[, i] = msg
    }
  return(out)
  }

bt.apply.matrix <- function(b, xfun = Cl, ...) {
  out = b
  out[] = NA
  nsymbols = ncol(b)
  xfun = match.fun(xfun)
  
  for(i in 1:nsymbols) {
    msg = try(xfun(coredata(b[, i]), ...) , silent = T)
    if (class(msg)[1] == 'try-error')
      warning(i, msg, '\n')
    else
      out[, i] = msg
    }
  return(out)
  }

bt.apply.ex <- function(b, xfun = Cl, ..., periodicity = NULL, 
                        period.ends = NULL, apply.periodicity = NULL, 
                        apply.period.ends = NULL, fill.gaps = F) {
  
  temp = bt.apply.setup.helper(b$weight, xfun, periodicity, period.ends, 
                               apply.periodicity, apply.period.ends)
  period.ends = temp$period.ends
  apply.period.ends = temp$apply.period.ends
  map = temp$map
  out = b$weight
  out[] = NA
  symbolnames = b$symbolnames
  nsymbols = length(symbolnames)
  xfun = match.fun(xfun)
  
  if(is.null(apply.period.ends)) {
    if(is.null(period.ends))
      for(i in 1:nsymbols) {
        msg = try(xfun(coredata(b[[symbolnames[i]]]), ...), silent=T)
        if(class(msg)[1] != 'try-error')
          out[, i] = msg
        else
          warning(i, msg, '\n')
        }
    else
      for(i in 1:nsymbols) {
        msg = try(xfun(coredata(b[[symbolnames[i]]][period.ends, ]), ... ), 
                  silent = T)
        if(class(msg)[1] != 'try-error')
          out[period.ends, i] = msg
        else
          warning(i, msg, '\n')
        }
    } else {
      if(is.null(period.ends))
        for(i in 1:nsymbols) {
          x = coredata(b[[symbolnames[i]]])
          for(j in apply.period.ends) {
            msg = try(xfun(x[1:j, , drop = F], ...) , silent = T)
            if(class(msg)[1] != 'try-error')
              out[j, i] = msg
            else
              warning(i, msg, '\n')
            }
          }
      else
        for(i in 1:nsymbols) {
          x = coredata(b[[symbolnames[i]]][period.ends, ])
          for(j in apply.period.ends) {
            msg = try(xfun(x[1:map[j]], ...), silent = T)
            if(class(msg)[1] != 'try-error')
              out[j, i] = msg
            else
              warning(i, msg, '\n')
          }
        }
    }
  if(fill.gaps) bt.apply.matrix(out, ifna.prev) else out
  }

bt.apply.matrix.ex <- function(b, xfun = Cl, ..., periodicity = NULL, 
                               period.ends = NULL, apply.periodicity = NULL, 
                               apply.period.ends = NULL, fill.gaps = F) {
  temp = bt.apply.setup.helper(b, xfun, periodicity, period.ends, 
                               apply.periodicity, apply.period.ends)
  period.ends = temp$period.ends
  apply.period.ends = temp$apply.period.ends
  map = temp$map
  out = b
  out[] = NA
  nsymbols = ncol(b)
  xfun = match.fun(xfun)
  
  if(is.null(apply.period.ends)) {
    if(is.null(period.ends))
      for(i in 1:nsymbols) {
        msg = try(xfun(coredata(b[, i]), ...) , silent = T)
        if(class(msg)[1] != 'try-error')
          out[, i] = msg
        else
          warning(i, msg, '\n')
        }
    else
      for(i in 1:nsymbols) {
        msg = try(xfun( coredata(b[period.ends, i]), ...) , silent = T)
        if(class(msg)[1] != 'try-error')
          out[period.ends,i] = msg
        else
          warning(i, msg, '\n')
        }
    } else {
      if(is.null(period.ends))
        for(i in 1:nsymbols) {
          x = coredata(b[, i])
          for(j in apply.period.ends) {
            msg = try(xfun(x[1:j], ...), silent = T)
            if(class(msg)[1] != 'try-error')
              out[j, i] = msg
            else
              warning(i, msg, '\n')
          }
        }
      else
        for(i in 1:nsymbols) {
          x = coredata(b[period.ends,i])
          for(j in apply.period.ends) {
            msg = try(xfun(x[1:map[j]], ...) , silent = T)
            if(class(msg)[1] != 'try-error')
              out[j, i] = msg
            else
              warning(i, msg, '\n')
          }
        }
      }
  if(fill.gaps) bt.apply.matrix(out, ifna.prev) else out
  }

bt.apply.setup.helper <- function(m, xfun, periodicity, period.ends, 
                                  apply.periodicity, apply.period.ends) {
  if(!is.null(periodicity) && is.null(period.ends))
    period.ends = endpoints(m, periodicity)
  if(!is.null(apply.periodicity) && is.null(apply.period.ends))
    apply.period.ends = endpoints(m, apply.periodicity)
  if(!is.null(apply.period.ends))
    apply.period.ends = apply.period.ends[apply.period.ends > 0]
  if(!is.null(period.ends))
    period.ends = period.ends[period.ends > 0]
  
  map = NULL
  if(!is.null(apply.period.ends) && !is.null(period.ends)) {
    map = array(NA, nrow(m))
    map[period.ends] = 1:len(period.ends)
    map = ifna.prev(map)
    map = ifna(map,1)
    }
  list(period.ends = period.ends, apply.period.ends = apply.period.ends, 
       map = map)
  }

bt.apply.ex2 <- function(b, xfun = Cl, ..., periodicity = NULL, 
                         period.ends = NULL, apply.periodicity = NULL, 
                         apply.period.ends = NULL, fill.gaps = F) {
  temp = bt.apply.setup.helper(b$weight, xfun, periodicity, period.ends, 
                               apply.periodicity, apply.period.ends)
  period.ends = temp$period.ends
  apply.period.ends = temp$apply.period.ends
  map = temp$map
  temp = b$weight
  temp[] = NA
  out = env(out = temp, n = 1, name = 'out')
  index = 1:nrow(temp)
  symbolnames = b$symbolnames
  nsymbols = length(symbolnames)
  xfun = match.fun(xfun)
  
  if(is.null(apply.period.ends)) {
    for(i in 1:nsymbols)
      if(is.null(period.ends))
        set.result.helper(b[[symbolnames[i]]], index, xfun, out, i, ...)
    else
      set.result.helper(b[[symbolnames[i]]][period.ends, ], period.ends, 
                        xfun, out, i, ...)
    } else {
      for(i in 1:nsymbols) {
        x = coredata(iif(is.null(period.ends), b[[symbolnames[i]]], 
                         b[[symbolnames[i]]][period.ends, ]))
        for(j in apply.period.ends)
          if(is.null(period.ends))
            set.result.helper(x[1:j, , drop = F], j, xfun, out, i, ...)
        else
          set.result.helper(x[1:map[j], , drop = F], j, xfun, out, i, ...)
      }
    }
  bt.apply.fill.gaps.helper(out, fill.gaps)
  }

bt.apply.matrix.ex2 <- function(b, xfun = Cl, ..., periodicity = NULL, 
                                period.ends = NULL, apply.periodicity = NULL, 
                                apply.period.ends = NULL, fill.gaps = F) {
  temp = bt.apply.setup.helper(b, xfun, periodicity, period.ends, 
                               apply.periodicity, apply.period.ends)
  period.ends = temp$period.ends
  apply.period.ends = temp$apply.period.ends
  map = temp$map
  temp = b
  temp[] = NA
  out = env(out = temp, n = 1, name = 'out')
  index = 1:nrow(temp)
  nsymbols = ncol(b)
  xfun = match.fun(xfun)
  
  if(is.null(apply.period.ends)) {
    for(i in 1:nsymbols)
      if(is.null(period.ends))
        set.result.helper(b[, i], index, xfun, out, i, ...)
    else
      set.result.helper(b[period.ends, i], period.ends, xfun, out, i, ...)
    } else {
      for(i in 1:nsymbols) {
        x = coredata(iif(is.null(period.ends), b[, i], b[period.ends, i]))
        for(j in apply.period.ends)
          if(is.null(period.ends)) {
            set.result.helper(x[1:j], j, xfun, out, i, ...)
            } else
              set.result.helper(x[1:map[j]], j, xfun, out, i, ...)
      }
    }
  bt.apply.fill.gaps.helper(out, fill.gaps)
  }

set.result.helper = function(x, j, xfun, out, i, ...) {
  msg = try(xfun(iif(is.xts(x), coredata(x), x), ...), silent = T)
  if (class(msg)[1] == 'try-error')
    warning(i, msg, '\n')
  else {
    nresult = iif(is.null(dim(msg)), 1, ncol(msg))
    if(nresult != out$n) {
      temp = out[[out$name[1]]]
      rm(list = out$name, envir = out)
      out$name = iif(is.null(dim(msg)), names(msg), colnames(msg))
      out$n = nresult
      for(result.name in out$name)
        out[[result.name]] = temp
      }
    if(out$n == 1)
      out$out[j,i] = msg
    else
      for(result.name in out$name)
        out[[result.name]][j, i] = iif(len(j) == 1, msg[result.name], 
                                      msg[, result.name])
    }
  }

bt.apply.fill.gaps.helper = function(out, fill.gaps) {
  if(out$n == 1) {
    if(fill.gaps)
      bt.apply.matrix(out$out, ifna.prev)
    else
      out$out
    } else {
      if(fill.gaps)
        for(result.name in out$name)
          out[[result.name]] = bt.apply.matrix(out[[result.name]], ifna.prev)
        rm(list = c('n', 'name'), envir = out)
        out
    }
  }

bt.apply.test = function() {
  load.packages('quantmod,quadprog,corpcor,lpSolve')
  tickers = spl('SPY,QQQ,EEM,IWM,EFA,TLT,IYR,GLD')
  data = env()
  getSymbols(tickers, src = 'yahoo', from = '1980-01-01', env = data, 
             auto.assign = T)
  for(i in ls(data)) data[[i]] = adjustOHLC(data[[i]], use.Adjusted = T)
  bt.prep(data, align = 'remove.na', dates = '1990::')
  t2 = bt.apply.ex2(data, function(x) ATR(HLC(x)))
  prices = data$prices
  t01 = bt.apply.matrix(prices, SMA, 100)
  t02 = bt.apply(data, function(x) SMA(Cl(x), 100))
  t11 = bt.apply.matrix.ex(prices, SMA, 100)
  t12 = bt.apply.ex(data, function(x) SMA(Cl(x), 100))
  t21 = bt.apply.matrix.ex2(prices, SMA, 100)
  t22 = bt.apply.ex2(data, function(x) SMA(Cl(x), 100))
  
  print(all.equal(t01, t02))
  print(all.equal(t01, t11))
  print(all.equal(t01, t12))
  print(all.equal(t01, t21))
  print(all.equal(t01, t22))
  
  t11 = bt.apply.matrix.ex(prices, SMA, 10, periodicity = 'months')
  t12 = bt.apply.ex(data, function(x) SMA(Cl(x), 10), periodicity = 'months')
  t21 = bt.apply.matrix.ex2(prices, SMA, 10, periodicity = 'months')
  t22 = bt.apply.ex2(data, function(x) SMA(Cl(x), 10), 
                     periodicity = 'months')
  
  print(all.equal(t11, t12))
  print(all.equal(t11, t21))
  print(all.equal(t11, t22))
  
  t11 = bt.apply.matrix.ex(prices, function(x) mean(mlast(x, 100)), 
                           apply.periodicity = 'quarters')
  t12 = bt.apply.ex(data, function(x) mean(mlast(Cl(x), 100)), 
                    apply.periodicity = 'quarters')
  t21 = bt.apply.matrix.ex2(prices, function(x) mean(mlast(x, 100)), 
                            apply.periodicity = 'quarters')
  t22 = bt.apply.ex2(data, function(x) mean(mlast(Cl(x), 100)), 
                     apply.periodicity = 'quarters')
  
  print(all.equal(t11, t12))
  print(all.equal(t11, t21))
  print(all.equal(t11, t22))
  
  t11 = bt.apply.matrix.ex(prices, function(x) mean(mlast(x, 10)), 
                           periodicity = 'months', 
                           apply.periodicity = 'quarters')
  t12 = bt.apply.ex(data, function(x) mean(mlast(Cl(x), 10)), 
                    periodicity = 'months', apply.periodicity = 'quarters')
  t21 = bt.apply.matrix.ex2(prices, function(x) mean(mlast(x, 10)), 
                            periodicity = 'months', 
                            apply.periodicity = 'quarters')
  t22 = bt.apply.ex2(data, function(x) mean(mlast(Cl(x), 10)), 
                     periodicity = 'months', apply.periodicity = 'quarters')
  
  print(all.equal(t11, t12))
  print(all.equal(t11, t21))
  print(all.equal(t11, t22))
  
  load.packages('rbenchmark')
  test01 = function() t01 = bt.apply.matrix(prices, SMA, 100)
  test02 = function() t02 = bt.apply(data, function(x) SMA(Cl(x), 100))
  test11 = function() t11 = bt.apply.matrix.ex(prices, SMA, 100)
  test12 = function() t12 = bt.apply.ex(data, function(x) SMA(Cl(x), 100))
  test21 = function() t21 = bt.apply.matrix.ex2(prices, SMA, 100)
  test22 = function() t22 = bt.apply.ex2(data, function(x) SMA(Cl(x), 100))
  
  library(rbenchmark)
  benchmark(test01(), test02(), test11(), test12(), test21(), test22(), 
            columns = c("test", "replications", "elapsed", "relative"), 
            order = "relative", replications = 50)
  }

exrem <- function(x) {
  temp = c(0, ifna(ifna.prev(x), 0))
  itemp = which(temp != mlag(temp))
  x[] = NA
  x[(itemp - 1)] = temp[itemp]
  return(x)
  }

exrem.test <- function() {
  exrem(c(NA, 1, 1, 0, 1, 1, NA, 0))
  }

bt.exrem <- function(weight) {
  bt.apply.matrix(weight, exrem)
  }

bt.test <- function() {
  load.packages('quantmod')
  tickers = spl('SPY')
  data <- new.env()
  getSymbols(tickers, src = 'yahoo', from = '1970-01-01', env = data, 
             auto.assign = T)
  bt.prep(data, align = 'keep.all', dates = '1970::2011')
  prices = data$prices
  data$weight[] = 1
  buy.hold = bt.run(data)
  sma = bt.apply(data, function(x) { SMA(Cl(x), 200) })
  
  data$weight[] = NA
  data$weight[] = iif(prices >= sma, 1, 0)
  sma.cross = bt.run(data, trade.summary = T)
  
  png(filename = 'plot1.png', width = 600, height = 500, units = 'px', 
      pointsize = 12, bg = 'white')
  plotbt.custom.report.part1(sma.cross, buy.hold)
  dev.off()
  
  png(filename = 'plot2.png', width = 1200, height = 800, units = 'px', 
      pointsize = 12, bg = 'white')
  plotbt.custom.report.part2(sma.cross, buy.hold)
  dev.off()
  
  png(filename = 'plot3.png', width = 600, height = 500, units = 'px', 
      pointsize = 12, bg = 'white')
  plotbt.custom.report.part3(sma.cross, buy.hold)
  dev.off()
  
  pdf(file = 'report.pdf', width=8.5, height = 11)
  plotbt.custom.report(sma.cross, buy.hold, trade.summary = T)
  dev.off()
  
  data$weight[] = NA
  data$weight$SPY = 1
  temp = bt.run(data)
  data$weight[] = NA
  data$weight$SPY = 2
  temp = bt.run(data)
  data$weight[] = NA
  data$weight$SPY = 1
  capital = 100000
  data$weight[] = (capital / prices) * data$weight
  temp = bt.run(data, type = 'share', capital = capital)
  data$weight[] = NA
  data$weight$SPY = 2
  capital = 100000
  data$weight[] = (capital / prices) * data$weight
  temp = bt.run(data, type = 'share', capital = capital)
  }

compute.cagr <- function(equity, nyears = NA) {
  if(is.numeric(nyears))
    as.double(last(equity, 1)^(1/nyears) - 1)
  else
    as.double(last(equity, 1)^(1/compute.nyears(equity)) - 1)
  }

compute.nyears <- function(x) {
  as.double(diff(as.Date(range(index.xts(x)))))/365
  }

compute.raw.annual.factor = function(x) {
  round(nrow(x) / compute.nyears(x))
  }

compute.annual.factor = function(x) {
  possible.values = c(252, 52, 26, 13, 12, 6, 4, 3, 2, 1)
  index = which.min(abs(compute.raw.annual.factor(x) - possible.values))
  round(possible.values[index])
  }

compute.sharpe <- function(x) {
  temp = compute.annual.factor(x)
  x = as.vector(coredata(x))
  return(sqrt(temp) * mean(x) / sd(x))
  }

compute.calmar <- function(x) {
  compute.cagr(x) / compute.max.drawdown(x)
  }

compute.R2 <- function(equity) {
  x = as.double(index.xts(equity))
  y = equity
  return(cor(y, x)^2)
  }

compute.DVR <- function(bt) {
  return(compute.sharpe(bt$ret) * compute.R2(bt$equity))
  }

compute.risk <- function(x) {
  temp = compute.annual.factor(x)
  x = as.vector(coredata(x))
  return(sqrt(temp) * sd(x))
  }

compute.drawdown <- function(x) {
  return(x / cummax(c(1, x))[-1] - 1)
  }

compute.max.drawdown <- function(x) {
  as.double(min(compute.drawdown(x)))
  }

compute.avg.drawdown <- function(x) {
  drawdown = c(0, compute.drawdown(coredata(x)), 0)
  dstart = which(drawdown == 0 & mlag(drawdown, -1) != 0)
  dend = which(drawdown == 0 & mlag(drawdown, 1) != 0)
  drawdowns = apply(cbind(dstart, dend), 1, function(x) 
    min(drawdown[x[1]:x[2]], na.rm = T))
  mean(drawdowns)
  }

compute.cdar <- function(x, probs = 0.05) {
  drawdown = c(0, compute.drawdown(coredata(x)), 0)
  dstart = which(drawdown == 0 & mlag(drawdown, -1) != 0)
  dend = which(drawdown == 0 & mlag(drawdown, 1) != 0)
  drawdowns = apply(cbind(dstart, dend), 1, function(x) 
    min(drawdown[x[1]:x[2]], na.rm = T))
  
  if(len(drawdowns) > 2)
    mean(drawdowns[drawdowns < quantile(drawdowns, probs=probs)])
  else
    min(drawdowns)
  }

compute.exposure <- function(weight) {
  sum(apply(weight, 1, function(x) sum(x != 0)) != 0) / nrow(weight)
  }

compute.var <- function(x, probs = 0.05) {
  quantile(coredata(x), probs = probs)
  }

compute.cvar <- function(x, probs = 0.05) {
  x = coredata(x)
  mean(x[x < quantile(x, probs=probs)])
  }

compute.stats <- function(data, fns, do.na.omit = T) {
  out = matrix(double(), len(fns), len(data))
  colnames(out) = names(data)
  rownames(out) = names(fns)
  if(do.na.omit)
    for(c in 1:len(data)) {
      for(r in 1:len(fns)) {
        out[r, c] = match.fun(fns[[r]])(fast.na.omit(data[[c]]))
      }
    }
  else
    for(c in 1:len(data)) {
      for(r in 1:len(fns)) {
        out[r, c] = match.fun(fns[[r]])(data[[c]])
      }
    }
  return(out)
  }

bt.simple <- function(data, signal, silent = F) {
  signal = Lag(signal, 1)
  signal = na.locf(signal, na.rm = FALSE)
  signal[is.na(signal)] = 0
  ret = ROC(Cl(data), type = 'discrete')
  ret[1] = 0
  n = nrow(ret)
  
  bt <- list()
  bt$ret = ret * signal
  bt$best = max(bt$ret)
  bt$worst = min(bt$ret)
  bt$equity = cumprod(1 + bt$ret)
  bt$cagr = bt$equity[n] ^ (1/nyears(data)) - 1
  if(!silent) {
    cat('', spl('CAGR,Best,Worst'), '\n', sep = '\t')
    cat('', sapply(cbind(bt$cagr, bt$best, bt$worst), function(x) 
      round(100 * x, 1)), '\n', sep = '\t')
    }
  return(bt)
  }

bt.simple.test <- function() {
  load.packages('quantmod')
  data = getSymbols('SPY', src = 'yahoo', from = '1980-01-01', 
                    auto.assign = F)
  signal = rep(1, nrow(data))
  buy.hold = bt.simple(data, signal)
  sma = SMA(Cl(data), 200)
  signal = ifelse(Cl(data) > sma, 1, 0)
  sma.cross = bt.simple(data, signal)
  dates = '2000::2009'
  buy.hold.equity <- buy.hold$equity[dates] / 
    as.double(buy.hold$equity[dates][1])
  sma.cross.equity <- sma.cross$equity[dates] / 
    as.double(sma.cross$equity[dates][1])
  
  chartSeries(buy.hold.equity, TA = c(addTA(sma.cross.equity, on = 1, 
                                            col = 'red')), theme = 'white', 
              yrange = range(buy.hold.equity, sma.cross.equity))
  }

bt.apply.min.weight <- function(weight, long.min.weight = 0.1, 
                                short.min.weight = long.min.weight) {
  if(is.null(dim(weight))) dim(weight) = c(1, len(weight))
  pos = apply(weight, 1, function(row) sum(row[row > 0]))
  neg = rowSums(weight) - pos
  
  pos.mat = iif(weight >= long.min.weight, weight, 0)
  neg.mat = iif(weight <= -short.min.weight, weight, 0)
  pos.mat = pos.mat * ifna(pos / rowSums(pos.mat), 1)
  neg.mat = neg.mat * ifna(neg / rowSums(neg.mat), 1)
  return(pos.mat + neg.mat)
  }

test.bt.apply.min.weight <- function() {
  data = c(0.1, 0.6, 0.2, 0.1, 0, -0.1, -0.6, -0.2, -0.1, 0)
  mm = matrix(data = data, nrow = 2, byrow = TRUE)
  print(bt.apply.min.weight(mm, 0.1))
  print(bt.apply.min.weight(mm, 0.2))
  
  data = c(0.1, 0.6, 0.2, 0.1, 0, -0.1, -0.6, -0.2, -0.1, 0)
  mm = matrix(data = data, nrow = 1, byrow = TRUE)
  print(bt.apply.min.weight(mm, 0.1))
  print(bt.apply.min.weight(mm, 0.2))
  
  data = c(0.1, 0.6, 0.2, 0.1, 0, -0.2, -0.5, -0.3, -0.1, 0)
  mm = matrix(data = data, nrow = 1, byrow = TRUE)
  print(bt.apply.min.weight(mm, 0.1))
  print(bt.apply.min.weight(mm, 0.2))
  }

bt.apply.round.weight <- function(weight, long.round.weight = 5/100, 
                                  short.round.weight = long.round.weight) {
  if(is.null(dim(weight))) dim(weight) = c(1, len(weight))
  
  pos = apply(weight, 1, function(row) sum(row[row > 0]))
  neg = rowSums(weight) - pos
  pos.mat = iif(weight >= 0, round(weight / long.round.weight) * 
                  long.round.weight, 0)
  neg.mat = iif(weight <= 0, round(weight / short.round.weight) * 
                  short.round.weight, 0)
  pos.mat = pos.mat * ifna(pos / rowSums(pos.mat), 1)
  neg.mat = neg.mat * ifna(neg / rowSums(neg.mat), 1)
  return(pos.mat + neg.mat)
  }

bt.start.dates <- function(b) {
  temp = lapply(b, function(x) index(x[1]))
  temp$dates = NULL
  temp$prices = NULL
  temp$weight = NULL
  temp$execution.price = NULL
  temp$symbolnames = NULL
  temp = temp[order(sapply(temp, function(x) x))]
  out = t(t(sapply(temp, function(x) as.character(x))))
  colnames(out) = 'Start'
  out
  }

bt.end.dates <- function(b) {
  temp = lapply(b, function(x) index(last(x)))
  temp$dates = NULL
  temp$prices = NULL
  temp$weight = NULL
  temp$execution.price = NULL
  temp$symbolnames = NULL
  temp = temp[order(sapply(temp, function(x) x))]
  out = t(t(sapply(temp, function(x) as.character(x))))
  colnames(out) = 'Start'
  out
  }

bt.append.today <- function(b, data.today) {
  date.column = find.names('Date', data.today)
  valid.index = which(!is.na(data.today[, date.column, with = F]))
  data.today = data.today[valid.index]
  data = make.stock.xts(read.xts(data.today, date.column = date.column, 
                                 format = '%m/%d/%Y', decreasing = NULL))
  tickers = data.today$Symbol
  Yesterday = data.today$Yesterday
  
  for(i in 1:len(tickers)) {
    if(is.null(b[[tickers[i]]])) next
    if(last(index(data[i, ])) > last(index(b[[tickers[i]]])))
      b[[tickers[i]]] = rbind(data[i, ], b[[tickers[i]]])
    }
  }

bt.run.share.ex <- function(b, prices = b$prices, clean.signal = T, 
                            trade.summary = F, do.lag = 1, silent = F, 
                            capital = 100000, commission = 0, 
                            weight = b$weight, dates = 1:nrow(b$prices), 
                            lot.size = c(), control = list(
                              round.lot = default.round.lot.control()), 
                            cashflow.control = list(), adjusted = T, 
                            dividend.control = list(
                              foreign.withholding.tax = NULL, 
                              invest = spl('cash,rebalance,update')), 
                            tax.control = NULL) {
  
  invest.type.def = list(none = 0, cash = 1, rebalance = 2, update = 3)
  cashflow.type.def = list(regular = 0, fee.rebate = 1)
  dummy = NA * prices[, 1]
  cashflows = list(n = len(cashflow.control), cash = dummy) 
  
  if(cashflows$n > 0) {
    cashflows$last.index = rep(0, cashflows$n)
    cashflows$type = rep('', cashflows$n)
    
    for(i in 1:cashflows$n) cashflows$type[i] = ifnull(
      cashflow.control[[i]]$type, 'regular')[1]
    
    cashflows$type = cashflow.type.def[[tolower(cashflows$type)]]
    cashflows$invest = rep('', cashflows$n)
    
    for(i in 1:cashflows$n) cashflows$invest[i] = 
      ifnull(cashflow.control[[i]]$invest, 'cash')[1]
    
    cashflows$invest = invest.type.def[[tolower(cashflows$invest)]]
    cashflows$fn = lapply(cashflow.control, function(x) 
      iif(is.null(x$cashflow.fn), x$cashflow.fn, match.fun(x$cashflow.fn)))
    
    cashflows$cash = matrix(NA, nrow(prices), cashflows$n)
    colnames(cashflows$cash) = names(cashflow.control)
    
    for(i in 1:cashflows$n)
      if(!is.null(cashflow.control[[i]]$cashflows)) {
        dummy[] = NA
        dummy[index(cashflow.control[[i]]$cashflows)] = 
          cashflow.control[[i]]$cashflows
        cashflows$cash[, i] = dummy
      }
    }
  
  cashflows$cash = coredata(ifna(cashflows$cash, 0))
  weight[is.nan(weight) | is.infinite(weight)] = NA
  weight[!is.na(weight) & is.na(prices)] = 0
  weight = iif(do.lag == 1, weight, mlag(weight, do.lag - 1))
  weight = coredata(weight)
  temp = bt.exrem(weight)
  
  if(clean.signal) {
    weight = temp
    } else {
      index = ifna(weight == 0, F)
      weight[index] = temp[index]
      }
  
  prices = coredata(prices)
  n = ncol(prices)
  trade = !is.na(weight)
  
  if(sum(trade) > 0) {
    execution.price = coredata(b$execution.price)
    prices[trade] = iif(is.na(execution.price[trade]), prices[trade], 
                        execution.price[trade])
  }
  
  prices[] = ifna( bt.apply.matrix(prices, ifna.prev), 1)
  if(!is.list(commission)) 
    commission = list(cps = commission, fixed = 0.0, percentage = 0.0)
  
  lot.size = map2vector(lot.size, colnames(prices), 1)
  dividend.control$foreign.withholding.tax = 
    map2vector(dividend.control$foreign.withholding.tax, colnames(prices), 0) 
  dividend.control$invest = ifnull(dividend.control$invest, 'cash')
  dividend.control$invest = 
    invest.type.def[[tolower(dividend.control$invest[1])]]
  
  bt = bt.run.share.ex.internal(b, prices, capital, commission, weight, 
                                lot.size, control, cashflows, adjusted, 
                                dividend.control, tax.control)
  bt$ret = make.xts(bt$ret, index(b$prices))
  bt$weight = make.xts(bt$weight, index(b$prices))
  bankrupt = which(bt$ret <= -1)
  
  if(len(bankrupt) > 0) bt$ret[bankrupt[1]:n] = -1
  bt$equity = cumprod(1 + bt$ret)
  bt$dates.index = dates2index(b$prices, dates)
  bt = bt.run.trim.helper(bt, bt$dates.index)
  
  if(trade.summary) bt$trade.summary = bt.trade.summary(b, bt)
  if(!silent) {
    cat('Latest weights :\n')
    print(round(100 * last(bt$weight), 2))
    cat('\n')
    cat('Performance summary :\n')
    cat('', spl('CAGR,Best,Worst'), '\n', sep = '\t')
    cat('', sapply(cbind(bt$cagr, bt$best, bt$worst), function(x) 
      round(100 * x, 1)), '\n', sep = '\t')
    cat('\n')
    }
  bt
  }

default.tax.control = function(nonqualified = c()) {
  list(capital.gain = list(short.term.tax = 35/100, long.term.tax = 15/100, 
                           wash.sale.min.holding.period = 30, 
                           long.term.min.holding.period = 365), 
       dividend = list(qualified.tax = 20/100, 
                       qualified.min.holding.period = 60, 
                       nonqualified.tax = 35/100, 
                       nonqualified = nonqualified))
  }

bt.run.share.ex.internal <- function(b, prices, capital, commission, weight, 
                                     lot.size, control, cashflows, adjusted, 
                                     dividend.control, tax.control) {
  invest.type.def = list(none = 0, cash = 1, rebalance = 2, update = 3)
  cashflow.type.def = list(regular = 0, fee.rebate = 1)
  weight1 = ifna(bt.apply.matrix(weight, ifna.prev), 0)
  trade.cashflow = rowSums(cashflows$cash != 0) > 0
  trade = !is.na(weight)
  trade.index = rowSums(trade) > 0
  
  if(!adjusted) {
    dividends = coredata(bt.apply(b, function(x) x[, 'Dividend']))
    splits = coredata(bt.apply(b, function(x) x[, 'Split']))
    trade.dividend = rowSums(mlag(weight1) != 0 & dividends > 0, 
                             na.rm = T) > 0
    trade.split = rowSums(mlag(weight1) != 0 & splits > 0, na.rm = T) > 0
    event.index = which(trade.index | trade.cashflow | trade.dividend | 
                          trade.split)
    } else
      event.index = which(trade.index | trade.cashflow)
  nperiods = nrow(prices)
  n = ncol(prices)
  
  if(!is.null(tax.control)) {
    holdings = env(n.trades = rep(0, n), share = matrix(0, 100, n), 
                   date = matrix(0, 100, n), price = matrix(0, 100, n))
    
    wash.sale = env(n.trades = rep(0, n), date = matrix(0, 100, n), 
                    share = matrix(0, 100, n), 
                    loss.per.share = matrix(0, 100, n), 
                    long.term = matrix(T, 100, n))
    
    tax = env(long.term.cap = rep(0, nperiods), 
              short.term.cap = rep(0, nperiods), 
              qualified.div = rep(0, nperiods), 
              non.qualified.div = rep(0, nperiods))
    
    tax.control$dividend$nonqualified = 
      map2vector(tax.control$dividend$nonqualified, colnames(prices), F)
    }
  
  cash.wt = rep(capital, nperiods)
  event.type = div = com = cashflow = fee.rebate = rep(0, nperiods)
  event.type.def = list(none = 0, trade = 1, split = 2, dividend = 3, 
                        cashflow = 4)
  share.wt = matrix(0, nperiods, n)
  colnames(share.wt) = colnames(prices)
  info = env(cash = cash.wt, share = share.wt)
  
  if(!is.null(tax.control)) {
    info$dates = b$dates
    info$tax = tax
    info$tax.control = tax.control
    info$holdings = holdings
    }
  
  last.trade = 0
  weight.last = weight1[1, ]
  for(i in event.index) {
    trade.invest.type = iif(trade.index[i], invest.type.def$rebalance, 
                            invest.type.def$none)
    trade.today = trade.index[i]
    
    if(last.trade > 0) {
      index = (last.trade + 1) : i
      n.index = len(index)
      share.wt[index,] = rep.row(info$share[last.trade, ], n.index)
      info$share[index,] = rep.row(info$share[last.trade,], n.index)
      cash.wt[index] = info$cash[last.trade]
      info$cash[index] = info$cash[last.trade]
      weight.last = weight1[i - 1, ]
      }
    
    if(!adjusted) {
      if(trade.dividend[i]) {
        if(!is.null(tax.control)) {
          tax.update.dividends(tax.control, dividend.control, holdings, tax, 
                               info$share[i, ], dividends[i, ], i, b$dates)
          }
        
        asset.cashflow = 
          sum(info$share[i, ] * dividends[i, ] * 
                iif(info$share[i, ] < 0, 1, 1 - 
                      dividend.control$foreign.withholding.tax))
        
        info$cash[i] = info$cash[i] + asset.cashflow
        cash.wt[i] = cash.wt[i] + asset.cashflow
        div[i] = asset.cashflow
        event.type[i] = event.type.def$dividend
        
        if(dividend.control$invest == invest.type.def$rebalance | 
           dividend.control$invest == invest.type.def$update) {
          trade.index[i] = T
          if(trade.invest.type == invest.type.def$none) 
            trade.invest.type = dividend.control$invest
        }
      }
      
      if(trade.split[i]) {
        for(a in which(info$share[i, ] != 0 & splits[i, ] > 0))
          info$share[i, a] = share.wt[i, a] = info$share[i, a] / splits[i, a]
        event.type[i] = event.type.def$split
        
        if(!is.null(tax.control)) {
          for(a in which(info$share[i, ] !=0 & splits[i, ] > 0)) {
            n.trades = 1:holdings$n.trades[a]
            holdings$share[n.trades, a] = holdings$share[n.trades, a] / 
              splits[i, a]
            holdings$price[n.trades, a] = holdings$price[n.trades, a] * 
              splits[i, a]
          }
          
          for(a in which(wash.sale$n.trades > 0 & splits[i, ] > 0)) {
            n.trades = 1:wash.sale$n.trades[a]
            wash.sale$share[n.trades, a] = wash.sale$share[n.trades, a] / 
              splits[i, a]
            wash.sale$loss.per.share[n.trades, a] = 
              wash.sale$loss.per.share[n.trades, a] * splits[i, a]
          }
        }
      }
    }
    
    if(trade.cashflow[i]) {
      for(c in (1:cashflows$n)[cashflows$cash[i, ] != 0]) {
        if(!is.null(cashflows$fn[[c]])) {
          cashflows$cash[i, c] = cashflows$fn[[c]](info, i, 
                                                   cashflows$last.index[c])
          }
        
        info$cash[i] = info$cash[i] + cashflows$cash[i, c]
        if(cashflows$type[c] == cashflow.type.def$regular)
          cashflow[i] = cashflow[i] + cashflows$cash[i, c]
        else
          fee.rebate[i] = fee.rebate[i] + cashflows$cash[i, c]
        
        event.type[i] = event.type.def$cashflow
        cashflows$last.index[c] = i
        if(cashflows$invest[c] == invest.type.def$rebalance | 
           cashflows$invest[c] == invest.type.def$update) {
          trade.index[i] = T
          
          if(trade.invest.type == invest.type.def$none) 
            trade.invest.type = cashflows$invest[c]
        }
      }
    }
    
    if(trade.index[i]) {
      weight.change.index = iif(trade.today, !is.na(weight[i, ]), rep(T, n)) 
      if(trade.invest.type == invest.type.def$rebalance)
        out = bt.run.share.ex.allocate(
          weight.new = weight1[i, ], weight.prev = weight.last, 
          weight.change.index = weight.change.index, 
          price = prices[i, ], share = info$share[i, ], cash = info$cash[i], 
          commission, lot.size, control = control$round.lot, 
          cashflow = cashflow[i] + fee.rebate[i] + div[i])
      
      if(trade.invest.type == invest.type.def$update) {
        out = bt.run.share.ex.invest(
          weight.new = weight1[i, ], weight.prev = weight.last, 
          weight.change.index = weight.change.index, 
          price = prices[i, ], share = info$share[i, ], cash = info$cash[i], 
          cashflow = cashflow[i] + fee.rebate[i] + div[i], commission, 
          lot.size, control = control$round.lot)
      }
      
      if(!is.null(tax.control) && sum(info$share[i,] != out$share) > 0) {
        if(any(info$share[i, ] != sapply(1:n, function(a) 
          sum(iif(holdings$n.trades[a] > 0, 
                  holdings$share[1:holdings$n.trades[a], a], 0)))))
          
          cat('Wrong Holdings', info$share[i, ][5], '\n')
        tax.update.holdings(tax.control, holdings, tax, wash.sale, 
                            info$share[i, ], out$share, prices[i, ], i, 
                            b$dates)
        } 
      info$share[i,] = out$share
      info$cash[i] = out$cash
      com[i] = out$com
      event.type[i] = event.type.def$trade
    }
    last.trade = i
  }
  
  if(last.trade > 0 & last.trade < nperiods) {
    index = (last.trade + 1):nperiods
    n.index = len(index)
    share.wt[index,] = rep.row(info$share[last.trade, ], n.index)
    info$share[index, ] = rep.row(info$share[last.trade, ], n.index)
    cash.wt[index] = info$cash[last.trade]
    info$cash[index] = info$cash[last.trade]
  }
  
  bt = list(type = 'share', capital = capital, share = info$share, 
            cash = info$cash, value = info$cash + 
              rowSums(info$share * prices), com = com, div = div, 
            weight = share.wt * prices / 
              (cash.wt + rowSums(share.wt * prices)), 
            event.type = factor(event.type, 
                                as.character(unlist(event.type.def)), 
                                names(event.type.def)))
  
  if(cashflows$n > 0) {
    bt$cashflows = cashflows$cash
    bt$cashflow = cashflow
    bt$fee.rebate = fee.rebate
    }
  
  if(!is.null(tax.control)) {
    bt$long.term.cap = tax$long.term.cap
    bt$short.term.cap = tax$short.term.cap
    bt$qualified.div = tax$qualified.div
    bt$non.qualified.div = tax$non.qualified.div
  }
  
  value = c(capital, bt$value)
  bt$ret = (value / mlag(value) - 1)[-1]
  
  if(sum(abs(cashflow)) > 0) {
    index = cashflow < 0
    bt$ret[index] = (c(capital, bt$value - cashflow) / 
                       mlag(value) - 1)[-1][index]
    index = cashflow > 0
    value1 = c(capital, bt$value + mlag(cashflow, -1))
    bt$ret[index] = (value / mlag(value1) - 1)[-1][index]
    }
  bt
  }

bt.run.share.ex.n.days = function(index.today, index.hist, dates) {
  as.numeric(dates[index.today] - dates[index.hist])
  }

tax.update.dividends = function(tax.control, dividend.control, holdings, tax, 
                                share, dividend, index, dates) {
  qualified.div = dividend * (1 - dividend.control$foreign.withholding.tax / 
                                tax.control$dividend$qualified.tax)
  non.qualified.div = dividend * 
    (1 - dividend.control$foreign.withholding.tax / 
       tax.control$dividend$nonqualified.tax)
  qualified.div.adj = non.qualified.div.adj = 0
  pos.cashflow = share > 0 & dividend > 0
  
  for(a in which(pos.cashflow & !tax.control$dividend$nonqualified)) {
    n.trades = 1:holdings$n.trades[a]
    
    trade.days = bt.run.share.ex.n.days(index, holdings$date[n.trades, a], 
                                        dates)
    
    nonqualified.trade.days = trade.days < 
      tax.control$dividend$qualified.min.holding.period
    
    if(sum(nonqualified.trade.days) == 0) next
    qualified.div.adj = qualified.div.adj - qualified.div[a] * 
      sum(holdings$share[n.trades, a][nonqualified.trade.days])
    
    non.qualified.div.adj = non.qualified.div.adj + non.qualified.div[a] * 
      sum(holdings$share[n.trades, a][nonqualified.trade.days])
    }
  
  tax$qualified.div[index] = tax$qualified.div[index] + qualified.div.adj + 
    sum((share * qualified.div)[pos.cashflow & 
                                  !tax.control$dividend$nonqualified])
  
  tax$non.qualified.div[index] = tax$non.qualified.div[index] + 
    non.qualified.div.adj + sum(
      (share * non.qualified.div)[pos.cashflow & 
                                    tax.control$dividend$nonqualified])
  
  for(a in which(share < 0 & dividend > 0)) {
    n.trades = 1:holdings$n.trades[a]
    holdings$price[n.trades, a] = holdings$price[n.trades, a] - dividend[a]
    }
  }

record.wash.sale = function(a, n.trades, pnl, price, trigger, holdings, 
                            wash.sale) {
  wash.sale.index = pnl[n.trades] < 0
  n.wash.sale.index = sum(wash.sale.index)
  
  if(n.wash.sale.index > 0) {
    if(wash.sale$n.trades[a] + n.wash.sale.index > nrow(wash.sale$share)) {
      n = ncol(wash.sale$date)
      
      wash.sale$date = rbind(wash.sale$date, matrix(0, 100, n))
      wash.sale$share = rbind(wash.sale$share, matrix(0, 100, n))
      
      wash.sale$loss.per.share = rbind(wash.sale$loss.per.share, 
                                       matrix(0, 100, n))
      
      wash.sale$long.term = rbind(wash.sale$long.term, matrix(T, 100, n))
      }
    
    n1 = wash.sale$n.trades[a] + 1:n.wash.sale.index
    wash.sale$date[n1, a] = holdings$date[n.trades, a][wash.sale.index]
    wash.sale$share[n1, a] = holdings$share[n.trades, a][wash.sale.index]
    
    wash.sale$loss.per.share[n1, a] = 
      (price[a] - holdings$price[n.trades, a])[wash.sale.index]
    
    wash.sale$long.term[n1,a] = trigger[n.trades][wash.sale.index]
    wash.sale$n.trades[a] = wash.sale$n.trades[a] + n.wash.sale.index
    }
  }

check.wash.sale = function(a, dates, tax, tax.control, holdings, wash.sale) {
  if(wash.sale$n.trades[a] == 0) return()
  i = holdings$n.trades[a]
  index = holdings$date[i, a]
  n.trades = 1:wash.sale$n.trades[a]
  
  trade.days = bt.run.share.ex.n.days(index , wash.sale$date[n.trades, a], 
                                      dates)
  trigger = trade.days <= 
    tax.control$capital.gain$wash.sale.min.holding.period
  n.trigger = sum(trigger)
  
  if(n.trigger == 0) {
    wash.sale$n.trades[a] = 0
    return()
  }
  share1 = abs(holdings$share[i, a])
  run.share = cumsum(abs(wash.sale$share[n.trades, a][trigger]))
  n1 = which( run.share > share1)
  
  if(len(n1) == 0 || mlast(run.share) == share1) {
    loss = wash.sale$loss.per.share[n.trades, a][trigger] * 
      wash.sale$share[n.trades, a][trigger]
    
    holdings$price[i, a] = holdings$price[i, a] - sum(loss) / 
      holdings$share[i, a]
    wash.sale$n.trades[a] = 0
    
    tax$long.term.cap[index] = tax$long.term.cap[index] - 
      sum(iif(wash.sale$long.term[n.trades, a][trigger], loss, 0))
    
    tax$short.term.cap[index] = tax$short.term.cap[index] - 
      sum(iif(wash.sale$long.term[n.trades, a][trigger], 0, loss))
    return()
  }
  
  n1 = n1[1]
  trigger.index = which(trigger)
  
  if(run.share[n1] == share1) {
    loss = wash.sale$loss.per.share[n.trades, a][trigger][1:n1] * 
      wash.sale$share[n.trades, a][trigger][1:n1]
    
    holdings$price[i, a] = holdings$price[i, a] - sum(loss) / 
      holdings$share[i, a]
    
    tax$long.term.cap[index] = tax$long.term.cap[index] - 
      sum(iif(wash.sale$long.term[n.trades, a][trigger][1:n1], loss, 0))
    
    tax$short.term.cap[index] = tax$short.term.cap[index] - 
      sum(iif(wash.sale$long.term[n.trades, a][trigger][1:n1], 0, loss))
    n1 = n1 + 1
    } else {
      share.left = run.share[n1] - share1
      last.index = trigger.index[n1]
      
      wash.sale$share[last.index, a] = wash.sale$share[last.index, a] - 
        share.left
      loss = wash.sale$loss.per.share[n.trades, a][trigger][1:n1] * 
        wash.sale$share[n.trades, a][trigger][1:n1]
      
      holdings$price[i, a] = holdings$price[i, a] - sum(loss) / 
        holdings$share[i, a]
      
      tax$long.term.cap[index] = tax$long.term.cap[index] - 
        sum(iif(wash.sale$long.term[n.trades, a][trigger][1:n1], loss, 0))
      
      tax$short.term.cap[index] = tax$short.term.cap[index] - 
        sum(iif(wash.sale$long.term[n.trades, a][trigger][1:n1], 0, loss))
      wash.sale$share[last.index, a] = share.left
    }
  
  if(n1 > 1) {
    from.index = trigger.index[n1]:wash.sale$n.trades[a]
    to.index = 1:len(from.index)
    
    wash.sale$date[to.index, a] = wash.sale$date[from.index, a]
    wash.sale$share[to.index, a] = wash.sale$share[from.index, a]
    
    wash.sale$loss.per.share[to.index, a] = 
      wash.sale$loss.per.share[from.index, a]
    
    wash.sale$long.term[to.index, a] = wash.sale$long.term[from.index, a]
    wash.sale$n.trades[a] = len(to.index)
  }
}

tax.update.holdings = function(tax.control, holdings, tax, wash.sale, share0, 
                               share1, price, index, dates) {
  n = len(price)
  if(any(share0 != sapply(1:n, function(a) 
    sum(iif(holdings$n.trades[a] > 0, holdings$share[1:holdings$n.trades[a], 
                                                     a], 0)))))
    cat('Mismatch holding shares', index, '\n')
  
  for(a in (1:n)[share0 != share1]) {
    n.trades = 1:holdings$n.trades[a]
    if( sign(share0[a] * share1[a]) <= 0) {
      if(share0[a] != 0) {
        pnl = holdings$share[n.trades,a] * 
          (price[a] - holdings$price[n.trades, a])
        
        trade.days = 
          bt.run.share.ex.n.days(index, holdings$date[n.trades, a], dates)
        
        trigger = trade.days > 
          tax.control$capital.gain$long.term.min.holding.period
        
        tax$long.term.cap[index] = tax$long.term.cap[index] + 
          sum(iif(trigger, pnl, 0))
        
        tax$short.term.cap[index] = tax$short.term.cap[index] + 
          sum(iif(trigger, 0, pnl))
        
        record.wash.sale(a, n.trades, pnl, price, trigger, holdings, 
                         wash.sale)
        holdings$n.trades[a] = 0
      }
      
      if(share1[a] != 0) {
        holdings$share[1, a] = share1[a]
        holdings$price[1, a] = price[a]
        holdings$date[1, a] = index
        holdings$n.trades[a] = 1
        check.wash.sale(a, dates, tax, tax.control, holdings, wash.sale)
        
        } else
          holdings$n.trades[a] = 0
        } else {
          
          if(abs(share1[a]) > abs(share0[a])) {
            if(holdings$n.trades[a] + 1 > nrow(holdings$share)) {
              holdings$share = rbind(holdings$share, matrix(0, 100, n))
              holdings$price = rbind(holdings$price, matrix(0, 100, n))
              holdings$date  = rbind(holdings$date, matrix(0, 100, n))
            }
            
            n1 = holdings$n.trades[a] + 1
            holdings$share[n1,a] = share1[a] - share0[a]
            holdings$price[n1,a] = price[a]
            holdings$date[n1, a]  = index
            holdings$n.trades[a] = n1
            check.wash.sale(a, dates, tax, tax.control, holdings, wash.sale)
            }
          
          if(abs(share1[a]) < abs(share0[a])) {
            remove.share = share0[a] - share1[a]
            pnl = holdings$share[n.trades, a] * 
              (price[a] - holdings$price[n.trades, a])
            
            trade.days = 
              bt.run.share.ex.n.days(index, holdings$date[n.trades, a], 
                                     dates)
            
            trigger = trade.days > 
              tax.control$capital.gain$long.term.min.holding.period
            
            run.share = cumsum(holdings$share[n.trades, a])
            n1 = which(abs(run.share) >= abs(remove.share))[1]
            
            if(run.share[n1] == remove.share) {
              tax$long.term.cap[index] = tax$long.term.cap[index] + 
                sum(iif(trigger, pnl, 0)[1:n1])
              
              tax$short.term.cap[index] = tax$short.term.cap[index] + 
                sum(iif(trigger, 0, pnl)[1:n1])
              
              record.wash.sale(a, 1:n1, pnl, price, trigger, holdings, 
                               wash.sale)
              n1 = n1 + 1
              } else {
                share.left = run.share[n1] - remove.share
                pnl[n1] = pnl[n1] - (holdings$share[n1, a] - share.left) * 
                  (price[a] - holdings$price[n1, a])
                
                tax$long.term.cap[index] = tax$long.term.cap[index] + 
                  sum(iif(trigger, pnl, 0)[1:n1])
                
                tax$short.term.cap[index] = tax$short.term.cap[index] + 
                  sum(iif(trigger, 0, pnl)[1:n1])
                
                holdings$share[n1,a] = holdings$share[n1, a] - share.left
                record.wash.sale(a, 1:n1, pnl, price, trigger, holdings, 
                                 wash.sale)
                
                holdings$share[n1,a] = share.left
                }
            
            if(n1 > 1) {
              from.index = n1:holdings$n.trades[a]
              to.index = 1:len(from.index)
              holdings$share[to.index, a] = holdings$share[from.index, a]
              holdings$price[to.index, a] = holdings$price[from.index, a]
              holdings$date[to.index, a] = holdings$date[from.index, a]
              holdings$n.trades[a] = len(to.index)
            }
          }
        }
    if( share1[a] != sum(iif(holdings$n.trades[a] > 0, 
                             holdings$share[1:holdings$n.trades[a], a], 0)))
      cat('a', a, index, '\n')
    }
  
  if( any(share1 != sapply(1:n, function(a) 
    sum(iif(holdings$n.trades[a] > 0, holdings$share[1:holdings$n.trades[a], 
                                                     a], 0)))))
    cat('Mismatch holding shares', index, '\n')
  }

default.round.lot.control = function() {
  list(select = c('best.match', 'minimum.turnover'), diff.target = 5/100)
  }

bt.run.share.ex.invest = function(weight.new, weight.prev, 
                                  weight.change.index, price, share, cash, 
                                  cashflow, commission, lot.size, silent = T, 
                                  control = default.round.lot.control()) {
  if(cashflow >= 0) {
    return(list(share = share, cash = cash, com = 0))
    } else {
      current.cash = sum(price * share) + cash - sum(price * abs(share))
      current.cash = (sum(price * share) + cash) - 
        sum((price * share)[share > 0])
      
      if(current.cash >= 0)
        return(list(share = share, cash = cash, com = 0))
      }
  
  if(F) {
    n = len(share)
    out = bt.run.share.ex.allocate(weight.new = weight.new, 
                                   weight.prev = rep(0, n), 
                                   weight.change.index = rep(T, n), 
                                   price = price, share = rep(0, n), 
                                   cash = abs(cashflow), commission, 
                                   lot.size, control = control)
    
    if(cashflow >= 0)
      return(list(share = share + out$share, cash = (cash - cashflow) + 
                    out$cash, com = out$com))
    else {
      out = bt.run.share.ex.allocate(weight.new = weight.new, 
                                     weight.prev = rep(0, n), 
                                     weight.change.index = rep(T, n), 
                                     price = price, share = rep(0, n), 
                                     cash = abs(cashflow) + 5 * out$com, 
                                     commission, lot.size, control = control) 
      return(list(
        share = share - out$share, cash = cash + sum(share * price) - 
          (sum((share - out$share) * price) + out$com), com = out$com))
    }
  }
  
  out = bt.run.share.ex.allocate(weight.new = weight.new, 
                                 weight.prev = weight.prev, 
                                 weight.change.index = weight.change.index, 
                                 price = price, share = share, cash = cash, 
                                 commission, lot.size, control = control, 
                                 cashflow = cashflow) 
  
  if(cashflow >= 0) {
    if(any(abs(out$share) < abs(share))) {
      weight.change.index[abs(out$share) < abs(share)] = F
      
      out = bt.run.share.ex.allocate(
        weight.new = weight.new, weight.prev = weight.prev, 
        weight.change.index = weight.change.index, price = price, 
        share = share, cash = cash, commission, lot.size, control = control)
      }
    } else {
      if(any(abs(out$share) > abs(share))) {
        weight.change.index[abs(out$share) > abs(share)] = F
        out = bt.run.share.ex.allocate(
          weight.new = weight.new, weight.prev = weight.prev, 
          weight.change.index = weight.change.index, price = price, 
          share = share, cash = cash, commission, lot.size, 
          control = control, cashflow = cashflow)
      }
    }
  out
  }

bt.run.share.ex.allocate = function(weight.new, weight.prev, 
                                    weight.change.index, price, share, cash, 
                                    commission, lot.size, silent = T, 
                                    control = default.round.lot.control(), 
                                    cashflow = 0) {
  value = sum(price * share) + cash
  compute.commission = function(share.prev, share.new, price, commission) {
    if(is.null(dim(share.new))) {
      share.diff = abs(share.prev - share.new)
      return(
        sum(share.diff) * commission$cps + sum(sign(share.diff)) * 
          commission$fixed + sum(price * share.diff) * commission$percentage)
      }
    
    share.prev = rep.row(share.prev, nrow(share.new))
    price = rep.row(price, nrow(share.new))
    share.diff = abs(share.prev - share.new)
    
    rowSums(share.diff) * commission$cps + rowSums(sign(share.diff)) * 
      commission$fixed + rowSums(price * share.diff) * commission$percentage
    }
  
  compute.cash = function(value, share, price, com) {
    if(is.null(dim(share)))
      value - sum(price * share) - com
    else {
      price = rep.row(price, nrow(share))
      value - rowSums(price * share) - com
    }
  }
  
  compute.weight.diff = function(target, share, cash) {
    if(is.null(dim(share)))
      sum(abs(target - c(share * price, cash) / (sum(price * share) + cash)))
    else {
      price = rep.row(price, nrow(share))
      target = rep.row(target, nrow(share))
      rowSums(abs(target - cbind(share * price, cash) / 
                    (rowSums(price * share) + cash)))
    }
  }
  
  allocate = function(value, share) {
    new.total.weight = sum(abs(weight.new[weight.change.index]))
    if(new.total.weight == 0)
      share[weight.change.index] = 0
    else {
      allocate.value = value * sum(abs(weight.new)) - 
        sum(abs(share * price)[!weight.change.index])
      
      share[weight.change.index] = allocate.value * 
        (weight.new / price)[weight.change.index] / new.total.weight
      }
    share
    }
  
  allocate.lot = function(value, share, lot.size) {
    if(len(lot.size) == 0) return(allocate(value, share))
    new.total.weight = sum(abs(weight.new[weight.change.index]))
    
    if(new.total.weight == 0) {
      shares = rep.row(share, 2)
      shares[2, weight.change.index] = 0
      } else {
        allocate.value = value * sum(abs(weight.new)) - 
          sum(abs(share * price)[!weight.change.index])
        
        lot.size = lot.size[weight.change.index]
        w = weight.new[weight.change.index] / new.total.weight
        p = price[weight.change.index]
        shares = rep.row(share, 3)
        shares[2, weight.change.index] = 
          round.lot.basic(w, p, allocate.value, lot.size)
        
        shares[3, weight.change.index] = 
          round.lot.basic.base(w, p, allocate.value, lot.size)
        }
    shares
    }
  
  new.share = allocate(value, share)
  com = compute.commission(share, new.share, price, commission)
  if(com > 0 || len(lot.size) > 0) {
    share1 = allocate.lot(value - 2 * com, share, lot.size)
    if(cashflow < 0)
      if((value - sum((price * share)[share > 0])) < 0)
        share1 = share1[-1, , drop = F]
      
      com1 = compute.commission(share, share1, price, commission)
      cash1 = compute.cash(value, share1, price, com1)
      target = c(weight.new, 1 - sum(weight.new))
      diff1 = compute.weight.diff(target, share1, cash1)
      j = which.min(diff1)
      
      if( control$select[1] == 'minimum.turnover' ) {
        j1 = which(diff1 - diff1[j] <= control$diff.target)
        j = j1[which.min(com1[j1])]
        }
      
      new.share = share1[j, ]
      com = com1[j]
      }
  
  new.cash = value - sum(price * new.share) - com
  if(!silent) {
    cat('Old[T,V,C]', sum(price * share) + cash, sum(price * share), cash, 
        '\n', 
        'New [T,V,C,COM]', sum(price * new.share) + new.cash + com, 
        sum(price * new.share), new.cash, com, '\n')
    
    cat('Old Weight', weight.new, '\n', 
        'New Weight', price * new.share / (sum(price * new.share) + 
                                             new.cash), '\n')
  }
  list(share = new.share, cash = new.cash, com = com)
}

bt.run.share.ex.allocate.test = function() {
  commission = list(cps = 0.0, fixed = 0.0, percentage = 0.0)
  commission = list(cps = 0.01, fixed = 1.0, percentage = 0.0)
  weight.prev = c(0,0) / 10
  share = c(0, 0)
  price = c(1, 2)
  cash = 100
  lot.size=c()
  weight.new = c(10, 0) / 10
  weight.change.index = c(T, T)
  a = bt.run.share.ex.allocate(weight.new,weight.prev, weight.change.index, 
                               price, share, cash, commission, lot.size, F)
weight.new = c(-10,0) / 10
weight.change.index = c(T, T)
a = bt.run.share.ex.allocate(weight.new,weight.prev,weight.change.index,price,share,cash,commission,lot.size,F)
weight.new = c(13,-3) / 10
weight.change.index = c(T, T)
a = bt.run.share.ex.allocate(weight.new,weight.prev,weight.change.index,price,share,cash,commission,lot.size, F)
weight.new = c(2,8) / 10
weight.change.index = c(T, T)
a = bt.run.share.ex.allocate(weight.new,weight.prev,weight.change.index,price,share,cash,commission,lot.size)
weight.new = c(0,8) / 10
weight.change.index = c(T, T)
a = bt.run.share.ex.allocate(weight.new,weight.prev,weight.change.index,price,share,cash,commission,lot.size)
weight.new = c(0,8) / 10
weight.change.index = c(F, T)
a = bt.run.share.ex.allocate(weight.new,weight.prev,weight.change.index,price,share,cash,commission,lot.size)
weight.new = c(-10,0) / 10
weight.change.index = c(T, T)
a = bt.run.share.ex.allocate(weight.new,weight.prev,weight.change.index,price,share,cash,commission,lot.size)
weight.new = c(-10,0) / 10
weight.change.index = c(T, F)
a = bt.run.share.ex.allocate(weight.new,weight.prev,weight.change.index,price,share,cash,commission,lot.size)
weight.new = c(13,-3) / 10
weight.change.index = c(T, T)
a = bt.run.share.ex.allocate(weight.new,weight.prev,weight.change.index,price,share,cash,commission,lot.size)
weight.new = c(-10,10) / 10
weight.change.index = c(T, T)
a = bt.run.share.ex.allocate(weight.new,weight.prev,weight.change.index,price,share,cash,commission,lot.size)
weight.new = c(10,10) / 10
weight.change.index = c(T, T)
a = bt.run.share.ex.allocate(weight.new,weight.prev,weight.change.index,price,share,cash,commission,lot.size)
weight.new = c(2,8) / 10
weight.prev = c(0,8) / 10
weight.change.index = c(T, F)
price = c(1,2)
share = c(0, 40)
cash = 20
a = bt.run.share.ex.allocate(weight.new,weight.prev,weight.change.index,price,share,cash,commission,lot.size)
weight.new = c(2,8) / 10
weight.prev = c(0,8) / 10
weight.change.index = c(T, T)
price = c(1,2)
share = c(0, 40)
cash = 20
a = bt.run.share.ex.allocate(weight.new,weight.prev,weight.change.index,price,share,cash,commission,lot.size)
}
round.lot = function(weight, price, capital, lot.size) {
weight = coredata(ifna(weight, 0))
price = coredata(ifna(price, 1))
lot.size = ifna(iif( len(lot.size) == 1, rep(lot.size, ncol(weight)), lot.size), 1)
round.lot.basic(weight, price, capital, lot.size)
}
round.lot.basic.base = function(weight, price, capital, lot.size) {
sign(weight) * floor(abs(weight * capital / price / lot.size)) * lot.size
}
round.to = function(x, to) {
sign(x) * floor(abs(x) / to) * to
}
round.lot.basic = function(weight, price, capital, lot.size) {
share = abs(weight * capital) / price
share1 = floor(share / lot.size) * lot.size
discrepancy = (share - share1) * price
cash = sum(discrepancy)
lot.cash = price * lot.size
min.lot.cash = min(lot.cash)
index = (1:len(weight))[cash >= lot.cash & floor(abs(share)) != 0]
for(i in order(discrepancy[index], decreasing=T)) {
if(cash < min.lot.cash) break
j = index[i]
if(cash < lot.cash[j]) next
share1[j] = share1[j] + lot.size[j]
cash = cash - lot.cash[j]
}
sign(weight) * share1
}
round.lot.basic.test = function() {
weight = c(1, 1, 1, 1) / 4
price = c(1.345, 2.4, 3.5, 4.6)
capital = 100
lot.size = c(1, 1, 1, 1)
w = round.lot.basic(weight, price, capital, lot.size)
w
sum(abs(w * price))
weight = c(1, -1, 3, -1) / 4
w = round.lot.basic(weight, price, capital, lot.size)
w
sum(abs(w * price))
w = (weight * capital) / price
w
sum(abs(w * price))
}
bt.run.share.ex.test = function() {
load.packages('quantmod')
tickers = 'SPY'
tickers = 'SPY,XLY,XLP,XLE,XLF,XLV,XLI,XLB,XLK,XLU'
data <- new.env()
getSymbols.extra(tickers, src = 'yahoo', from = '1970-01-01', env = data, auto.assign = T, set.symbolnames=T)
for(i in data$symbolnames) data[[i]] = adjustOHLC(data[[i]], use.Adjusted=T)
bt.prep(data, align='remove.na', fill.gaps = T)
prices = data$prices
n = ncol(prices)
nperiods = nrow(prices)
period.ends = date.ends(data$prices,'months')
models = list()
commission = list(cps = 0.01, fixed = 1.0, percentage = 0.0)
data$weight[] = NA
data$weight[period.ends,] = ntop(prices[period.ends,], n)
models$test = bt.run.share(data, clean.signal=F, silent=F, commission=commission)
data$weight[] = NA
data$weight[period.ends,] = ntop(prices[period.ends,], n)
models$test.ex = bt.run.share.ex(data, clean.signal=F, silent=F, commission=commission)
data$weight[] = NA
data$weight[period.ends,] = ntop(prices[period.ends,], n)
models$test.ex.lot = bt.run.share.ex(data, clean.signal=F, silent=F, commission=commission, lot.size=50)
data$weight[] = NA
data$weight[period.ends,] = ntop(prices[period.ends,], n)
models$test.ex.lot.turnover = bt.run.share.ex(data, clean.signal=F, silent=F, commission=commission,
lot.size=50, control = list(round.lot = list(select = 'minimum.turnover', diff.target = 5/100))
)
last(models$test.ex.lot$share[period.ends,], 20)
strategy.performance.snapshoot(models, T)
layout(1:3)
plotbt.transition.map(models$test$weight, 'BuyHold')
plotbt.transition.map(models$test.ex$weight, 'BuyHold.ex')
plotbt.transition.map(models$test.ex.lot$weight, 'BuyHold.ex')
layout(1:3)
plot(models$test.ex.lot$value, type='l')
plot(models$test.ex.lot$cash[-c(1:20)], type='l')
plot(models$test.ex.lot$com, type='l')
pdf(file = 'report1.pdf', width=8.5, height=11)
strategy.performance.snapshoot(models, data=data)
dev.off()
models = list()
commission = list(cps = 0.01, fixed = 10.0, percentage = 0.0)
obj = portfolio.allocation.helper(data$prices,
period.ends = period.ends, lookback.len = 250, silent=T,
min.risk.fns = list(
EW=equal.weight.portfolio,
RP=risk.parity.portfolio(function(ia) ia$risk),
MV=min.var.portfolio,
Sharpe.RP=risk.parity.portfolio(function(ia) ia$risk / ia$expected.return)
)
)
for(i in names(obj$weights)) {
data$weight[] = NA
data$weight[period.ends,] = obj$weights[[i]]
models[[paste0(i)]] = bt.run.share(data, clean.signal=F, silent=T, commission=commission)
data$weight[] = NA
data$weight[period.ends,] = obj$weights[[i]]
models[[paste0(i,'.ex')]] = bt.run.share.ex(data, clean.signal=F, silent=T, commission=commission)
data$weight[] = NA
data$weight[period.ends,] = obj$weights[[i]]
models[[paste0(i,'.ex.lot')]] = bt.run.share.ex(data, clean.signal=F, silent=T, commission=commission, lot.size=50)
data$weight[] = NA
data$weight[period.ends,] = obj$weights[[i]]
models[[paste0(i,'.ex.lot.turnover')]] = bt.run.share.ex(data, clean.signal=F, silent=T, commission=commission,
lot.size=50, control = list(round.lot = list(select = 'minimum.turnover', diff.target = 5/100)))
}
range(models$MV.ex.lot$share)
strategy.performance.snapshoot(models, T)
pdf(file = 'report2.pdf', width=8.5, height=11)
strategy.performance.snapshoot(models, data=data)
dev.off()
weight = rep(1/n, n)
price = coredata(last(prices))
share = rep(0, n)
cash = 100000
commission = list(cps = 0.01, fixed = 1.0, percentage = 0.0)
lot.size = rep(100, n)
bt.run.share.ex.allocate(weight, weight, rep(T, n),
price, share, cash, commission, lot.size)
}
bt.unadjusted.add.div.split = function(
data.raw,
yahoo.round.up.nearest.cents=F,
infer.div.split.from.adjusted=F
) {
if( !exists('symbolnames', data.raw, inherits = F) )
tickers = ls(data.raw)
else
tickers = data.raw$symbolnames
for(ticker in spl(tickers)) {
price = data.raw[[ticker]]
price$Dividend = price$Split = 0
if(infer.div.split.from.adjusted) {
close = coredata(Cl(price))
adjusted = coredata(Ad(price))
implied.split = close / mlag(close) * mlag(adjusted) / adjusted
isplit.index = ifna(implied.split < 0.9 | implied.split > 1.2,F)
isplit = implied.split[isplit.index]
isplit = round(100 * isplit) / 100
implied.div = mlag(close) * adjusted / mlag(adjusted) - close
idiv.index = ifna(implied.div > 1e-3, F) & !isplit.index
idiv = implied.div[idiv.index]
idiv = round(1e3 * idiv) / 1e3
price$Dividend[idiv.index] = idiv
price$Split[isplit.index] = isplit
} else {
dividend = getDividends(ticker, from = '1900-01-01')
split = getSplits(ticker, from = '1900-01-01')
split = split[split > 0]
dividend = dividend[dividend > 0]
if(is.xts(split) && is.xts(dividend) && nrow(split) > 0 && nrow(dividend) > 0)
dividend = dividend * 1/adjRatios(splits=merge(split, index(dividend)))[,1]
dividend = dividend[index(price)]
split = split[index(price)]
if( is.xts(dividend) && nrow(dividend) > 0 )
if( nrow(price[index(dividend)]) != nrow(dividend) )
stop(paste('Missing Price date for dividend. Symbol =', ticker))
else
price$Dividend[index(dividend)] = dividend
if( is.xts(split) && nrow(split) > 0 )
if( nrow(price[index(split)]) != nrow(split) )
stop(paste('Missing Price date for split. Symbol =', ticker))
else
price$Split[index(split)] = split
}
if(yahoo.round.up.nearest.cents) {
map.col = unlist(find.names('Close,Open,High,Low,Adjusted', price, F))
price[,map.col] = ceiling(100 * price[,map.col]) / 100
}
data.raw[[ticker]] = price
}
}
test.implied.div.split = function(ticker) {
load.packages('quantmod')
ticker = 'IBM'
data <- new.env()
getSymbols.extra(ticker, src = 'yahoo', from = '1970-01-01', env = data, auto.assign = T, set.symbolnames=T)
dividend = getDividends(ticker, from = '1900-01-01')
split = getSplits(ticker, from = '1900-01-01')
if(is.xts(split) && is.xts(dividend) && nrow(split) > 0 && nrow(dividend) > 0)
dividend1 = dividend * 1/adjRatios(splits=merge(split, index(dividend)))[,1]
else
dividend1 = dividend
close = Cl(data$IBM)
adjusted = Ad(data$IBM)
implied.split = close / mlag(close) * mlag(adjusted) / adjusted
isplit.index = implied.split < 0.8 | implied.split > 1.2
isplit = implied.split[isplit.index]
isplit = round(100 * isplit) / 100
cbind(isplit['1970::'], split['1970::'])
implied.div = mlag(close) * adjusted / mlag(adjusted) - close
idiv.index = implied.div > 1e-3
idiv = implied.div[idiv.index & !isplit.index]
idiv = round(1e3 * idiv) / 1e3
len(idiv['1970::'])
len(dividend1['1970::'])
setdiff( index(dividend1['1970::']), index(idiv['1970::']))
setdiff( index(idiv['1970::']), index(dividend1['1970::']) )
cbind(idiv['1970::'], dividend1['1970::'])
tickers = dow.jones.components()
for(ticker in tickers) {
data <- new.env()
getSymbols.extra(ticker, src = 'yahoo', from = '1970-01-01', env = data, auto.assign = T)
dividend = getDividends(ticker, from = '1900-01-01')
split = getSplits(ticker, from = '1900-01-01')
split = split[split > 0]
dividend = dividend[dividend > 0]
if(is.xts(split) && is.xts(dividend) && nrow(split) > 0 && nrow(dividend) > 0)
dividend1 = dividend * 1/adjRatios(splits=merge(split, index(dividend)))[,1]
else
dividend1 = dividend
close = Cl(data[[ticker]])
adjusted = Ad(data[[ticker]])
implied.split = close / mlag(close) * mlag(adjusted) / adjusted
isplit.index = ifna(implied.split < 0.9 | implied.split > 1.1, F)
isplit = implied.split[isplit.index]
isplit = round(100 * isplit) / 100
if(len(isplit)>0)
cat(ticker, 'SPL', len(isplit['1970::']) - len(split['1970::']), max(round(isplit['1970::'],3) - round(split['1970::'],3)), '\n')
else
cat(ticker, 'SPL', len(isplit['1970::']) - len(split['1970::']), '\n')
implied.div = mlag(close) * adjusted / mlag(adjusted) - close
idiv.index = ifna(implied.div > 1e-3, F)
idiv = implied.div[idiv.index & !isplit.index]
idiv = round(1e3 * idiv) / 1e3
len(idiv['1970::'])
len(dividend1['1970::'])
cat(ticker, 'DIV', len(idiv['1970::']) - len(dividend1['1970::']),
len(setdiff( index(dividend1['1970::']), index(idiv['1970::']))),
len(setdiff( index(idiv['1970::']), index(dividend1['1970::']) )),
max(round(idiv['1970::'],3)- round(dividend1['1970::'],3)), '\n')
setdiff( index(dividend1['1970::']), index(idiv['1970::']))
setdiff( index(idiv['1970::']), index(dividend1['1970::']) )
}
}
bt.run.share.unadjusted.test.data = function() {
load.packages('quantmod')
tickers = 'IBM'
data = env()
getSymbols.extra(tickers, src = 'yahoo', from = '1970-01-01', env = data, auto.assign = T, set.symbolnames=T)
data.raw = env(data)
bt.unadjusted.add.div.split(data.raw, yahoo.round.up.nearest.cents = T)
ticker = 'IBM'
adjusted = data.raw[[ticker]]$Adjusted / mlag(data.raw[[ticker]]$Adjusted) - 1
adjusted = ifna(adjusted,0)
prod(1 + adjusted)
split = iif(data.raw[[ticker]]$Split > 0, 1 / data.raw[[ticker]]$Split, 1)
unadjusted = (data.raw[[ticker]]$Close * split +  data.raw[[ticker]]$Dividend) / mlag(data.raw[[ticker]]$Close) - 1
unadjusted = ifna(unadjusted,0)
prod(1 + unadjusted)
index = which(round(adjusted - unadjusted,4) != 0)
cbind(round(adjusted - unadjusted, 4), data.raw[[ticker]]$Split, data.raw[[ticker]]$Dividend)[index]
plota.matplot(cbind(cumprod(1 + adjusted),cumprod(1 + unadjusted)))
index.max = which.max(abs(adjusted - unadjusted))
cbind(adjusted, unadjusted, round(adjusted - unadjusted, 4), data.raw[[ticker]]$Split, data.raw[[ticker]]$Dividend)[index.max]
data.raw[[ticker]][(index.max-1):index.max,]
(65.875 + 0.3025) / 68.250 - 1
12.04279 / 12.25577 - 1
16.47 /17.06 -1
8.4485 / 7.8843 - 1
yhist = read.xts(hist.quotes.url('IBM', '1992-11-01', '1992-11-30', 'yahoo'))
ghist = read.xts(hist.quotes.url('IBM', '1992-11-01', '1992-11-30', 'google'), format='%d-%b-%y')
qhist = read.xts(hist.quotes.url('IBM', '1992-11-01', '1992-11-30', 'quotemedia'))
dividends = getDividends('IBM', from = '1900-01-01')
splits = getSplits('IBM', from = '1900-01-01')
dividends['1992:11::1992:11']
if(is.xts(splits) && is.xts(dividends) && nrow(splits) > 0 && nrow(dividends) > 0)
dividends = dividends * 1/adjRatios(splits=merge(splits, index(dividends)))[,1]
dividends['1992:11::1992:11']
dividend = dividends[index(yhist)]
yhist = yhist[,spl('Close,Adj_Close')]
colnames(yhist) = spl('Yahoo.Close,Yahoo.Adjusted')
yhist$Dividend = 0
yhist$Dividend[index(dividend)] = dividend
yhist
ghist = ghist[,'Close']
colnames(ghist) = spl('Google.Adjusted')
qhist = qhist[,c('close', 'adjclose')]
colnames(qhist) = spl('Quotemedia.Close,Quotemedia.Adjusted')
temp = cbind(yhist, ghist, qhist)
temp[,spl('Yahoo.Close,Dividend,Quotemedia.Close')]
to.return = function(x) round(100*(x/mlag(x)-1),3)
Yahoo.Return = to.return(temp$Yahoo.Close + temp$Dividend)
Yahoo.Return1 = to.return(ceiling(100*temp$Yahoo.Close)/100 + temp$Dividend)
Yahoo.Return.Adjusted = to.return(temp$Yahoo.Adjusted)
Google.Return.Adjusted = to.return(temp$Google.Adjusted)
Quotemedia.Return = to.return(temp$Quotemedia.Close + temp$Dividend)
Quotemedia.Return.Adjusted = to.return(temp$Quotemedia.Adjusted)
ret = cbind(Yahoo.Return, Yahoo.Return1, Yahoo.Return.Adjusted, Google.Return.Adjusted, Quotemedia.Return, Quotemedia.Return.Adjusted)
t(apply(ret,1,range))
t(diff(apply(ret,1,range)))
ret['1992:11:05']
txt = '
Date 		Open High Low Close Volume
Nov-2-1992 67.00 69.00 67.00 68.88 2,322,100
Nov-3-1992 68.50 69.88 68.50 69.13 2,375,200
Nov-4-1992 69.00 69.63 68.13 68.25 2,079,800
Nov-5-1992 67.13 67.25 65.63 65.88 2,136,200
Nov-6-1992 65.38 66.50 65.00 66.25 2,642,300
Nov-9-1992 66.25 67.63 66.25 67.50 2,216,400
Nov-10-1992 67.63 68.00 65.88 65.88 2,187,100
Nov-11-1992 65.63 65.75 64.50 65.00 3,145,100
Nov-12-1992 65.13 65.38 64.13 64.13 3,133,000
Nov-13-1992 64.88 65.13 64.00 64.88 1,851,300
Nov-16-1992 64.75 65.50 64.63 64.88 1,765,100
Nov-17-1992 64.88 65.00 64.00 64.25 2,020,700
Nov-18-1992 64.13 64.38 62.75 63.13 2,707,100
Nov-19-1992 63.00 63.13 61.00 61.25 3,307,600
Nov-20-1992 61.38 62.63 60.88 62.25 3,715,200
Nov-23-1992 62.38 63.88 62.25 63.25 2,220,200
Nov-24-1992 63.75 65.50 63.50 64.88 2,847,100
Nov-25-1992 65.38 66.00 65.13 65.38 1,788,700
Nov-27-1992 65.88 66.25 65.25 66.00 1,229,500
Nov-30-1992 67.88 68.63 67.50 68.25 3,239,000
'
IBM = read.xts(txt,sep=' ', format='%b-%d-%Y')
IBM = IBM[,'Close']
colnames(IBM) = spl('IBM.Close')
IBM$Dividend = 0
IBM$Dividend['1992:11:05'] = 1.21
IBM.Return = to.return(IBM$IBM.Close + IBM$Dividend)
ret = cbind(IBM.Return, Yahoo.Return, Yahoo.Return1, Yahoo.Return.Adjusted, Google.Return.Adjusted, Quotemedia.Return, Quotemedia.Return.Adjusted)
ret['1992:11:05']
ret$IBM.Close - ret$Yahoo.Close.1
setdiff(index, which(data.raw[[ticker]]$Dividend > 0))
temp.adjusted = adjusted
temp.adjusted[index] = 0
prod(1 + temp.adjusted)
temp.unadjusted = unadjusted
temp.unadjusted[index] = 0
prod(1 + temp.unadjusted)
plota.matplot(cbind(cumprod(1 + temp.adjusted),cumprod(1 + temp.unadjusted)))
}
bt.make.trade.event.summary.table = function(bt, to.text=F) {
index = bt$event.type != 'none'
out = data.frame(
Type = bt$event.type,
bt$share,
Cash=bt$cash,
Com=bt$com,
Div=bt$div,
Value=bt$value
)[index,]
rownames(out) = format(index(bt$equity)[index], '%Y-%m-%d')
if(to.text) to.nice(out,0)
else out
}
bt.make.cashflow.event.summary.table = function(bt, to.text=F) {
if( is.null(bt$cashflow) ) return()
index = rowSums(bt$cashflows != 0) > 0
out = data.frame(
Cashflow = bt$cashflow,
Fee.Rebate = bt$fee.rebate,
bt$cashflows
)[index,]
rownames(out) = format(index(bt$equity)[index], '%Y-%m-%d')
if(to.text) to.nice(out,0)
else out
}
bt.run.share.unadjusted.test = function() {
load.packages('quantmod')
tickers = 'IBM'
data = env()
getSymbols.extra(tickers, src = 'yahoo', from = '1970-01-01', env = data, auto.assign = T, set.symbolnames=T)
data.raw = env(data)
data.raw1 = env(data)
for(i in data$symbolnames) data[[i]] = adjustOHLC(data[[i]], use.Adjusted=T)
bt.prep(data, align='remove.na', fill.gaps = T)
prices = data$prices
n = ncol(prices)
nperiods = nrow(prices)
period.ends = date.ends(data$prices,'months')
models = list()
commission = list(cps = 0.01, fixed = 1.0, percentage = 0.0)
data$weight[] = NA
data$weight[period.ends,] = 1
models$test = bt.run.share(data, clean.signal=F, silent=F, commission=commission)
data$weight[] = NA
data$weight[period.ends,] = 1
models$test.ex = bt.run.share.ex(data, clean.signal=F, silent=F, commission=commission,
lot.size=1)
bt.unadjusted.add.div.split(data.raw, yahoo.round.up.nearest.cents = T)
bt.prep(data.raw, align='remove.na', fill.gaps = T)
prices = data.raw$prices
n = ncol(prices)
nperiods = nrow(prices)
period.ends = date.ends(data.raw$prices,'months')
commission = list(cps = 0.01, fixed = 1.0, percentage = 0.0)
data.raw$weight[] = NA
data.raw$weight[period.ends,] = 1
models$test.unadjusted = bt.run.share.ex(data.raw, clean.signal=F, silent=F, commission=commission,
lot.size=50, adjusted = F)
data.raw = data.raw1
bt.unadjusted.add.div.split(data.raw, yahoo.round.up.nearest.cents = T, infer.div.split.from.adjusted=T)
bt.prep(data.raw, align='remove.na', fill.gaps = T)
prices = data.raw$prices
n = ncol(prices)
nperiods = nrow(prices)
period.ends = date.ends(data.raw$prices,'months')
commission = list(cps = 0.01, fixed = 1.0, percentage = 0.0)
data.raw$weight[] = NA
data.raw$weight[period.ends,] = 1
models$test.unadjusted1 = bt.run.share.ex(data.raw, clean.signal=F, silent=F, commission=commission,
lot.size=50, adjusted = F)
strategy.performance.snapshoot(models, T)
layout(1:2)
plotbt.transition.map(models$test.ex$weight)
plotbt.transition.map(models$test.unadjusted$weight)
mlast(bt.make.trade.event.summary.table(models$test.unadjusted), 20)
load.packages('quantmod')
tickers = 'SPY'
tickers = 'SPY,XLY,XLP,XLE,XLF,XLV,XLI,XLB,XLK,XLU'
tickers = 'AAPL,IBM,VTI,IEV,EWJ,EEM,RWX,DBC,GLD,TLT,IEF,SHY'
data = env()
getSymbols.extra(tickers, src = 'yahoo', from = '1970-01-01', env = data, auto.assign = T, set.symbolnames=T)
data.raw = env(data)
for(i in data$symbolnames) data[[i]] = adjustOHLC(data[[i]], use.Adjusted=T)
bt.prep(data, align='remove.na', fill.gaps = T)
prices = data$prices
n = ncol(prices)
nperiods = nrow(prices)
period.ends = date.ends(data$prices,'months')
models = list()
commission = list(cps = 0.01, fixed = 1.0, percentage = 0.0)
data$weight[] = NA
data$weight[period.ends,] = ntop(prices[period.ends,], n)
models$test = bt.run.share(data, clean.signal=F, silent=F, commission=commission)
data$weight[] = NA
data$weight[period.ends,] = ntop(prices[period.ends,], n)
models$test.ex = bt.run.share.ex(data, clean.signal=F, silent=F, commission=commission,
lot.size=50)
strategy.performance.snapshoot(models, T)
plotbt.transition.map(models$test.ex$weight)
bt.unadjusted.add.div.split(data.raw)
bt.prep(data.raw, align='remove.na', fill.gaps = T)
prices = data.raw$prices
n = ncol(prices)
nperiods = nrow(prices)
period.ends = date.ends(data.raw$prices,'months')
commission = list(cps = 0.01, fixed = 1.0, percentage = 0.0)
data.raw$weight[] = NA
data.raw$weight[period.ends,] = ntop(prices[period.ends,], n)
models$test.unadjusted = bt.run.share.ex(data.raw, clean.signal=F, silent=F, commission=commission,
lot.size=50, adjusted = F)
data.raw$weight[] = NA
data.raw$weight[period.ends,] = ntop(prices[period.ends,], n)
models$test.unadjusted1 = bt.run.share.ex(data.raw, clean.signal=F, silent=F, commission=commission,
lot.size=50,
adjusted = F,
control = list(round.lot = list(select = 'minimum.turnover', diff.target = 5/100)),
dividend.control = list(foreign.withholding.tax = 30/100)
)
strategy.performance.snapshoot(models, T)
layout(1:2)
plotbt.transition.map(models$test.unadjusted$weight)
plotbt.transition.map(models$test.unadjusted1$weight)
mlast(bt.make.trade.event.summary.table(models$test.unadjusted), 20)
pdf(file = 'report.u.pdf', width=8.5, height=11)
strategy.performance.snapshoot(models, data=data)
dev.off()
}
bt.run.share.ex.test.cashflow = function() {
load.packages('quantmod')
tickers = 'SPY'
tickers = 'SPY,XLY,XLP,XLE,XLF,XLV,XLI,XLB,XLK,XLU'
data <- new.env()
getSymbols.extra(tickers, src = 'yahoo', from = '1970-01-01', env = data, auto.assign = T, set.symbolnames=T)
for(i in data$symbolnames) data[[i]] = adjustOHLC(data[[i]], use.Adjusted=T)
bt.prep(data, align='remove.na', fill.gaps = T)
prices = data$prices
n = ncol(prices)
nperiods = nrow(prices)
period.ends = date.ends(data$prices,'months')
models = list()
commission = list(cps = 0.01, fixed = 1.0, percentage = 0.0)
data$weight[] = NA
data$weight[period.ends,] = ntop(prices[period.ends,], n)
models$test = bt.run.share(data, clean.signal=F, silent=F, commission=commission)
data$weight[] = NA
data$weight[period.ends,] = ntop(prices[period.ends,], n)
models$test.ex = bt.run.share.ex(data, clean.signal=F, silent=F, commission=commission)
data$weight[] = NA
data$weight[period.ends,] = ntop(prices[period.ends,], n)
models$test.ex.lot = bt.run.share.ex(data, clean.signal=F, silent=F, commission=commission, lot.size=50)
data$weight[] = NA
data$weight[period.ends,] = ntop(prices[period.ends,], n)
models$test.ex.lot.cashflow = bt.run.share.ex(data, clean.signal=F, silent=F, commission=commission,
lot.size=50,
cashflow.control = list(
monthly.income = list(
cashflows = event.at(prices, 'quarter', 1000, offset=0),
invest = 'cash',
type = 'regular'
)
)
)
mlast(bt.make.trade.event.summary.table(models$test.ex.lot.cashflow), 20)
mlast(bt.make.cashflow.event.summary.table(models$test.ex.lot.cashflow), 20)
matplot(cbind(
models$test.ex.lot$value,
models$test.ex.lot.cashflow$value
), type='l')
pdf(file = 'report.c.pdf', width=8.5, height=11)
strategy.performance.snapshoot(models, data=data)
dev.off()
}
bt.run.share.ex.invest.test = function() {
commission = list(cps = 0.01, fixed = 1.0, percentage = 0.0)
weight.prev = c(10,0) / 10
share = c(100, 0)
price = c(1,2)
cash = 0
lot.size=c()
weight.new = c(10,0) / 10
weight.change.index = c(T, T)
cashflow = -10
cash = cash + cashflow
a = bt.run.share.ex.invest(weight.new,weight.prev,weight.change.index,price,share,cash,cashflow,commission,lot.size,F)
a
}
bt.run.share.ex.test.tax = function() {
load.packages('quantmod')
tickers = 'SPY'
tickers = 'SPY,XLY,XLP,XLE,XLF,XLV,XLI,XLB,XLK,XLU'
tickers = 'AAPL,IBM,VTI,IEV,EWJ,EEM,RWX,DBC,GLD,TLT,IEF,SHY'
data <- new.env()
getSymbols.extra(tickers, src = 'yahoo', from = '1970-01-01', env = data, auto.assign = T, set.symbolnames=T)
data.raw = env(data)
for(i in data$symbolnames) data[[i]] = adjustOHLC(data[[i]], use.Adjusted=T)
bt.prep(data, align='remove.na', fill.gaps = T)
prices = data$prices
n = ncol(prices)
nperiods = nrow(prices)
period.ends = date.ends(data$prices,'months')
models = list()
commission = list(cps = 0.01, fixed = 1.0, percentage = 0.0)
weights = ntop(prices[period.ends,], n)
data$weight[] = NA
data$weight[period.ends,] = weights
models$test = bt.run.share(data, clean.signal=F, silent=F, commission=commission)
data$weight[] = NA
data$weight[period.ends,] = weights
models$test.ex = bt.run.share.ex(data, clean.signal=F, silent=F, commission=commission)
data$weight[] = NA
data$weight[period.ends,] = weights
models$test.ex.lot = bt.run.share.ex(data, clean.signal=F, silent=F, commission=commission, lot.size=50)
data$weight[] = NA
data$weight[period.ends,] = weights
models$test.ex.lot.tax = bt.run.share.ex(data, clean.signal=F, silent=F, commission=commission,
lot.size=50,
control = list(round.lot = list(select = 'minimum.turnover', diff.target = 5/100)),
tax.control = default.tax.control(),
cashflow.control = list(
taxes = list(
cashflows = event.at(prices, 'year', offset=60),
cashflow.fn = tax.cashflows,
invest = 'cash',
type = 'fee.rebate'
)
)
)
mlast(bt.make.trade.event.summary.table(models$test.ex.lot.tax), 20)
mlast(bt.make.cashflow.event.summary.table(models$test.ex.lot.tax), 20)
pdf(file = 'report.t.pdf', width=8.5, height=11)
strategy.performance.snapshoot(models, data=data)
dev.off()
bt.unadjusted.add.div.split(data.raw)
bt.prep(data.raw, align='remove.na', fill.gaps = T)
data.raw$weight[] = NA
data.raw$weight[period.ends,] = weights
models$test.unadjusted = bt.run.share.ex(data.raw, clean.signal=F, silent=F, commission=commission,
lot.size=50, adjusted = F)
data.raw$weight[] = NA
data.raw$weight[period.ends,] = weights
models$test.unadjusted1 = bt.run.share.ex(data.raw, clean.signal=F, silent=F, commission=commission,
lot.size=50,
adjusted = F,
control = list(round.lot = list(select = 'minimum.turnover', diff.target = 5/100)),
dividend.control = list(foreign.withholding.tax = 30/100)
)
data$weight[] = NA
data$weight[period.ends,] = weights
models$test.unadjusted.tax = bt.run.share.ex(data.raw, clean.signal=F, silent=F, commission=commission,
lot.size=50,
control = list(round.lot = list(select = 'minimum.turnover', diff.target = 5/100)),
adjusted = F,
tax.control = default.tax.control(),
cashflow.control = list(
taxes = list(
cashflows = event.at(prices, 'year', offset=60),
cashflow.fn = tax.cashflows,
invest = 'cash',
type = 'fee.rebate'
)
)
)
mlast(bt.make.trade.event.summary.table(models$test.unadjusted.tax), 20)
mlast(bt.make.cashflow.event.summary.table(models$test.unadjusted.tax), 20)
models$test1 = models$test.ex.lot.tax
models$test2 = models$test.unadjusted.tax
look.at.taxes(models$test1)['2007']
look.at.taxes(models$test2)['2007']
tax.summary(models$test1)
tax.summary(models$test2)
tax.summary(models$test1, function(x) as.numeric(format(x,'%Y%m')))[1:14,]
tax.summary(models$test2, function(x) as.numeric(format(x,'%Y%m')))[1:14,]
data$prices[period.ends,][1:14,1:3]
data.raw$prices[period.ends,][1:14,1:3]
bt.make.trade.event.summary.table(models$test1)[1:14,]
bt.make.trade.event.summary.table(models$test2)[1:36,]
if(F) {
models$test.unadjusted.tax$long.term.cap
models$test.unadjusted.tax$short.term.cap
models$test.unadjusted.tax$qualified.div
models$test.unadjusted.tax$non.qualified.div
}
pdf(file = 'report.t.u.pdf', width=8.5, height=11)
strategy.performance.snapshoot(models, data=data)
dev.off()
}
event.at = function(x, period = 'month', amount = 1, period.ends = date.ends(x,period), offset = 1) {
nperiods = nrow(x)
index = period.ends + offset
index[index > nperiods] = nperiods
index[index < 1] = 1
cashflow = x[index,1]
cashflow[] = amount
cashflow
}
look.at.taxes = function(m) {
temp = data.frame(m[spl('long.term.cap,short.term.cap,qualified.div,non.qualified.div')])
temp = make.xts(cbind(temp, total=rowSums(temp)), data$dates)
temp[temp$total!=0]
}
tax.summary = function(m, by.fn=date.year) {
temp = aggregate(
m[spl('long.term.cap,short.term.cap,qualified.div,non.qualified.div')],
list(year=by.fn(data$dates)),
sum
)
cbind(temp, total=rowSums(temp))
}
tax.cashflows = function(info, index, last.index) {
ii = date.year(info$dates) == date.year(info$dates[index]) - 1
if(sum(ii) == 0) return(0)
if( is.null(info$tax.copy) ) {
info$tax.copy = env(
long.term.cap = info$tax$long.term.cap,
short.term.cap = info$tax$short.term.cap
)
} else {
info$tax.copy$long.term.cap[ii] = info$tax$long.term.cap[ii]
info$tax.copy$short.term.cap[ii] = info$tax$short.term.cap[ii]
}
i = max(which(ii))
long.term.cap = sum(info$tax.copy$long.term.cap[1:i])
short.term.cap = sum(info$tax.copy$short.term.cap[1:i])
tax.cashflows.cap.helper = function(neg, pos) {
if( -neg > pos) {
neg = neg + pos
pos = 0
} else {
pos = pos + neg
neg = 0
}
return(list(neg = neg, pos = pos))
}
if(long.term.cap < 0 && short.term.cap > 0) {
temp = tax.cashflows.cap.helper(long.term.cap, short.term.cap)
long.term.cap = temp$neg
short.term.cap = temp$pos
} else if(long.term.cap > 0 && short.term.cap < 0) {
temp = tax.cashflows.cap.helper(short.term.cap, long.term.cap)
long.term.cap = temp$pos
short.term.cap = temp$neg
}
tax = 0
info$tax.copy$long.term.cap[1:i] = 0
if(long.term.cap >= 0)
tax = tax + long.term.cap * info$tax.control$capital.gain$long.term.tax
else
info$tax.copy$long.term.cap[i] = long.term.cap
info$tax.copy$short.term.cap[1:i] = 0
if(short.term.cap >= 0)
tax = tax + short.term.cap * info$tax.control$capital.gain$short.term.tax
else
info$tax.copy$short.term.cap[i] = short.term.cap
qualified.div = sum(info$tax$qualified.div[ii])
non.qualified.div = sum(info$tax$non.qualified.div[ii])
tax = tax + qualified.div * info$tax.control$dividend$qualified.tax
tax = tax + non.qualified.div * info$tax.control$dividend$nonqualified.tax
-tax
}
bt.run.share.fast <- function
(
b,
clean.signal = T,
do.lag = 1,
capital = 100000,
lot.size = c()
)
{
weight = b$weight
weight[is.nan(weight) | is.infinite(weight)] = NA
weight[!is.na(weight) & is.na(b$prices)] = 0
weight = iif( do.lag == 1, weight, mlag(weight, do.lag - 1) )
weight = coredata(weight)
if(F) {
temp = bt.exrem(weight)
if(clean.signal) {
weight = temp
} else {
index = ifna(weight == 0, F)
weight[index] = temp[index]
}
}
prices = coredata(prices)
n = ncol(prices)
nperiods = nrow(prices)
trade = !is.na(weight)
trade.index = which(rowSums(trade) > 0)
cash.wt = cash = rep(capital, nperiods)
share.wt = share = matrix(0, nperiods, n)
last.trade = 0
lot.size = map2vector(lot.size, colnames(prices), 1)
lot.size = rep(1,n)
for(i in trade.index) {
if(last.trade > 0) {
index = (last.trade + 1) : i
n.index = len(index)
share.wt[index,] = rep.row(share[last.trade,], n.index)
cash.wt[index] = cash[last.trade]
share[index,] = rep.row(share[last.trade,], n.index)
cash[index] = cash[last.trade]
}
p = prices[i,]
p[is.na(p)] = 1
w = weight[i,]
w[is.na(w)] = 0
value = cash[i] + sum(p * share[i,])
share[i,] = round.lot.basic(w, p, value, lot.size)
cash[i] = value - sum(share[i,] * p)
last.trade = i
}
if( last.trade > 0 & last.trade < nperiods) {
index = (last.trade + 1) : nperiods
n.index = len(index)
share.wt[index,] = rep.row(share[last.trade,], n.index)
cash.wt[index] = cash[last.trade]
share[index,] = rep.row(share[last.trade,], n.index)
cash[index] = cash[last.trade]
}
bt = list(type = 'share', capital = capital, share=share)
bt$value = cash + rowSums(share * prices, na.rm=T)
bt$weight = share.wt * prices / (cash.wt + rowSums(share.wt * prices, na.rm=T))
value = c(capital, bt$value)
bt$ret = (value / mlag(value) - 1)[-1]
bt$equity = cumprod(1 + bt$ret)
bt
}
bt.exrem.time.exit <- function(signal, nlen, create.weight = T) {
signal[is.na(signal)] = FALSE
signal.index = which(signal)
nsignal.index = len(signal.index)
nperiods = len(signal)
signal.index.exit = iif(signal.index + nlen - 1 > nperiods, nperiods, signal.index + nlen)
if(!create.weight) {
for(i in 1:nsignal.index) {
if( signal[ signal.index[i] ] ) {
signal[ (signal.index[i]+1) : signal.index.exit[i] ] = FALSE
}
}
return(signal)
} else {
temp = signal * NA
for(i in 1:nsignal.index) {
if( signal[ signal.index[i] ] ) {
signal[ (signal.index[i]+1) : signal.index.exit[i] ] = FALSE
temp[ signal.index.exit[i] ] = 0
}
}
temp[signal] = 1
return(temp)
}
}
bt.min.holding.period <- function(x, nlen) {
x = coredata(x)
enter = x != 0
enter[is.na(enter)] = FALSE
enter.index = which(enter)
for(t in enter.index)
if( enter[ t ] ) {
index = t + nlen
enter[ t : index ] = FALSE
x[ t : index ] = x[t]
}
return(x)
}
bt.time.stop <- function(weight, nlen)
{
bt.apply.matrix(weight, bt.ts.time.stop, nlen)
}
bt.price.stop <- function(b, price, pstop)
{
out = b
out[] = NA
nsymbols = ncol(b)
if(is.null(dim(pstop))) pstop = rep.row(pstop, nrow(b))
for( i in 1:nsymbols )
out[,i] = bt.ts.price.stop(coredata(b[,i]), coredata(price[,i]), coredata(pstop[,i]))
return(out)
}
bt.time.price.stop <- function(b, nlen, price, pstop)
{
out = b
out[] = NA
nsymbols = ncol(b)
if(is.null(dim(pstop))) pstop = rep.row(pstop, nrow(b))
for( i in 1:nsymbols )
out[,i] = bt.ts.time.price.stop(coredata(b[,i]), nlen, coredata(price[,i]), coredata(pstop[,i]))
return(out)
}
bt.ts.trade.index <- function(x)
{
enter = x != 0
enter[is.na(enter)] = FALSE
enter[length(x)] = FALSE
enter.index = which(enter)
temp = ifna.prev(x)
temp0 = mlag(temp)
exit = temp0 != 0 & temp != temp0
exit[ !exit ] = NA
exit = ifna.prevx.rev(exit)
list(enter = enter, enter.index = enter.index, exit = exit)
}
bt.ts.enter.state <- function(x)
{
enter = x != 0
enter[is.na(enter)] = FALSE
enter
}
bt.ts.time.stop <- function(x, nlen)
{
temp = bt.ts.trade.index(x)
enter = temp$enter
enter.index = temp$enter.index
exit = temp$exit
for(t in enter.index)
if( enter[ t ] )
if( exit[ t ] < t + nlen )
enter[ t : exit[ t ] ] = FALSE
else {
enter[ t : (t + nlen) ] = FALSE
x[ (t + nlen) ] = 0
}
return(x)
}
time.stop.test <- function() {
bt.ts.time.stop(c(1,1,1,0,1,1,NA,1,0),2)
bt.ts.time.stop(c(1,0,1,1,1,1,1,1,1),3)
}
bt.ts.price.stop <- function(x, price, pstop)
{
price = coredata(price)
pstop = coredata(pstop)
if(length(pstop) == 1) pstop = rep(pstop, len(x))
dummy = 1:length(x)
temp = bt.ts.trade.index(x)
enter = temp$enter
enter.index = temp$enter.index
exit = temp$exit
for(t in enter.index)
if( enter[ t ] ) {
if( x[ t ] > 0 )
temp = price[ t : exit[ t ] ] < price[ t ] - pstop[ t ]
else
temp = price[ t : exit[ t ] ] > price[ t ] + pstop[ t ]
if( any(temp, na.rm=T) ) {
iexit = t - 1 + dummy[temp][1]
enter[ t : iexit ] = FALSE
x[ iexit ] = 0
} else
enter[ t : exit[ t ] ] = FALSE
}
return(x)
}
price.stop.test <- function() {
bt.ts.price.stop(c(1,1,1,1,1,1,NA,1,0),
c(1,1,0.9,0.7,1,1,1,1,0),
0.2
)
bt.ts.price.stop(-c(1,1,1,1,1,1,NA,1,0),
c(1,1,0.9,1.7,1,1,1,1,0),
0.2
)
}
bt.ts.time.price.stop <- function(x, nlen, price, pstop)
{
price = coredata(price)
pstop = coredata(pstop)
if(length(pstop) == 1) pstop = rep(pstop, len(x))
dummy = 1:length(x)
temp = bt.ts.trade.index(x)
enter = temp$enter
enter.index = temp$enter.index
exit = temp$exit
for(t in enter.index)
if( enter[ t ] ) {
if( x[ t ] > 0 )
temp = price[ t : exit[ t ] ] < price[ t ] - pstop[ t ]
else
temp = price[ t : exit[ t ] ] > price[ t ] + pstop[ t ]
if( any(temp, na.rm=T) ) {
iexit = t - 1 + dummy[temp][1]
if( iexit < t + nlen ) {
enter[ t : iexit ] = FALSE
x[ iexit ] = 0
} else {
enter[ t : (t + nlen) ] = FALSE
x[ (t + nlen) ] = 0
}
} else
if( exit[ t ] < t + nlen )
enter[ t : exit[ t ] ] = FALSE
else {
enter[ t : (t + nlen) ] = FALSE
x[ (t + nlen) ] = 0
}
}
return(x)
}
time.price.stop.test <- function() {
bt.ts.time.price.stop(c(1,1,1,1,1,1,NA,1,0),
4,
c(1,1,0.9,0.7,1,1,1,1,0),
0.2
)
bt.ts.time.price.stop(-c(1,1,1,1,1,1,NA,1,0),
4,
c(1,1,0.9,1.7,1,1,1,1,0),
0.2
)
}
custom.stop.fn <- function(x, price, stop.fn, ...)
{
price = coredata(price)
if(is.character(stop.fn)) stop.fn = match.fun(stop.fn)
dummy = 1:length(x)
temp = bt.ts.trade.index(x)
enter = temp$enter
enter.index = temp$enter.index
exit = temp$exit
for(t in enter.index)
if( enter[ t ] ) {
temp = stop.fn(x[ t ], price, t, exit[ (t + 1) ], ...)
if( any(temp, na.rm=T) ) {
iexit = t - 1 + dummy[temp][1]
enter[ t : iexit ] = FALSE
x[ iexit ] = 0
} else
enter[ t : exit[ t ] ] = FALSE
}
return(x)
}
custom.stop.fn.list <- function(x, price, stop.fn, ...)
{
price = coredata(price)
if(is.character(stop.fn)) stop.fn = match.fun(stop.fn)
dummy = 1:length(x)
temp = bt.ts.trade.index(x)
enter = temp$enter
enter.index = temp$enter.index
exit = temp$exit
for(t in enter.index)
if( enter[ t ] ) {
out = stop.fn(x[ t ], price, t, exit[ (t + 1) ], ...)
temp = out$state
if( any(temp, na.rm=T) ) {
iexit = t - 1 + dummy[temp][1]
if(out$clean.signal) enter[ t : iexit ] = FALSE
x[ iexit ] = out$value
} else
enter[ t : exit[ t ] ] = FALSE
}
return(x)
}
custom.stop.fn.full <- function(x, price, stop.fn, ...)
{
price = coredata(price)
if(is.character(stop.fn)) stop.fn = match.fun(stop.fn)
temp = bt.ts.trade.index(x)
enter = temp$enter
enter.index = temp$enter.index
exit = temp$exit
for(t in enter.index)
if( enter[ t ] ) {
out = stop.fn(x, price, t, exit[ (t + 1) ], ...)
x = out$x
if(out$clean.signal) enter[ t : out$tlast ] = FALSE
}
return(x)
}
custom.trailing.stop.test <- function(weight, price, tstart, tend, sma, nstop) {
index = tstart : tend
if(weight > 0) {
temp = price[ index ] < cummax(0.9 * sma[ index ])
temp = temp | price[ index ] > cummax(1.1 * sma[ index ])
} else {
temp = price[ index ] > cummax(1.1 * sma[ index ])
}
if( tend - tstart > nstop ) temp[ (nstop + 1) ] = T
return( temp )
}
custom.stop.fn.test <- function() {
signal = c(1,1,1,1,1,1,NA,1,0)
price = c(1,1,0.9,0.7,1,1,1,1,0)
custom.stop.fn(signal, price,
custom.trailing.stop.test,
sma = ifna(SMA(price, 2), price),
nstop = 20
)
signal = -c(1,1,1,1,1,1,NA,1,0)
price = c(1,1,0.9,1.7,1,1,1,1,0)
custom.stop.fn(signal, price,
custom.trailing.stop.test,
sma = ifna(SMA(price, 2), price),
nstop = 4
)
}
plotbt.custom.report <- function
(
...,
dates = NULL,
main = '',
trade.summary = FALSE,
x.highlight = NULL
)
{
ilayout =
'1,1
1,1
2,2
3,3
4,6
4,6
5,7
5,8'
plota.layout(ilayout)
models = variable.number.arguments( ... )
plotbt(models, dates = dates, main = main, plotX = F, log = 'y', LeftMargin = 3, x.highlight = x.highlight)
mtext('Cumulative Performance', side = 2, line = 1)
plotbt(models[1], plottype = '12M', dates = dates, plotX = F, LeftMargin = 3, x.highlight = x.highlight)
mtext('12 Month Rolling', side = 2, line = 1)
plotbt(models[1], dates = dates, xfun = function(x) { 100 * compute.drawdown(x$equity) }, LeftMargin = 3, x.highlight = x.highlight)
mtext('Drawdown', side = 2, line = 1)
model = models[[1]]
name=ifnull(names(models),'')[1]
plotbt.transition.map(model$weight, x.highlight = x.highlight, name=name)
temp = plotbt.monthly.table(model$equity, smain=name)
plotbt.holdings.time(model$weight, smain=name)
if ( !is.null(model$trade.summary) ) {
plot.table( list2matrix(bt.detail.summary(model, model$trade.summary)), keep_all.same.cex = TRUE, smain=name)
} else {
plot.table( list2matrix(bt.detail.summary(model)), keep_all.same.cex = TRUE, smain=name)
}
if( len(models) > 1 ) plotbt.strategy.sidebyside(models)
if ( trade.summary & !is.null(model$trade.summary)) {
ntrades = min(20, nrow(model$trade.summary$trades))
temp = last(model$trade.summary$trades, ntrades)
if( ntrades == 1 ) temp = model$trade.summary$trades
print( temp )
print( model$trade.summary$stats )
layout(c(1,rep(2,10)))
make.table(1,1)
a = matrix(names(models)[1],1,1)
cex = plot.table.helper.auto.adjust.cex(a)
draw.cell(a[1],1,1, text.cex=cex,frame.cell=F)
plot.table( temp )
}
}
plotbt.custom.report.part1 <- function
(
...,
dates = NULL,
main = '',
trade.summary = FALSE,
x.highlight = NULL
)
{
layout(1:3)
models = variable.number.arguments( ... )
model = models[[1]]
plotbt(models, dates = dates, main = main, plotX = F, log = 'y', LeftMargin = 3, x.highlight = x.highlight)
mtext('Cumulative Performance', side = 2, line = 1)
plotbt(models[1], plottype = '12M', dates = dates, plotX = F, LeftMargin = 3, x.highlight = x.highlight)
mtext('12 Month Rolling', side = 2, line = 1)
plotbt(models[1], dates = dates, xfun = function(x) { 100 * compute.drawdown(x$equity) }, LeftMargin = 3, x.highlight = x.highlight)
mtext('Drawdown', side = 2, line = 1)
}
plotbt.custom.report.part2 <- function
(
...,
dates = NULL,
main = '',
trade.summary = FALSE,
x.highlight = NULL
)
{
models = variable.number.arguments( ... )
model = models[[1]]
name=ifnull(names(models),'')[1]
if( len(models) > 1 )
ilayout =
'1,1,3,4
2,2,5,5
2,2,6,6'
else
ilayout =
'1,1,1,3
2,2,4,4
2,2,5,5'
plota.layout(ilayout)
plotbt.transition.map(model$weight, x.highlight = x.highlight, name=name)
temp = plotbt.monthly.table(model$equity, smain=name)
if( len(models) > 1 )
plotbt.holdings.time(model$weight, smain=name)
plot.table(to.percent(t(last(models[[1]]$weight))), smain=name)
if ( !is.null(model$trade.summary) ) {
plot.table( list2matrix(bt.detail.summary(model, model$trade.summary)), keep_all.same.cex = TRUE, smain=name)
} else {
plot.table( list2matrix(bt.detail.summary(model)), keep_all.same.cex = TRUE, smain=name)
}
if( len(models) > 1 )
plotbt.strategy.sidebyside(models)
else
plotbt.holdings.time(model$weight, smain=name)
}
plotbt.custom.report.part3 <- function
(
...,
dates = NULL,
main = '',
trade.summary = FALSE
)
{
models = variable.number.arguments( ... )
model = models[[1]]
if ( trade.summary & !is.null(model$trade.summary)) {
ntrades = min(20, nrow(model$trade.summary$trades))
temp = last(model$trade.summary$trades, ntrades)
if( ntrades == 1 ) temp = model$trade.summary$trades
print( temp )
print( model$trade.summary$stats )
layout(c(1,rep(2,10)))
make.table(1,1)
a = matrix(names(models)[1],1,1)
cex = plot.table.helper.auto.adjust.cex(a)
draw.cell(a[1],1,1, text.cex=cex,frame.cell=F)
plot.table( temp )
}
}
bt.detail.summary <- function
(
bt,
trade.summary = NULL
)
{
out.all = list()
out = list()
out$Period = join( format( range(index.xts(bt$equity)), '%b%Y'), ' - ')
out$Cagr = compute.cagr(bt$equity)
out$Sharpe = compute.sharpe(bt$ret) / 100
out$DVR = compute.DVR(bt) / 100
out$Volatility = compute.risk(bt$ret)
out$MaxDD = compute.max.drawdown(bt$equity)
out$AvgDD = compute.avg.drawdown(bt$equity)
if( !is.null(trade.summary) ) {
out$Profit.Factor = trade.summary$stats['profitfactor', 'All']
}
out$VaR = compute.var(bt$ret)
out$CVaR = compute.cvar(bt$ret)
out$Exposure = compute.exposure(bt$weight)
out.all$System = lapply(out, function(x) if(is.double(x)) round(100*x,2) else x)
if( !is.null(bt$trade.summary) ) trade.summary = bt$trade.summary
out = list()
if( !is.null(trade.summary) ) {
out$Win.Percent = trade.summary$stats['win.prob', 'All']
out$Avg.Trade = trade.summary$stats['avg.pnl', 'All']
out$Avg.Win = trade.summary$stats['win.avg.pnl', 'All']
out$Avg.Loss = trade.summary$stats['loss.avg.pnl', 'All']
out = lapply(out, function(x) if(is.double(x)) round(100*x,1) else x)
out$Best.Trade = max(as.double(trade.summary$trades[, 'return']))
out$Worst.Trade = min(as.double(trade.summary$trades[, 'return']))
out$WinLoss.Ratio = round( -trade.summary$stats['win.avg.pnl', 'All']/trade.summary$stats['loss.avg.pnl', 'All'] , 2)
out$Avg.Len = round(trade.summary$stats['len', 'All'],2)
out$Num.Trades = trade.summary$stats['ntrades', 'All']
}
out.all$Trade = out
out = list()
out$Win.Percent.Day = sum(bt$ret > 0, na.rm = T) / len(bt$ret)
out$Best.Day = bt$best
out$Worst.Day = bt$worst
month.ends = endpoints(bt$equity, 'months')
mret = ROC(bt$equity[month.ends,], type = 'discrete')
out$Win.Percent.Month = sum(mret > 0, na.rm = T) / len(mret)
out$Best.Month = max(mret, na.rm = T)
out$Worst.Month = min(mret, na.rm = T)
year.ends = endpoints(bt$equity, 'years')
mret = ROC(bt$equity[year.ends,], type = 'discrete')
out$Win.Percent.Year = sum(mret > 0, na.rm = T) / len(mret)
out$Best.Year = max(mret, na.rm = T)
out$Worst.Year = min(mret, na.rm = T)
out.all$Period = lapply(out, function(x) if(is.double(x)) round(100*x,1) else x)
return(out.all)
}
engineering.returns.kpi <- function
(
bt,
trade.summary = NULL
)
{
if( !is.null(bt$trade.summary) ) trade.summary = bt$trade.summary
out = list()
out$Period = join( format( range(index(bt$equity)), '%b%Y'), ' - ')
out$Cagr = compute.cagr(bt$equity)
out$Sharpe = compute.sharpe(bt$ret) / 100
out$DVR = compute.DVR(bt) / 100
out$R2 = compute.R2(bt$equity) / 100
out$Volatility = compute.risk(bt$ret)
out$MaxDD = compute.max.drawdown(bt$equity)
out$Exposure = compute.exposure(bt$weight)
if( !is.null(trade.summary) ) {
out$Win.Percent = trade.summary$stats['win.prob', 'All']
out$Avg.Trade = trade.summary$stats['avg.pnl', 'All']
out$Profit.Factor = trade.summary$stats['profitfactor', 'All'] / 100
}
out = lapply(out, function(x) if(is.double(x)) round(100*x,2) else x)
if( !is.null(trade.summary) ) out$Num.Trades = trade.summary$stats['ntrades', 'All']
return( list(System=out))
}
plotbt.strategy.sidebyside <- function
(
... ,
perfromance.metric = spl('System,Trade,Period'),
perfromance.fn = 'bt.detail.summary',
return.table = FALSE,
make.plot = TRUE
)
{
models = variable.number.arguments( ... )
out = list()
for( i in 1:len(models) ) {
out[[ names(models)[i] ]] = match.fun(perfromance.fn)(models[[ i ]])[[ perfromance.metric[1] ]]
}
temp = list2matrix(out, keep.names=F)
if(make.plot) plot.table( temp, smain = perfromance.metric[1] )
if(return.table) return(temp)
}
plotbt <- function
(
...,
dates = NULL,
plottype = spl('line,12M'),
xfun=function(x) { x$equity },
main = NULL,
plotX = T,
log = '',
x.highlight = NULL,
LeftMargin = 0
)
{
models = variable.number.arguments( ... )
plottype = plottype[1]
n = length(models)
temp = list()
for( i in 1:n ) {
msg = try( match.fun(xfun)( models[[i]] ) , silent = TRUE)
if (class(msg)[1] != 'try-error') {
temp[[i]] = msg
}
}
nlag = max( 1, compute.annual.factor(temp[[1]]) )
yrange=c();
for( i in 1:n ) {
itemp = temp[[i]]
if(!is.null(dates)) {
itemp = itemp[dates]
if(itemp[1] != 0) itemp = itemp / as.double(itemp[1])
}
if( plottype == '12M' ) {
itemp = 100 * (itemp / mlag(itemp, nlag ) - 1)
}
temp[[i]] = itemp
yrange = range(yrange, itemp ,na.rm = T)
}
plota(temp[[1]], main = main, plotX = plotX, type = 'l', col = 1,
ylim = yrange,log = log, LeftMargin = LeftMargin, x.highlight = x.highlight)
if( n > 1 ) {
for( i in 2:n ) plota.lines(temp[[i]], col = i)
}
if( plottype == '12M' ) legend('topright', legend = '12 Month Rolling', bty = 'n')
plota.legend(names(models), paste('', 1:n, sep=''), temp)
}
plotbt.transition.map <- function
(
weight,
name = '',
col = rainbow(ncol(weight), start=0, end=.9),
x.highlight = NULL,
sort.asssets = T
)
{
par(mar=c(2, 4, 1, 1), cex = 0.8, cex.main=0.8,cex.sub=0.8,cex.axis=0.8,cex.lab=0.8)
weight[is.na(weight)] = 0
if(sort.asssets) weight = weight[, sort.list(colSums(weight!=0), decreasing=T)]
plota.stacked(index.xts(weight), weight, col = col, type='s', flip.legend=T, main = iif(nchar(name) > 0, paste('Transition Map for', name), ''), x.highlight = x.highlight)
}
plotbt.holdings <- function
(
weight,
smain = format(index.xts(last(weight)), '%d-%b-%Y')
)
{
par(mar=c(2, 2, 2, 2), cex = 0.8, cex.main=0.8,cex.sub=0.8,cex.axis=0.8,cex.lab=0.8)
icols=rainbow(ncol(weight), start=0, end=.9)
weight = weight[, sort.list(colSums(weight!=0, na.rm=T), decreasing=T), drop=F]
temp = 100 * as.vector(last(weight))
atemp = abs(temp)
if(sum(atemp)>0) {
pie(atemp, labels = paste(round(temp,0), '% ', colnames(weight), sep=''),
col = icols, cex =0.8,
main = paste('Allocation for ', smain, sep='')
)
}
}
plotbt.holdings.time <- function(weight, smain='')
{
weight = as.matrix( apply(abs(weight), 2, sum, na.rm = T) )
if( sum(abs(weight)) > 0 ) plotbt.holdings( t(weight) / sum(abs(weight), na.rm = T), smain = paste0(smain, ' in time'))
}
plotbt.monthly.table <- function(equity, make.plot = TRUE, smain = '')
{
equity = map2monthly(equity)
dates = index.xts(equity)
equity = coredata(equity)
if(T) {
month.ends = date.month.ends(dates)
year.ends =  date.year.ends(dates[month.ends])
year.ends = month.ends[year.ends]
nr = len(year.ends) + 1
} else {
month.ends = unique(c(endpoints(dates, 'months'), len(dates)))
month.ends = month.ends[month.ends>0]
year.ends =  unique(c(endpoints(dates[month.ends], 'years'), len(month.ends)))
year.ends = year.ends[year.ends>0]
year.ends = month.ends[year.ends]
nr = len(year.ends) + 1
}
temp = matrix( double(), nr, 12 + 2)
rownames(temp) = c(date.year(dates[year.ends]), 'Avg')
colnames(temp) = spl('Jan,Feb,Mar,Apr,May,Jun,Jul,Aug,Sep,Oct,Nov,Dec,Year,MaxDD')
index = c(1, year.ends)
for(iyear in 2:len(index)) {
iequity = equity[ index[(iyear-1)] : index[iyear] ]
iequity = ifna( ifna.prev(iequity), 0)
temp[(iyear-1), 'Year'] = last(iequity, 1) / iequity[1] -1
temp[(iyear-1), 'MaxDD'] = min(iequity / cummax(iequity) - 1, na.rm = T)
}
index = month.ends
monthly.returns = c(NA, diff(equity[index]) / equity[index[-len(index)]])
index = date.month(range(dates[index]))
monthly.returns = c( rep(NA, index[1]-1), monthly.returns, rep(NA, 12-index[2]) )
temp[1:(nr - 1), 1:12] = matrix(monthly.returns, ncol=12, byrow = T)
temp = ifna(temp, NA)
temp[nr,] = apply(temp[-nr,], 2, mean, na.rm = T)
if(make.plot) {
highlight = temp
highlight[] = iif(temp > 0, 'lightgreen', iif(temp < 0, 'red', 'white'))
highlight[nr,] = iif(temp[nr,] > 0, 'green', iif(temp[nr,] < 0, 'orange', 'white'))
highlight[,13] = iif(temp[,13] > 0, 'green', iif(temp[,13] < 0, 'orange', 'white'))
highlight[,14] = 'yellow'
}
temp[] = plota.format(100 * temp, 1, '', '')
if(make.plot) plot.table(temp, highlight = highlight, smain = smain)
return(temp)
}
list2matrix <- function
(
ilist,
keep.names = TRUE
)
{
if ( is.list( ilist[[1]] ) ) {
inc = 1
if( keep.names ) inc = 2
out = matrix('', nr = max(unlist(lapply(ilist, len))), nc = inc * len(ilist) )
colnames(out) = rep('', inc * len(ilist))
for( i in 1:len(ilist) ) {
nr = len(ilist[[i]])
colnames(out)[inc * i] = names(ilist)[i]
if(nr > 0){
if( keep.names ) {
out[1:nr,(2*i-1)] = names(ilist[[i]])
} else {
rownames(out) = names(ilist[[i]])
}
out[1:nr,inc*i] = unlist(ilist[[i]])
}
}
return(out)
} else {
return( as.matrix(unlist(ilist)) )
}
}
setup.cluster <- function(expr = NULL, varlist = NULL, envir = .GlobalEnv) {
load.packages('parallel')
cores = detectCores()
Sys.unsetenv("R_PROFILE_USER")
cl = makeCluster(cores)
temp = clusterEvalQ(cl, {
library(quantmod)
library(SIT)
NULL
})
if(!is.null(expr))
temp = clusterCall(cl, eval, substitute(expr), env = .GlobalEnv)
if(!is.null(varlist))
clusterExport(cl=cl, spl(varlist),envir=envir)
cl
}
compute.cor <- function
(
data,
method = c("pearson", "kendall", "spearman")
)
{
cor(data, use='complete.obs',method=method[1])
}
proxy.test <- function(data.all, names = ls(data.all), price.fn=Ad)
{
data = new.env()
data$symbolnames = names
for(n in data$symbolnames)
data[[n]] = make.stock.xts( price.fn( data.all[[n]] ) )
bt.prep(data, align='remove.na', fill.gaps=T)
prices = data$prices
layout(1)
plota.matplot(scale.one(prices))
rets = (prices/mlag(prices)-1)[-1,]
temp = cor(rets, use='complete.obs', method='pearson')
diag(temp) = NA
temp[lower.tri(temp)] = NA
temp = temp[-nrow(temp),,drop=F]
temp[] = plota.format(100 * temp, 0, '', '%')
out = temp
temp = compute.stats( as.list(rets),
list(
Mean=function(x) 252*mean(x,na.rm=T),
StDev=function(x) sqrt(252)*sd(x,na.rm=T)
)
)
temp[] = plota.format(100 * temp, 1, '', '%')
out = rbind(out,NA,temp)
print(out)
}
plot12month.rolling.spread <- function(data.all, names = ls(data.all), price.fn=Ad)
{
data = new.env()
data$symbolnames = names[1:2]
for(n in data$symbolnames)
data[[n]] = make.stock.xts( price.fn( data.all[[n]] ) )
bt.prep(data, align='remove.na', fill.gaps=T)
prices = data$prices
rets.12m.rolling = 100 * (prices / mlag(prices, 252) - 1)
spread.12m.rolling = rets.12m.rolling[,1] - rets.12m.rolling[,2]
layout(1)
plota(spread.12m.rolling, type='l',
main = paste('12 Month Rolling Returns Spread % for', names[1], 'and', names[2]))
abline(h=0, col='gray')
}
proxy.overlay.plot <- function(data.all, names = ls(data.all), price.fn=Ad)
{
data = new.env()
data$symbolnames = names
for(n in data$symbolnames)
data[[n]] = make.stock.xts( price.fn( data.all[[n]] ) )
bt.prep(data, align='keep.all', fill.gaps=T)
prices = data$prices
prices = scale.one(prices, T)
layout(1)
plota.matplot(prices)
}
proxy.prices <- function(data, names = ls(data)) {
n.names = len(names)
temp = list()
layout(1:n.names)
for(n in names) {
plota.matplot(cbind(Cl(data[[n]]),Ad(data[[n]])),main=n)
temp[[ paste(n, 'Price') ]] = Cl(data[[n]])
temp[[ paste(n, 'Total') ]] = Ad(data[[n]])
}
temp = compute.stats( lapply(temp, function(x) ifna(x/mlag(x) -1,NA)),
list(
Mean=function(x) 252*mean(x,na.rm=T),
StDev=function(x) sqrt(252)*sd(x,na.rm=T)
)
)
temp[] = plota.format(100 * temp, 1, '', '%')
print(temp)
}
proxy.map <- function(raw.data, tickers)
{
data <- new.env()
tickers = spl(tickers)
tickers = tickers[order(sapply(tickers, nchar),decreasing =T)]
getSymbols.extra(tickers, src = 'yahoo', from = '1980-01-01', env = data, raw.data = raw.data, set.symbolnames = T, auto.assign = T)
for(i in data$symbolnames) data[[i]] = adjustOHLC(data[[i]], use.Adjusted=T)
bt.prep(data, align='keep.all', fill.gaps=T)
layout(1)
plota.matplot(data$prices)
}
proxy.example.test <- function() {
load.packages('quantmod')
tickers = spl('GSG,DBC')
data = new.env()
getSymbols(tickers, src = 'yahoo', from = '1970-01-01', env = data, auto.assign = T)
for(i in ls(data)) data[[i]] = adjustOHLC(data[[i]], use.Adjusted=T)
temp = extract.table.from.webpage( join(readLines("TRJ_CRB")), 'EODValue' )
temp = join( apply(temp, 1, join, ','), '\n' )
data$CRB_1 = make.stock.xts( read.xts(temp, format='%m/%d/%y' ) )
data$CRB_2 = make.stock.xts( read.xts("prfmdata.csv", format='%m/%d/%Y' ) )
jpeg(filename = 'plot1.jpg', width = 500, height = 500, units = 'px', pointsize = 12)
proxy.test(data)
dev.off()
jpeg(filename = 'plot2.jpg', width = 500, height = 500, units = 'px', pointsize = 12)
proxy.overlay.plot(data)
dev.off()
load.packages('quantmod')
tickers = spl('IYR,VGSIX,RWO')
data = new.env()
getSymbols(tickers, src = 'yahoo', from = '1970-01-01', env = data, auto.assign = T)
for(i in ls(data)) data[[i]] = adjustOHLC(data[[i]], use.Adjusted=T)
jpeg(filename = 'plot3.jpg', width = 500, height = 500, units = 'px', pointsize = 12)
proxy.test(data)
dev.off()
jpeg(filename = 'plot4.jpg', width = 500, height = 500, units = 'px', pointsize = 12)
proxy.overlay.plot(data)
dev.off()
}
find.tokens <- function
(
txt,
marker,
pos = 1,
pos.start = T
)
{
marker = spl(marker)
for(i in 1:len(marker)) {
if( pos < 2 )
pos1 = regexpr(marker[i], txt)
else
pos1 = regexpr(marker[i], substr(txt, pos, nchar(txt)))
if( pos1 < 0 )
return(pos1)
else {
if( pos < 2 ) pos = pos1
else pos = pos1 + pos - 1
}
pos = pos + attr(pos1, 'match.length')
}
if( pos.start ) pos = pos - attr(pos1, 'match.length')
return(pos)
}
extract.token <- function
(
txt,
smarker,
emarker,
pos = 1,
keep.marker = F
)
{
pos1 = 1
if (nchar(smarker) > 0)
pos1 = find.tokens(txt, smarker, pos, pos.start = keep.marker)
if( pos1 < 0 ) return("")
pos1.marker = iif(keep.marker, pos1 + nchar(last(spl(smarker))), pos1)
pos2 = nchar(txt)
if (nchar(emarker) > 0)
pos2 = find.tokens(txt, emarker, pos1.marker, pos.start = !keep.marker) - 1
if( pos2 < 0 ) return("")
return(substr(txt,pos1,pos2))
}
remove.tags <- function
(
temp
)
{
temp = gsub(pattern = '<.*?>', replacement = '', temp, perl = TRUE)
temp = gsub(pattern = '\r', replacement = '', temp, perl = TRUE)
temp = gsub(pattern = '\n', replacement = '', temp, perl = TRUE)
temp = gsub(pattern = '\t', replacement = '', temp, perl = TRUE)
temp = gsub(pattern = '&nbsp;', replacement = '', temp, perl = TRUE)
temp = gsub(pattern = '&amp;', replacement = '', temp, perl = TRUE)
temp = gsub(pattern = '&raquo;', replacement = '', temp, perl = TRUE)
temp = gsub(pattern = '&#37;', replacement = '%', temp, perl = TRUE)
return(temp)
}
replace.token <- function
(
txt,
smarker,
emarker,
replacement,
pos = 1
)
{
token = extract.token(txt, smarker, emarker, pos, keep.marker = T)
if(nchar(token) == 0)
txt
else
replace.token(gsub(pattern = token, replacement = replacement, txt), smarker, emarker, replacement)
}
clean.table <- function
(
temp
)
{
temp = trim(temp)
temp[nchar(temp)==0] = NA
temp = temp[ncol(temp) > rowSums(is.na(temp)),,drop=F]
temp[,nrow(temp) > colSums(is.na(temp)),drop=F]
}
extract.table.from.webpage <- function
(
txt,
marker=NA,
has.header=T,
end.marker=NA
)
{
tryCatch({
pos1=1
if(!is.na(marker)) {
marker = spl(marker)
if(len(marker) > 0 && nchar(marker[1]) > 0)
for(i in 1:len(marker))
pos1 = regexpr(marker[i], substr(txt, pos1, nchar(txt))) + pos1
}
pos0 = tail(gregexpr('<table', substr(txt, 1, pos1))[[1]], 1)
if(pos0 == -1) pos0 = tail(gregexpr('<tbody', substr(txt, 1, pos1))[[1]], 1)
if(pos0 == -1) pos0 = pos1
pos2 = head(gregexpr('</table', substr(txt, pos1, nchar(txt)))[[1]], 1)
if(pos2 == -1) pos2 = head(gregexpr('</tbody', substr(txt, pos1, nchar(txt)))[[1]], 1)
if(pos2 == -1) pos2 = nchar(txt)+1
temp =  substr(txt, pos0, pos1 + pos2 - 2)
temp = gsub(pattern = '<br>', replacement = '', temp, perl = TRUE)
temp = gsub(pattern = '</tr>', replacement = ';row;', temp, perl = TRUE)
temp = gsub(pattern = '</td>', replacement = ';col;', temp, perl = TRUE)
temp = gsub(pattern = '</th>', replacement = ';col;', temp, perl = TRUE)
if(!is.na(end.marker)) {
marker = spl(end.marker)
if(len(marker) > 0 && nchar(marker[1]) > 0)
for(i in 1:len(marker))
temp = gsub(pattern = marker[i], replacement = ';row;', temp, perl = TRUE)
}
temp = gsub(pattern = '<.*?>', replacement = '', temp, perl = TRUE)
temp = gsub(pattern = '\r', replacement = '', temp, perl = TRUE)
temp = gsub(pattern = '\n', replacement = '', temp, perl = TRUE)
temp = gsub(pattern = '\t', replacement = '', temp, perl = TRUE)
temp = gsub(pattern = '&nbsp;', replacement = '', temp, perl = TRUE)
temp = gsub(pattern = '&amp;', replacement = '', temp, perl = TRUE)
temp = gsub(pattern = '&raquo;', replacement = '', temp, perl = TRUE)
temp = lapply( strsplit(temp, ';row;'), strsplit, ';col;')
n = max( sapply(temp[[1]], function(x) len(x)) )
temp = t( sapply(temp[[1]], function(x) x[1:n]) )
if(has.header) {
colnames(temp) = trim(temp[(has.header + 0), ])
temp = temp[-c(1:(has.header + 0)), ,drop=F]
}
}, error = function(ex) {
temp <<- txt
}, finally = {
return(temp)
})
}
extract.table.from.webpage.test <- function()
{
load.packages('quantmod')
Symbol = 'IBM'
url = paste('http://finance.yahoo.com/q/ks?s=', Symbol, sep = '')
txt = join(readLines(url))
temp = extract.table.from.webpage(txt, 'Market Cap', has.header = F)
temp = rbind(c('', Symbol), temp)
data = getSymbols(Symbol, from = '1980-01-01', auto.assign = FALSE)
y = data['2010::2011']
sma50 = SMA(Cl(y), 50)
png(filename = 'plot1.png', width = 500, height = 500, units = 'px', pointsize = 12, bg = 'white')
layout(c(1,1,2,3,3))
plota(y, type = 'candle', main = Symbol, plotX = F)
plota.lines(sma50, col='blue')
plota.legend(c(Symbol,'SMA 50'), 'green,blue', list(y,sma50))
y = plota.scale.volume(y)
plota(y, type = 'volume')
plot.table(temp)
dev.off()
}
PricingZeroCouponBond <- function
(
yield,
timetomaturity,
parvalue = 100
)
{
parvalue / ( 1 + yield ) ^ timetomaturity
}
processTBill <- function
(
yields,
timetomaturity = 1/4,
frequency = 365
)
{
yield = coredata(yields) / 100
pr = sapply( yield, function(x) PricingZeroCouponBond(x, timetomaturity) )
pr = ROC(pr, type='discrete')
pr[1] = 0
ir = (1+mlag(yield, nlag=1))^(1 / frequency)-1
ir[1] = 0
tr = pr + ir
close.price = cumprod(1 + pr)
adjusted.price = cumprod(1 + tr)
out = as.xts( cbind(close.price, adjusted.price), index(yields) )
colnames(out) = spl('Close,Adjusted')
return(out)
}
processTBill.test <- function()
{
quantmod::getSymbols("GS1", src = "FRED")
ir = (1 + mlag(GS1) / 100) ^ (1/12) - 1
ir[1] = 0
out = processTBill(GS1, timetomaturity = 1,12)
plota(cumprod(1 + ir), type='l', log = 'y')
plota.lines(Ad(out), type='l', col='red')
SHY = getSymbols('SHY', src='yahoo', auto.assign = FALSE)
tbill.m = quantmod::getSymbols('GS3', src='FRED', auto.assign = FALSE)
tbill.d = quantmod::getSymbols('DGS3', src='FRED', auto.assign = FALSE)
timetomaturity = 3
compute.raw.annual.factor(tbill.d)
compute.raw.annual.factor(tbill.m)
tbill.m = processTBill(tbill.m, timetomaturity = timetomaturity, 12)
tbill.d[] = ifna.prev(tbill.d)
tbill.d = processTBill(tbill.d, timetomaturity = timetomaturity,261)
dates = '2003::'
tbill.m = tbill.m[dates,2]
tbill.m = tbill.m / as.double(tbill.m[1])
tbill.d = tbill.d[dates,2]
tbill.d = tbill.d / as.double(tbill.d[1])
SHY = Ad(SHY[dates,])
SHY = SHY / as.double(SHY[1])
plota(tbill.d, type='l')
plota.lines(tbill.m, type='s', col='blue')
plota.lines(SHY, type='l', col='red')
plota.legend('Daily 3YR T-Bills,Monthly 3YR T-Bills,SHY','black,blue,red')
}
get.CRB <- function(...)
{
load.packages('gtools,gdata')
url = paste('http://www.jefferies.com/html/ProductsServices/SalesTrading/Commodities/scripts/genExcel.pl?Index=RJCRB_Total&StartDate=19940101&EndDate=', format(Sys.Date(), '%Y%m%d'), sep='')
temp = read.xls(url, ...)
temp = as.matrix(temp[-c(1:7),])
out = repmat(as.double(temp[,2]), 1, 6)
colnames(out) = spl('Open,High,Low,Close,Volume,Adjusted')
out[, 'Volume'] = 0
out = make.xts( out, as.POSIXct(temp[,1], tz = Sys.getenv('TZ'), format='%m/%d/%y'))
indexClass(out) = 'Date'
return(out)
}
get.CRB.test <- function()
{
CRB = get.CRB()
load.packages('quantmod')
tickers = spl('GSG,DBC')
getSymbols(tickers, src = 'yahoo', from = '1970-01-01')
out = na.omit(merge(Cl(CRB), Cl(GSG), Cl(DBC)))
colnames(out) = spl('CRB,GSG,DBC')
temp = out / t(repmat(as.vector(out[1,]),1,nrow(out)))
layout(1:2)
plota(temp, ylim=range(temp))
plota.lines(temp[,1],col=1)
plota.lines(temp[,2],col=2)
plota.lines(temp[,3],col=3)
plota.legend(colnames(temp),1:3)
temp = cor(temp / mlag(temp)- 1, use='complete.obs', method='pearson')
temp[] = plota.format(100 * temp, 0, '', '%')
plot.table(temp)
layout(1:3)
plota.matplot(CRB[,c('Close','Adjusted')])
plota.matplot(DBC[,c('DBC.Close','DBC.Adjusted')])
plota.matplot(GSG[,c('GSG.Close','GSG.Adjusted')])
layout(1)
comm = extend.data(DBC, CRB, scale=T)
plota(comm, type='l', col=1)
plota.lines(CRB*0.078, type='l', lwd=5, col=col.add.alpha(2,150))
plota.lines(DBC, type='l', lwd=5, col=col.add.alpha(3,150))
plota.lines(comm, type='l', col=1)
plota.legend('comm,CRB,DBC', 1:3, list(comm,CRB,DBC))
}
dow.jones.components <- function()
{
url = 'http://money.cnn.com/data/dow30/'
txt = join(readLines(url))
temp = gsub(pattern = '">', replacement = '<td>', txt, perl = TRUE)
temp = gsub(pattern = '</a>', replacement = '</td>', temp, perl = TRUE)
temp = extract.table.from.webpage(temp, 'Volume', has.header = T)
trim(temp[,'Company'])
}
dow.jones.components.0 <- function()
{
url = 'http://finance.yahoo.com/q/cp?s=^DJI+Components'
txt = join(readLines(url))
temp = extract.table.from.webpage(txt, 'Volume', has.header = T)
temp[, 'Symbol']
}
dow.jones.components.1 <- function()
{
load.packages('readxl,httr')
dir.create(paste(getwd(), 'temp', sep='/'), F)
GET('http://www.djaverages.com/?go=export-components&symbol=DJI', write_disk('temp/DJI.xls', overwrite=TRUE))
temp = read_excel('temp/DJI.xls')
temp$Ticker
}
nasdaq.100.components <- function()
{
url = 'http://www.nasdaq.com/quotes/nasdaq-100-stocks.aspx?render=download'
temp = read.csv(url, header=TRUE, stringsAsFactors=F)
colnames(temp) = trim(colnames(temp))
tickers = temp[, 'Symbol']
return(tickers)
}
sector.spdr.components <- function(sector.etf = 'XLE')
{
url = paste('http://www.sectorspdr.com/sectorspdr/IDCO.Client.Spdrs.Holdings/Export/ExportCsv?symbol=', sector.etf, sep='')
temp = read.csv(url, skip=1, header=TRUE, stringsAsFactors=F)
tickers = temp[, 'Symbol']
return(tickers)
}
sp500.components <- function()
{
url = 'http://en.wikipedia.org/wiki/List_of_S%26P_500_companies'
txt = join(readLines(url))
temp = extract.table.from.webpage(txt, 'Ticker', has.header = T)
tickers = temp[, 'Ticker symbol']
sector = temp[, 'GICS Sector']
return(list(tickers=tickers, sector=sector))
}
sp100.components <- function()
{
url = 'http://www.barchart.com/stocks/sp100.php'
txt = join(readLines(url))
temp = extract.table.from.webpage(txt, 'Components', has.header = T)
i.start = grep('Name', temp[,2])
tickers = trim(temp[-c(1:i.start), 1])
return(tickers)
}
ftse100.components <- function()
{
url = 'http://uk.ishares.com/en/rc/products/ISF/all-holdings/'
txt = join(readLines(url))
txt = gsub('&#37;','%',txt)
temp = extract.table.from.webpage(txt, 'Security', has.header = T)
temp = trim(temp)
colnames(temp) = temp[1,]
temp = temp[-1,]
holdings = temp
page.label = ''
ticker2ISIN = c()
for(i in 1:100) {
cat(i,'\n')
url = paste('http://www.londonstockexchange.com/exchange/prices-and-markets/stocks/indices/constituents-indices.html?index=UKX&page=', i, sep='')
txt = join(readLines(url))
pos = regexpr('Page [0-9]+ of [0-9]+', txt, ignore.case = T)
page.label.new = substr(txt, pos, pos + attr(pos, 'match.length')-1)
if(page.label == page.label.new) break
page.label = page.label.new
temp.table = extract.table.from.webpage(txt, 'Price', has.header = T)
colnames(temp.table)[1] = 'tickers'
temp = gsub(pattern = '<a', replacement = '<td>', txt, perl = TRUE)
temp = gsub(pattern = '</a>', replacement = '</td>', temp, perl = TRUE)
temp = extract.table.from.webpage(temp, 'Price', has.header = T)
pos = regexpr('fourWayKey=', temp[,2])
ISIN = as.vector(sapply(1:nrow(temp), function(j)
substr(temp[j,2], pos[j] + attr(pos, 'match.length')[j], pos[j] + attr(pos, 'match.length')[j] + 12 - 1)
))
ticker2ISIN = rbind(ticker2ISIN, cbind(temp.table[,spl('ticker,Name,Price'), drop=F], ISIN))
}
ISIN = intersect(holdings[,'ISIN'],ticker2ISIN[,'ISIN'])
holdings = cbind(holdings[match(ISIN, holdings[,'ISIN']), ],
ticker2ISIN[match(ISIN, ticker2ISIN[,'ISIN']), spl('ticker,Name,Price')])
return(apply(holdings, 2, list))
}
us.ishares.components <- function(Symbol = 'DVY', date = NULL, debug = F)
{
url = paste('http://us.ishares.com/product_info/fund/holdings/', Symbol, '.htm?periodCd=d', sep='')
if( !is.null(date) )
url = paste('http://us.ishares.com/product_info/fund/holdings/', Symbol, '.htm?asofDt=', date.end(date), '&periodCd=m', sep='')
txt = join(readLines(url))
temp = remove.tags(extract.token(txt, 'Holdings Detail', 'Holdings subject to change'))
date = as.Date(spl(trim(temp),' ')[3], '%m/%d/%Y')
temp = extract.table.from.webpage(txt, 'Symbol', has.header = T)
colnames(temp) = trim(colnames(temp))
temp = trim(temp)
tickers = temp[, 'Symbol']
keep.index = nchar(tickers)>1
weights = as.double(temp[keep.index, '% Net Assets']) / 100
tickers = tickers[keep.index]
out = list(tickers = tickers, weights = weights, date = date)
if(debug) out$txt = txt
out
}
google.search <- function
(
query
)
{
url = paste("http://google.com/search?ie=utf-8&oe=utf-8&q=", URLencode(query), "&num=10&gws_rd=cr", sep='')
txt = join(readLines(url))
tokens = spl(txt, '<li class="g">')
if(len(tokens) < 2) return(NULL)
records = matrix('', nrow=len(tokens)-1,nc=2)
colnames(records) = c('label','url')
for(i in 2:len(tokens)) {
token = tokens[i]
token = extract.token(token, '<a href=', '</a>', keep.marker = T)
url = extract.token(token, 'url\\?q=', '&amp;sa=U&amp;')
label = remove.tags(token)
records[i-1,] = c(label,url)
}
return(records)
}
getQuote.google <- function(tickers) {
url = paste('http://finance.google.com/finance/info?client=ig&q=', join(tickers,','), sep='')
txt = join(readLines(url))
temp = gsub(':', ',', txt)
temp = scan(text = temp, what='', sep=',', quiet=T)
temp = matrix(trim(temp), nr=len(temp)/len(tickers), byrow=F)
index = match(spl('t,l,lt'), tolower(temp[,1]))+1
names(index) = spl('ticker,last,date')
last = as.double(temp[index['last'],])
date = strptime(temp[index['date'],],format=' %b %d, %H,%M')
out = data.frame(last,date)
rownames(out) = temp[index['ticker'],]
out
}
getQuote.google.xml <- function(tickers) {
url = paste('http://www.google.com/ig/api?', paste('stock=',tickers, '&', sep='', collapse=''), sep='')
txt = join(readLines(url))
temp = txt
temp = gsub('<finance.*?>', '', temp, perl = TRUE)
temp = gsub('</finance>', '', temp, perl = TRUE)
temp = gsub('<xml.*?>', '', temp, perl = TRUE)
temp = gsub('</xml.*?>', '', temp, perl = TRUE)
temp = gsub('<\\?xml.*?>', '', temp, perl = TRUE)
temp = gsub('data=', '', temp, perl = TRUE)
temp = gsub('/><', ' ', temp)
temp = gsub('>', '', temp)
temp = gsub('<', '', temp)
temp = scan(text = temp, what='', sep=' ', quiet=T)
temp = matrix(trim(temp), nr=len(temp)/len(tickers), byrow=F)
cnames = spl('trade_date_utc,trade_time_utc,symbol,last,high,low,volume,open,avg_volume,market_cap,y_close')
index = match(cnames, tolower(temp[,1]))+1
names(index) = cnames
date = strptime(paste(temp[index['trade_date_utc'],], temp[index['trade_time_utc'],]), format='%Y%m%d %H%M%S',tz='UTC')
date = as.POSIXct(date, tz = Sys.getenv('TZ'))
out = data.frame(t(temp[index[-c(1:3)],]))
colnames(out) = cnames[-c(1:3)]
rownames(out) = temp[index['symbol'],]
out
}
getSymbol.intraday.google <- function
(
Symbol,
Exchange,
interval = 60,
period = '1d'
)
{
url = paste('http://www.google.com/finance/getprices?q=', Symbol,
'&x=', Exchange,
'&i=', interval,
'&p=', period,
'&f=', 'd,o,h,l,c,v', sep='')
load.packages('data.table')
out = fread(url, stringsAsFactors=F)
if(ncol(out) < 5) {
cat('Error getting data from', url, '\n')
return(NULL)
}
setnames(out, spl('Date,Open,High,Low,Close,Volume'))
date = out$Date
date.index = substr(out$Date,1,1) == 'a'
date = as.double(gsub('a','',date))
temp = NA * date
temp[date.index] = date[date.index]
temp = ifna.prev(temp)
date = temp + date * interval
date[date.index] = temp[date.index]
class(date) = c("POSIXt", "POSIXct")
date = date - (as.double(format(date[1],'%H')) - 9)*60*60
make.xts(out[,-1,with=F], date)
}
getQuote.yahoo.today <- function(Symbols) {
require('data.table')
what = yahooQF(names = spl('Name,Symbol,Last Trade Time,Last Trade Date,Open,Days High,Days Low,Last Trade (Price Only),Volume,Previous Close'))
names = spl('Name,Symbol,Time,Date,Open,High,Low,Close,Volume,Yesterday')
all.symbols = lapply(seq(1, len(Symbols), 100), function(x) na.omit(Symbols[x:(x + 99)]))
out = c()
for(i in 1:len(all.symbols)) {
url = paste('http://download.finance.yahoo.com/d/quotes.csv?s=',
join( trim(all.symbols[[i]]), ','),
'&f=', what[[1]], sep = '')
txt = join(readLines(url),'\n')
data = fread(paste0(txt,'\n'), stringsAsFactors=F, sep=',')
setnames(data,names)
setkey(data,'Symbol')
out = rbind(out, data)
}
out
}
extend.GLD <- function(GLD) {
extend.data(GLD, bundes.bank.data.gold(), scale=T)
}
extend.data <- function
(
current,
hist,
scale = F
)
{
colnames(current) = sapply(colnames(current), function(x) last(spl(x,'\\.')))
colnames(hist) = sapply(colnames(hist), function(x) last(spl(x,'\\.')))
close.index = find.names('Close', hist)
if(len(close.index)==0) close.index = 1
adjusted.index = find.names('Adjusted', hist)
if(len(adjusted.index)==0) adjusted.index = close.index
if(scale) {
cur.close.index = find.names('Close', current)
if(len(cur.close.index)==0) cur.close.index = 1
cur.adjusted.index = find.names('Adjusted', current)
if(len(cur.adjusted.index)==0) cur.adjusted.index = cur.close.index
common = merge(current[,cur.close.index], hist[,close.index], join='inner')
scale = as.numeric(common[1,1]) / as.numeric(common[1,2])
if( close.index == adjusted.index )
hist = hist * scale
else {
hist[,-adjusted.index] = hist[,-adjusted.index] * scale
common = merge(current[,cur.adjusted.index], hist[,adjusted.index], join='inner')
scale = as.numeric(common[1,1]) / as.numeric(common[1,2])
hist[,adjusted.index] = hist[,adjusted.index] * scale
}
}
hist = hist[format(index(current[1])-1,'::%Y:%m:%d'),,drop=F]
if( ncol(hist) != ncol(current) )
hist = rep.col(hist[,adjusted.index], ncol(current))
else
hist = hist[, colnames(current)]
colnames(hist) = colnames(current)
rbind( hist, current )
}
extend.data.proxy <- function(data, data.proxy = NULL, proxy.filename = 'data.proxy.Rdata') {
if(is.null(data.proxy) && file.exists(proxy.filename))
load(file=proxy.filename)
if(!is.null(data.proxy))
for(n in ls(data.proxy))
if( !is.null(data[[n]]) )
data[[n]] = extend.data(data[[n]], data.proxy[[n]], scale=T)
}
create.leveraged = function(hist, leverage=2) {
rets = 1 + leverage * (hist / mlag(hist) - 1)
rets[1,] = 1
bt.apply.matrix(rets, cumprod)
}
create.leveraged.test = function() {
tickers = spl('TMF,UBT,TLT')
data = new.env()
getSymbols(tickers, src = 'yahoo', from = '1970-01-01', env = data, auto.assign = T)
test2 = extend.data(data$UBT, create.leveraged(data$TLT, leverage=2), scale=T)
test3 = extend.data(data$TMF, create.leveraged(data$TLT, leverage=3), scale=T)
proxy.test(list(TLT=data$TLT, UBT=test2, TMF=test3),price.fn=Ad)
test0 = create.leveraged(data$TLT, leverage=2)
proxy.test(list(UBT=data$UBT, EXTEND=test0),price.fn=Ad)
test0 = create.leveraged(data$TLT, leverage=3)
proxy.test(list(TMF=data$TMF, EXTEND=test0),price.fn=Ad)
}
bundes.bank.data <- function(symbol) {
url = paste('http://www.bundesbank.de/cae/servlet/CsvDownload?tsId=', symbol, '&its_csvFormat=en&mode=its', sep='')
temp = read.csv(url, skip=5, header=F, stringsAsFactors=F)
hist = make.xts(as.double(temp[,2]), as.POSIXct(temp[,1], tz = Sys.getenv('TZ'), format='%Y-%m-%d'))
indexClass(hist) = 'Date'
colnames(hist)='Close'
return( hist[!is.na(hist)] )
}
bundes.bank.data.gold <- function() {
bundes.bank.data('BBEX3.D.XAU.USD.EA.AC.C05')
}
fx.sauder.data <- function(start.year, end.year, base.cur, target.curs) {
url = paste('http://fx.sauder.ubc.ca/cgi/fxdata?b=', base.cur, join(paste('&c=', spl(target.curs), sep='')), '&rd=&fd=1&fm=1&fy=', start.year, '&ld=31&lm=12&ly=', end.year, '&y=daily&q=volume&f=csv&o=', sep='')
temp = read.csv(url, skip=1, header=T, stringsAsFactors=F)
hist = make.xts(as.matrix(temp[,-c(1:3)]), as.POSIXct(temp[,2], tz = Sys.getenv('TZ'), format='%Y/%m/%d'))
indexClass(hist) = 'Date'
colnames(hist) = gsub(paste('.', base.cur, sep=''), '', colnames(hist))
return( hist[!is.na(hist[,1]),] )
}
getSymbols.PI <- function
(
Symbols,
env = .GlobalEnv,
auto.assign = TRUE,
download = TRUE
)
{
temp.folder = paste(getwd(), 'temp', sep='/')
dir.create(temp.folder, F)
for (i in 1:len(Symbols)) {
if(download) {
url = paste('http://pitrading.com/free_eod_data/', Symbols[i], '.zip', sep='')
filename = paste(temp.folder, '/', Symbols[i], '.zip', sep='')
download.file(url, filename,  mode = 'wb')
unzip(filename, exdir=temp.folder)
}
filename = paste(temp.folder, '/', Symbols[i], '.txt', sep='')
temp = read.delim(filename, header=TRUE, sep=',')
out = make.xts(temp[,-1], as.POSIXct(temp[,1], tz = Sys.getenv('TZ'), format='%m/%d/%Y'))
indexClass(out) = 'Date'
out$Adjusted = out$Close
cat(i, 'out of', len(Symbols), 'Reading', Symbols[i], '\n', sep='\t')
if (auto.assign) {
assign(paste(gsub('\\^', '', Symbols[i]), sep='_'), out, env)
}
}
if (!auto.assign) {
return(out)
} else {
return(env)
}
}
getSymbols.fxhistoricaldata <- function
(
Symbols,
type = spl('hour,day'),
env = .GlobalEnv,
auto.assign = TRUE,
download = FALSE,
name.has.type = TRUE
)
{
type = type[1]
type0 = paste0(type,'_')
temp.folder = paste(getwd(), 'temp', sep='/')
dir.create(temp.folder, F)
for (i in 1:len(Symbols)) {
if(download) {
url = paste('http://www.fxhistoricaldata.com/download/', Symbols[i], '_', type, '.zip', sep='')
filename = paste(temp.folder, '/', Symbols[i], '_', type, '.zip', sep='')
download.file(url, filename,  mode = 'wb')
unzip(filename, exdir=temp.folder)
}
filename = paste(temp.folder, '/', Symbols[i], '_', type, '.csv', sep='')
temp = read.delim(filename, header=TRUE, sep=',')
colnames(temp) = gsub('[X\\.|\\.]', '', colnames(temp))
out = make.xts(temp[,spl('OPEN,LOW,HIGH,CLOSE')],
strptime(paste(temp$DATE, temp$TIME), format='%Y%m%d %H:%M:%S'))
cat(i, 'out of', len(Symbols), 'Reading', Symbols[i], '\n', sep='\t')
if (auto.assign) {
assign(paste0(gsub('\\^', '', Symbols[i]), iif(name.has.type,type0,'')), out, env)
}
}
if (!auto.assign) {
return(out)
} else {
return(env)
}
}
get.G10 <- function
(
type = spl('currency')
)
{
if( type[1] != 'currency') {
cat('Warning:', type[1], 'is not yet implemented in getG10 function\n')
return()
}
map = '
FX          FX.NAME
DEXUSAL     U.S./Australia
DEXUSUK     U.S./U.K.
DEXCAUS     Canada/U.S.
DEXNOUS     Norway/U.S.
DEXUSEU     U.S./Euro
DEXJPUS     Japan/U.S.
DEXUSNZ     U.S./NewZealand
DEXSDUS     Sweden/U.S.
DEXSZUS     Switzerland/U.S.
'
map = matrix(scan(text = map, what='', quiet=T), nc=2, byrow=T)
colnames(map) = map[1,]
map = data.frame(map[-1,], stringsAsFactors=F)
convert.index = grep('DEXUS',map$FX, value=T)
load.packages('quantmod')
data.fx <- new.env()
quantmod::getSymbols(map$FX, src = 'FRED', from = '1970-01-01', env = data.fx, auto.assign = T)
for(i in ls(data.fx)) data.fx[[i]] = na.omit(data.fx[[i]])
for(i in convert.index) data.fx[[i]] = 1 / data.fx[[i]]
bt.prep(data.fx, align='remove.na')
fx = bt.apply(data.fx, '[')
return(fx)
}
wealthsimple.portfolio = function(portfolio.number = 10) {
url = paste0('http://faq.wealthsimple.com/article/', 120+portfolio.number, '-how-has-the-risk-level-',portfolio.number,'-portfolio-performed')
txt = join(readLines(url))
temp = extract.table.from.webpage(txt, 'Breakdown', has.header = F)
temp = gsub(pattern = '%', replacement = '', temp)
temp  = trim(temp[,c(2,4)])
temp  = temp[!is.na(temp[,1]),]
value = as.numeric(temp[,2])
names(value) = temp[,1]
value
}
wealthsimple.portfolio.test = function() {
portfolios = list()
for(i in 1:10)
portfolios[[i]] = wealthsimple.portfolio(i)
portfolios = t(sapply(portfolios, identity))
plota.stacked(1:10, portfolios/100, flip.legend = T, type='s', xaxp=c(1,10,9), las=1,
main='Wealthsimple Transition Matrix', xlab='Risk Portfolio')
}
getSymbols.TB <- function(
env = .GlobalEnv,
auto.assign = TRUE,
download = FALSE,
type = c('Both', 'Futures', 'Forex'),
rm.index =  'PB',
clean = FALSE,
custom.adjustments = TRUE
)
{
if(download) {
download.file('http://www.tradingblox.com/Data/DataOnly.zip', 'DataOnly.zip')
}
temp.folder = paste(getwd(), 'temp', sep='/')
dir.create(temp.folder, F)
temp.folder = paste(getwd(), '/', 'temp', sep='')
if(clean) shell('del /F /S /Q temp\\*.*', wait = TRUE)
files = unzip('DataOnly.zip', exdir=temp.folder)
def1 = try(read.csv('http://www.tradingblox.com/tradingblox/CSIUA/FuturesInfo.txt',skip=1,header=FALSE, stringsAsFactors=F),TRUE)
if(inherits(def1, 'try-error')) def1 = read.csv('FuturesInfo.txt',skip=1,header=FALSE, stringsAsFactors=F)
def1 = def1[-match(rm.index, def1[,1]),]
def1[,3] = 'Futures'
def2 = try(read.csv('http://www.tradingblox.com/tradingblox/CSIUA/ForexInfo.txt',skip=1,header=FALSE, stringsAsFactors=F),TRUE)
if(inherits(def2, 'try-error')) def2 = read.csv('ForexInfo.txt',skip=1,header=FALSE, stringsAsFactors=F)
def2[,3] = 'Forex'
def = rbind(def1[,1:4], def2[,1:4])
if(type[1] == 'Futures') def = def1[,1:4]
if(type[1] == 'Forex') def = def2[,1:4]
for( i in 1:nrow(def) ) {
symbol = def[i,1]
filename = paste(temp.folder, '/', def[i,3], '/', def[i,4], sep='')
if(file.exists(filename)) {
fr <- read.csv(filename, header = FALSE)
fr <- make.xts(fr[,-1], as.Date(as.character(fr[,1]),'%Y%m%d'))
colnames(fr) <- spl('Open,High,Low,Close,Volume,OpenInterest,DeliveryMonth,Unadjusted')[1:ncol(fr)]
fr$Adjusted = fr$Close
if (auto.assign) assign(symbol, fr, env)
cat(i, 'out of', nrow(def), 'Reading', symbol, format(index.xts(fr)[1],'%Y%m%d'), format(index.xts(fr)[nrow(fr)],'%Y%m%d'), '\n', sep='\t')
} else {
cat('\t\t\t Missing data for ', symbol, '\n');
}
}
index = match(ls(env)[ na.omit(match(def[,1], ls(env))) ], def[,1])
temp = def[index,1]
names(temp) = def[index,1]
env$symbolnames = temp
temp = def[index,2]
names(temp) = def[index,1]
env$symbol.descriptions = temp
temp = def[index,3]
names(temp) = def[index,1]
env$symbol.groups = temp
names = trim(gsub(pattern = '\\(.*?\\)', replacement = '', env$symbol.descriptions, perl = TRUE))
names = trim(gsub('-NYMEX','',names,ignore.case =T))
names = trim(gsub('-COMEX','',names,ignore.case =T))
names = trim(gsub('-CBT','',names,ignore.case =T))
names = trim(gsub('-CME-','',names,ignore.case =T))
names = trim(gsub('-CME','',names,ignore.case =T))
names = trim(gsub('-NYCE','',names,ignore.case =T))
names = trim(gsub('-Globex','',names,ignore.case =T))
names = trim(gsub('-FINEX','',names,ignore.case =T))
names = trim(gsub('-CSCE','',names,ignore.case =T))
names = trim(gsub(' w/Prj A','',names,ignore.case =T))
env$symbol.descriptions.print = names
data = env
if( custom.adjustments ) {
for(i in data$symbolnames[data$symbol.groups != 'Forex']) {
spot = as.vector(data[[i]]$Unadjusted)
dspot = spot - mlag(spot)
futures = as.vector(data[[i]]$Adjusted)
dfutures = futures - mlag(futures)
index = which(round(dspot - dfutures,4) != 0 )
spot.adjust.roll = spot
spot.adjust.roll[(index-1)] = spot.adjust.roll[index] - dfutures[index]
reta = dfutures / mlag(spot.adjust.roll)
reta[1] = 0
n = len(spot)
new.series = cumprod(1 + reta)
data[[i]]$Adjusted = data[[i]]$Close = spot[n] * new.series / new.series[n]
}
}
if (!auto.assign) {
return(fr)
} else {
return(env)
}
}
getSymbols.TB.test = function() {
filename = 'temp/Futures/CL20_I0B.TXT'
i = 'CL'
data = env()
fr <- read.csv(filename, header = FALSE)
fr <- make.xts(fr[,-1], as.Date(as.character(fr[,1]),'%Y%m%d'))
colnames(fr) <- spl('Open,High,Low,Close,Volume,OpenInterest,DeliveryMonth,Unadjusted')[1:ncol(fr)]
fr$Adjusted = fr$Close
data[[i]] = fr
spot = as.vector(data[[i]]$Unadjusted)
dspot = spot - mlag(spot)
futures = as.vector(data[[i]]$Adjusted)
dfutures = futures - mlag(futures)
index = which(round(dspot - dfutures,4) != 0 )
spot.adjust.roll = spot
spot.adjust.roll[(index-1)] = spot.adjust.roll[index] - dfutures[index]
reta = (mlag(spot.adjust.roll) + futures - mlag(futures)) / mlag(spot.adjust.roll)
reta[1] = 1
n = len(spot)
new.series = cumprod(reta)
Close = spot[n] * new.series / new.series[n]
plot.data = as.xts(list(
Unadjusted = data[[i]]$Unadjusted,
Adjusted = data[[i]]$Adjusted,
Implied.Close = make.xts(Close, index(data[[i]]))
))
png(filename = 'plot_CL_2009.png', width = 600, height = 500, units = 'px', pointsize = 12, bg = 'white')
plota.matplot( scale.one(plot.data['2009']), main='Crude oil, CL - 2009')
dev.off()
}
get.fama.french.data <- function(
name = c('F-F_Research_Data_Factors', 'F-F_Research_Data_Factors'),
periodicity = c('days','weeks', 'months'),
force.download = FALSE,
clean = FALSE,
file.suffix = '_TXT'
)
{
warning('get.fama.french.data is depreciated as of Apr 25, 2016 please use data.ff function instead')
data.ff(name, periodicity, force.download, clean, file.suffix)
}
data.ff <- function(
name = c('F-F_Research_Data_Factors', 'F-F_Research_Data_Factors'),
periodicity = c('days','weeks', 'months'),
force.download = FALSE,
clean = FALSE,
file.suffix = '_TXT'
)
{
map = c(days = '_daily', weeks = '_weekly', months = '')
period = ifna(map[periodicity[1]], periodicity[1])
filename.zip = paste(name[1], period, file.suffix, '.zip', sep='')
filename.txt = paste(name[1], period, '.txt', sep='')
url = paste('http://mba.tuck.dartmouth.edu/pages/faculty/ken.french/ftp/', filename.zip, sep='')
if( !file.exists(filename.zip) || force.download)
download.file(url, filename.zip)
temp.folder = paste(getwd(), 'temp', sep='/')
dir.create(temp.folder, F)
temp.folder = paste(getwd(), '/', 'temp', sep='')
if(clean) shell('del /F /S /Q temp\\*.*', wait = TRUE)
files = unzip(filename.zip, exdir=temp.folder)
if(len(files) == 1) {
filename = paste(temp.folder, '/', filename.txt, sep='')
return( data.ff.internal.one.file(filename) )
}
data = env()
library(stringr)
names = str_trim(str_match(files,'.*/(.*)\\..*')[,2])
for(i in 1:len(files))
data[[ names[i] ]] = data.ff.internal.one.file(files[i])
data
}
data.ff.internal.one.file = function(filename) {
out = readLines(filename)
index = which(nchar(out) == 0)
data.index = grep('^[ 0-9\\.\\+-]+$', out)
temp.index = which(diff(data.index) > 1)
data.index = matrix(data.index[sort(c(1, temp.index, temp.index+1, len(data.index)))], nc=2, byrow=T)
data = list()
for(i in 1:nrow(data.index)) {
start.index = index[which( index > data.index[i,1] ) - 1][1] + 1
if(is.na(start.index)) start.index = index[len(index)] + 1
end.index = data.index[i,1] - 1
n.index = end.index - start.index + 1
name = 'data'
colnames = scan(text = out[start.index], what='', quiet=T)
if(n.index == 2) {
name = trim(out[start.index])
colnames = scan(text = out[end.index], what='', quiet=T)
colnames1 = scan(text = out[end.index+1], what='', quiet=T)
if(len(colnames) > len(colnames1)) {
cindex = which(diff(gregexpr(' ',out[end.index+1])[[1]]) > 1)
cindex = c(1, gregexpr(' ',out[end.index+1])[[1]][(cindex+1)], nchar(out[end.index])+1)
colnames = rep('', len(cindex)-1)
for(j in 2:len(cindex))
colnames[j-1] = substr(out[end.index], cindex[j-1], cindex[j]-1)
colnames = trim(colnames)
colnames = colnames[nchar(colnames) > 0]
}
} else if(n.index > 2) {
name = trim(out[start.index])
colnames0 = scan(text = out[(end.index-1)], what='', quiet=T)
colnames1 = scan(text = out[end.index], what='', quiet=T)
colnames = paste(rep(colnames0, each = len(colnames1) / len(colnames0)), colnames1, sep='.')
}
colnames = gsub('-', '.', colnames)
temp =  matrix(scan(filename, what = double(), quiet=T,
skip = (data.index[i,1]-1),
nlines = (data.index[i,2] - data.index[i,1]+1))
, nc=len(colnames)+1, byrow=T)
date.format = '%Y%m%d'
date.format.add = ''
date.format.n = nchar(paste(temp[1,1]))
if( date.format.n == 6 ) {
date.format.add = '01'
} else if( date.format.n == 4 ) {
date.format.add = '0101'
}
find.name = function(name,data, i=0) if( is.null(data[[name]]) ) name else find.name(paste(name,i+1), data, i+1)
name = find.name(name, data)
data[[name]] = make.xts(temp[,-1], as.Date(paste(temp[,1], date.format.add, sep=''),date.format))
colnames(data[[name]]) = colnames
}
return( data )
}
download.helper <- function(url,download) {
temp.folder = paste(getwd(), 'temp', sep='/')
dir.create(temp.folder, F)
filename = paste0(temp.folder, '/', basename(url))
if(download || !file.exists(filename))
try(download.file(url, filename, mode='wb'), TRUE)
filename
}
getSymbol.CBOE <- function
(
Symbol,
Month,
Year,
download = FALSE
)
{
m.codes = spl('F,G,H,J,K,M,N,Q,U,V,X,Z')
url = paste0("http://cfe.cboe.com/Publish/ScheduledTask/MktData/datahouse/CFE_",
m.codes[Month], substring(Year,3,4), '_', Symbol, '.csv')
filename = download.helper(url, download)
if(file.exists(filename) && file.info(filename)$size > 1)
read.xts(filename, format='%m/%d/%Y')
else
NULL
}
cboe.volatility.term.structure.SPX <- function(make.plot = T) {
url = 'http://www.cboe.com/data/volatilityindexes/volatilityindexes.aspx'
txt = join(readLines(url))
temp.table = extract.table.from.webpage(txt, 'Trade Date', has.header = T)
colnames(temp.table) = gsub(' ','.',trim(tolower(colnames(temp.table))))
temp.table = data.frame(temp.table)
temp.table$trade.date = as.POSIXct(temp.table$trade.date, format="%m/%d/%Y %I:%M:%S %p")
temp.table$expiration.date = as.Date(temp.table$expiration.date, "%d-%b-%y")
temp.table[,3] = as.numeric(as.character(temp.table[,3]))
temp.table[,4] = as.numeric(as.character(temp.table[,4]))
temp.table
if(make.plot) {
plot(temp.table$expiration.date, temp.table$vix, type = 'b',
main=paste('VIX Term Structure, generated ',  max(temp.table$trade.date)),
xlab = 'Expiration Month', ylab='VIX Index Level')
grid()
}
temp.table
}
load.VXX.CBOE <- function() {
index = "::2007-03-25"
fields = spl('Open,High,Low,Close,Settle')
dr <- function(index, date) sum(index>date)
data <- new.env()
futures = list()
i = 1
for(y in 2004:(1+date.year(Sys.Date())))
for(m in 1:12) {
temp = getSymbol.CBOE('VX', m, y)
if(is.null(temp)) next
temp = temp[temp$Settle > 0]
if(nrow(temp)==0) next
if(len(temp[index,1])> 0)
temp[index,fields] = temp[index,fields]/10
label = paste0(y*100+m)
dates = index(temp)
futures[[ i ]] = list()
futures[[ i ]]$data = temp
futures[[ i ]]$label = label
futures[[ i ]]$index = dates
futures[[ i ]]$settle.date = last(dates)
if(i==1)
futures[[ i ]]$dt = len(dates)
else
futures[[ i ]]$dt = dr(dates, futures[[ i-1 ]]$settle.date)
temp$i = i
temp$dt = futures[[ i ]]$dt
temp$dr = (len(dates)-1):0
data[[ label ]] = temp
i = i + 1
}
bt.prep(data, align='keep.all')
data
}
extract.VXX.CBOE <- function(data, field, index, exact.match=T) {
map  = 1:ncol(data$prices)
temp = bt.apply(data, function(x) x[,field])
temp = coredata(temp)
t(apply(temp, 1, function(x) {
if(exact.match) {
pos = map[!is.na(x)][1] - 1
x[(index + pos)]
} else {
pos = map[!is.na(x)][index]
x[pos]
}
}))
}
reconstruct.VXX.CBOE <- function(exact.match=T) {
data = load.VXX.CBOE()
dt = extract.VXX.CBOE(data, 'dt', 1, exact.match)[1,]
dr = extract.VXX.CBOE(data, 'dr', 1, exact.match)[1,]
x  = extract.VXX.CBOE(data, 'Settle', 1:2, exact.match)
w = cbind(dr / dt, (dt - dr) / dt)
val.cur = rowSums(x * mlag(w))
val.yest = rowSums(mlag(x) * mlag(w))
ret = val.cur / val.yest - 1
index = ifna(mlag(dr) == 0, F)
ret[index] = (x[,1] / mlag(x[,2]) - 1)[index]
Close = cumprod(1+ifna(ret,0))
VXX = make.xts(cbind(Close,x,dt,dr,ret), data$dates)
x  = extract.VXX.CBOE(data, 'Settle', 4:7, exact.match)
w = cbind(dr / dt, 1, 1, (dt - dr) / dt)
val.cur = rowSums(x * mlag(w))
val.yest = rowSums(mlag(x) * mlag(w))
ret = val.cur / val.yest - 1
index = ifna(mlag(dr) == 0, F)
ret[index] = (rowSums(x[,-4]) / mlag(rowSums(x[,-1])) - 1)[index]
Close = cumprod(1+ifna(ret,0))
VXZ = make.xts(cbind(Close,x,dt,dr,ret), data$dates)
list(VXX = VXX, VXZ = VXZ)
}
country.code = function
(
force.download = FALSE,
data.filename = 'country.code.Rdata'
)
{
if(!force.download && file.exists(data.filename)) {
load(file=data.filename)
return(temp)
}
url = 'http://www.nationsonline.org/oneworld/country_code_list.htm'
library(curl)
txt = rawToChar(curl_fetch_memory(url)$content)
temp = extract.table.from.webpage(txt, 'Country or Area Name')
temp = trim(temp[,c(2:5)])
colnames(temp)=spl('name,code2,code3,code')
save(temp,file=data.filename)
temp
}
data.ft.search.ticker = function
(
search.field = 'tsx',
sec.type = 'IN'
)
{
url = paste0('http://markets.ft.com/Research/Markets/Company-Search?searchField=', curl_escape(search.field), '&country=&secType=', sec.type)
library(curl)
h = new_handle()
req = curl_fetch_memory(url, h)
if(req$status_code != 200)
warning('error getting data, status_code:', req$status_code, 'for url:', url, 'content:', rawToChar(req$content))
txt = rawToChar(req$content)
temp = gsub(pattern = '<a', replacement = '<td', txt, perl = TRUE)
temp = gsub(pattern = '</a>', replacement = '</td>', temp, perl = TRUE)
temp = extract.table.from.webpage(temp, 'Symbol,Exchange,Country')
colnames(temp)[1:4] = spl('Name,Symbol,Exchange,Country')
temp[,1:4]
}
data.ft.index.members = function
(
ft.symbol = 'INX:IOM',
force.download = FALSE
)
{
data.filename =  paste0(gsub(':','_',ft.symbol),'.Rdata')
if(!force.download && file.exists(data.filename)) {
load(file=data.filename)
return(data)
}
url = paste0('http://markets.ft.com/research/Markets/Tearsheets/Constituents?s=', ft.symbol)
library(curl)
h = new_handle()
txt = rawToChar(curl_fetch_memory(url, h)$content)
temp = gsub(pattern = '<a', replacement = '<td', txt, perl = TRUE)
temp = gsub(pattern = '</a>', replacement = '</td>', temp, perl = TRUE)
temp = extract.table.from.webpage(temp, 'Equities')
token = extract.token(txt,'<div class="wsod-paging-key">','</div>')
nstep = str_match(token,' data-ajax-paging-end-row="([0-9]+)">')[2]
nstep = as.numeric(nstep)
nfound = str_match(token,' data-ajax-paging-total-rows="([0-9]+)">')[2]
nfound = as.numeric(nfound)
data = matrix('',nr=nfound,nc=5)
colnames(data) = spl('Name,Symbol,LastPrice,TodayChange,YearChange')
data[1:nstep,] = temp[,1:5]
token = extract.token(txt,'<div class="wsodHidden">','</div>')
token = spl(token,'<input')
library(stringr)
names = str_match(token,'data-ajax-param="([^"]*)"')[-1,2]
values = str_match(token,'value="([^"]*)"')[-1,2]
settings = as.list(sapply(1:len(names), function(i) { t=c(values[i]); names(t)=names[i]; t}))
settings$ResetPaging = 'false'
h = handle_setopt(h, referer=url)
for(istart in seq(nstep+1,nfound,by=nstep)) {
cat(istart, 'out of', nfound, '\n')
settings$startRow = paste(istart)
handle_setform(h,.list=settings)
url = 'http://markets.ft.com/Research/Remote/UK/Tearsheets/IndexConstituentsPaging'
txt = rawToChar(curl_fetch_memory(url, h)$content)
library(jsonlite)
txt=fromJSON(txt)$html
temp = gsub(pattern = '<a', replacement = '<td', txt, perl = TRUE)
temp = gsub(pattern = '</a>', replacement = '</td>', temp, perl = TRUE)
temp = extract.table.from.webpage(temp, has.header=F)
data[istart:min(istart+nstep-1,nfound),] = temp[,1:5]
}
save(data,file=data.filename)
data
}
get.FOMC.dates <- function
(
force.download = FALSE,
data.filename = 'fomc.Rdata'
)
{
warning('get.FOMC.dates is depreciated as of Apr 25, 2016 please use data.fomc function instead')
data.fomc(force.download, data.filename)
}
data.fomc <- function
(
force.download = FALSE,
data.filename = 'fomc.Rdata'
)
{
url = 'http://www.federalreserve.gov/monetarypolicy/fomccalendars.htm'
txt = join(readLines(url))
library(stringr)
sb = string.buffer()
for(year in 2009:(1 + date.year(Sys.Date()))) {
temp = extract.table.from.webpage(txt, paste(year,'FOMC Meetings'))
if(nrow(temp) == 0) next
temp = tolower(trim(temp[,1:2]))
temp = temp[nchar(temp[,1]) > 0,]
month = temp[,1]
day = temp[,2]
status = paste0(
ifna( str_match(day, '\\(.*\\)'), ''),
ifna( str_match(day, '\\*'), '')
)
status = gsub('\\(','',gsub('\\)','',status))
day = str_replace(day, '(\\(.*?\\))','')
day = str_replace(day, '\\*','')
day = paste(year,
sapply( iif(grepl('/',month), month, paste0(month,'/',month)), spl, '/'),
sapply( iif(grepl('-',day), day, paste0(day,'-',day)), spl, '-')
)
day = matrix(day, nc=2, byrow=T)
for(i in 1:len(status)) add(sb, day[i,1], day[i,2], status[i])
}
data = matrix(scan(what='',text= string(sb),sep=',', quiet=T),nc=3,byrow=T)
close(sb)
sb=NULL
first.year = min(as.numeric(substr(data[,1],1,4)))
recent.data = data
if(!force.download && file.exists(data.filename)) {
load(file=data.filename)
if( last(FOMC$day) == as.Date(last(recent.data[,2]),'%Y %B %d') )
return(FOMC)
}
sb = string.buffer()
for(year in 1936:(first.year-1)) {
cat(year,'\n')
url = paste0('http://www.federalreserve.gov/monetarypolicy/fomchistorical', year, '.htm')
txt = join(readLines(url))
tokens = spl(txt,'<div id="historical">')
for(token in tokens[-1])
add(sb, colnames(extract.table.from.webpage(token, 'year'))[1])
}
data = scan(what='',text= string(sb),sep='\n', quiet=T)
close(sb)
sb=NULL
year = substring(data,nchar(data)-3)
day = tolower(substring(data,1,nchar(data)-4))
status = paste0(
iif(grepl('conference call',day), 'conference call', ''),
iif(grepl('meeting',day), 'meeting', '')
)
day = gsub('conference call', '', gsub('conference calls','',day))
day = gsub('meeting', '', gsub('meetings','',day))
day = gsub(',', '-', gsub('and', '',day))
parse.token = function(year, token) {
parts = trim(spl(token,'-'))
n = len(parts)
if( n > 1 ) {
month = ifna.prev(iif(nchar(parts) > 3,
sapply(parts, function(x) spl(x, ' ')[1]),
NA))
parts = iif(nchar(parts) > 3, parts, paste(month, parts))
}
paste(year, parts[c(1,n)])
}
day = sapply(1:len(day), function(i) parse.token(year[i], day[i]))
all.data = rbind(cbind(t(day), status), recent.data)
FOMC = list(day = as.Date(all.data[,2],'%Y %B %d'), start.day = as.Date(all.data[,1],'%Y %B %d'), status=all.data[,3])
save(FOMC,file=data.filename)
FOMC
}
edgar.info <- function(ticker)
{
url = paste0('http://www.sec.gov/cgi-bin/browse-edgar?CIK=', ticker, '&Find=Search&owner=exclude&action=getcompany')
txt = join(readLines(url))
out = list()
temp = extract.table.from.webpage(txt, 'seriesDiv,Filings', has.header = T)
out$fillings= clean.table(temp)
temp = extract.token(txt, 'contentDiv,mailer,Mailing Address','</div>')
out$mailing = t(clean.table(extract.table.from.webpage(temp, has.header=F, end.marker='</span>')))
colnames(out$mailing) = 'Mailing Address'
temp = extract.token(txt, 'contentDiv,mailer,Business Address','</div>')
out$business = t(clean.table(extract.table.from.webpage(temp, has.header=F, end.marker='</span>')))
colnames(out$business) = 'Business Address'
temp = extract.token(txt, 'contentDiv,companyInfo,>','</div>')
temp = gsub('\\|', '</span>', replace.token(temp, '<br','>','</span>'))
temp = clean.table(extract.table.from.webpage(temp, has.header=F, end.marker='</span>'))
out$company = t(temp)
colnames(out$company) = 'Company Info'
out$sic = trim(spl(spl(temp[grep('SIC', temp)],':')[2],'-'))
return(out)
}
zacks.info <- function(ticker = 'IBM')
{
url = paste0('http://www.zacks.com/stock/research/', ticker, '/earnings-announcements')
txt = join(readLines(url))
out = list()
require(jsonlite)
for(i in spl('earnings,webcasts,revisions,splits,dividends,guidance')) {
data = extract.token(txt,paste0('<script>,window.app_data_', i, ',=,"data"'),'</script>')
data = fromJSON(paste('{"data"', data))
out[[i]] = data$data
}
out
}
quantumonline.info <- function
(
id,
type=c(
'cusip',
'symbol',
'sname'
)
)
{
url = paste0('http://www.quantumonline.com/search.cfm?tickersymbol=', id, '&sopt=', type[1], '&1.0.1=Search')
txt = join(readLines(url))
out = list()
out$main = extract.table.from.webpage(gsub('&nbsp;', ',', txt), "Company's Online Profile", has.header = F)
out$address = extract.table.from.webpage( txt, 'Address:', has.header = F)
return(out)
}
hist.quotes.url <- function
(
ticker = 'IBM',
from = '1900-01-01',
to = Sys.Date(),
src = spl('yahoo,google,quotemedia')
)
{
if(class(from) != 'Date') from = as.Date(from, '%Y-%m-%d')
if(class(to) != 'Date') to = as.Date(to, '%Y-%m-%d')
switch(src[1],
yahoo = paste('http://ichart.finance.yahoo.com/table.csv?',
's=', ticker,
'&a=', sprintf('%.2d', date.month(from) - 1),
format(from, '&b=%d&c=%Y'),
'&d=', sprintf('%.2d', date.month(to) - 1),
format(to, '&e=%d&f=%Y'),
'&g=d&q=q&y=0&z=file&x=.csv',
sep=''),
google = paste('http://finance.google.com/finance/historical?',
'q=', ticker,
'&startdate=', format(from, '%b+%d+%Y'),
'&enddate=', format(to, '%b+%d+%Y'),
'&output=csv',
sep=''),
quotemedia = paste('http://app.quotemedia.com/quotetools/getHistoryDownload.csv?webmasterId=501&',
'symbol=', ticker,
'&startMonth=', sprintf('%.2d', date.month(from) - 1),
format(from, '&startDay=%d&startYear=%Y'),
'&endMonth=', sprintf('%.2d', date.month(to) - 1),
format(to, '&endDay=%d&endYear=%Y'),
'&isRanged=true',
sep=''),
''
)
}
data.clean <- function
(
data,
min.ratio = 2.5,
min.obs = 3*252,
iqr.mult = 20
)
{
data$symbolnames = iif(is.null(data$symbolnames), ls(data), data$symbolnames)
if(min.obs > 0) {
index = names(which(sapply(data$symbolnames, function(x) as.numeric(count(Cl(data[[x]])))) < min.obs))
if (len(index) > 0) {
cat('Removing', index, 'have less than', min.obs, 'observations','\n')
rm(list=index, envir=data)
data$symbolnames = setdiff(data$symbolnames, index)
}
}
for(ticker in data$symbolnames)
data[[ticker]] = data.clean.helper(data[[ticker]], ticker, min.ratio, iqr.mult)
}
data.clean.helper <- function
(
data,
ticker,
min.ratio = 2.5,
iqr.mult = 20
)
{
data = data[Cl(data) > 0 & Ad(data) > 0]
nperiods = nrow(data)
price = Ad(data)
ratio = as.vector((price)/mlag(price))
index = which(ratio > min.ratio)
if(len(index) > 0)
for(i in index) {
cat('Abnormal price found for', ticker, format(index(data)[i],'%d-%b-%Y'),'Ratio :', round(ratio[i],1),'\n')
for(name in find.names('Open,Close,High,Low,Adjusted', data))
data[i:nperiods,name] = data[i:nperiods,name] / ratio[i]
}
price = Ad(data)
ret = as.vector((price)/mlag(price)) - 1
threshold = iqr.mult * IQR(ret, na.rm=T)
index = which(ret > threshold | ret < -threshold)
if(len(index) > 0)
for(i in index) {
cat('Abnormal price found for', ticker, format(index(data)[i],'%d-%b-%Y'),'based on IQR, Ratio :', round(ratio[i],1),'\n')
for(name in find.names('Open,Close,High,Low,Adjusted', data))
data[i:nperiods,name] = data[i:nperiods,name] / ratio[i]
}
price = Ad(data)
ratio = as.vector(mlag(price)/(price))
index = which(ratio > min.ratio)
if(len(index) > 0)
for(i in index) {
cat('Abnormal price found for', ticker, format(index(data)[i],'%d-%b-%Y'),'Inverse Ratio :', round(ratio[i],1),'\n')
for(name in find.names('Open,Close,High,Low,Adjusted', data))
data[i:nperiods,name] = data[i:nperiods,name] * ratio[i]
}
data
}
make.data.proxy <- function() {
load.packages('quantmod')
raw.data = env()
filename = 'data/TR_CC-CRB'
if(file.exists(filename)) {
temp = extract.table.from.webpage( join(readLines(filename)), 'EODValue' )
temp = join( apply(temp, 1, join, ','), '\n' )
raw.data$CRB = make.stock.xts( read.xts(temp, format='%m/%d/%y' ) )
}
filename = 'data/TB3M.Rdata'
if(!file.exists(filename)) {
TB3M = quantmod::getSymbols('DTB3', src='FRED', auto.assign = FALSE)
save(TB3M, file=filename)
}
load(file=filename)
TB3M[] = ifna.prev(TB3M)
raw.data$TB3M = make.stock.xts(processTBill(TB3M, timetomaturity = 1/4, 261))
filename = 'data/TB3Y.Rdata'
if(!file.exists(filename)) {
TB3Y = quantmod::getSymbols('DGS3', src='FRED', auto.assign = FALSE)
save(TB3Y, file=filename)
}
load(file=filename)
TB3Y[] = ifna.prev(TB3Y)
raw.data$TB3Y = make.stock.xts(processTBill(TB3Y, timetomaturity = 3, 261))
filename = 'data/TB10Y.Rdata'
if(!file.exists(filename)) {
TB10Y = quantmod::getSymbols('DGS10', src='FRED', auto.assign = FALSE)
save(TB10Y, file=filename)
}
load(file=filename)
TB10Y[] = ifna.prev(TB10Y)
raw.data$TB10Y = make.stock.xts(processTBill(TB10Y, timetomaturity = 10, 261))
filename = 'data/TB20Y.Rdata'
if(!file.exists(filename)) {
TB20Y = quantmod::getSymbols('GS20', src='FRED', auto.assign = FALSE)
save(TB20Y, file=filename)
}
load(file=filename)
TB20Y[] = ifna.prev(TB20Y)
raw.data$TB20Y = make.stock.xts(processTBill(TB20Y, timetomaturity = 20, 12))
filename = 'data/GOLD.Rdata'
if(!file.exists(filename)) {
GOLD = bundes.bank.data.gold()
save(GOLD, file=filename)
}
load(file=filename)
raw.data$GOLD = make.stock.xts(GOLD)
filename = 'data/NAREIT.xls'
if(!file.exists(filename)) {
url = 'http://returns.reit.com/returns/MonthlyHistoricalReturns.xls'
download.file(url, filename,  mode = 'wb')
}
load.packages('readxl')
temp = read_excel(filename, sheet='Index Data', skip=7)
NAREIT = make.xts(temp$Index, as.Date(temp$Date))
raw.data$NAREIT = make.stock.xts(NAREIT)
tickers = '
COM = DBC;GSG + CRB
RExUS = [RWX] + VNQ + VGSIX
RE = [RWX] + VNQ + VGSIX
RE.US = [ICF] + VGSIX
EMER.EQ = [EEM] + VEIEX
EMER.FI = [EMB] + PREMX
GOLD = [GLD] + GOLD,
US.CASH = [BIL] + TB3M,
SHY + TB3Y,
US.HY = [HYG] + VWEHX
US.BOND = [AGG] + VBMFX
INTL.BOND = [BWX] + BEGBX
JAPAN.EQ = [EWJ] + FJPNX
EUROPE.EQ = [IEV] + FIEUX
US.SMCAP = IWM;VB + NAESX
TECH.EQ = [QQQ] + ^NDX
US.EQ = [VTI] + VTSMX + VFINX
US.MID = [VO] + VIMSX
EAFE = [EFA] + VDMIX + VGTSX
MID.TR = [IEF] + VFITX
CORP.FI = [LQD] + VWESX
TIPS = [TIP] + VIPSX + LSGSX
LONG.TR = [TLT] + VUSTX
'
data.proxy = env()
getSymbols.extra(tickers, src = 'yahoo', from = '1970-01-01', env = data.proxy, raw.data = raw.data, auto.assign = T)
data.proxy.raw = raw.data
save(data.proxy.raw, file='data/data.proxy.raw.Rdata',compress='gzip')
save(data.proxy, file='data/data.proxy.Rdata',compress='gzip')
}
load.aqr.data = function
(
data.set = 'betting-against-beta-equity-factors', #'time-series-momentum-factors'
frequency = c('monthly','daily'),
sheet = 1,
force.download = F,
last.col2extract = 'Global'
)
{
warning('load.aqr.data is depreciated as of Apr 25, 2016 please use data.aqr function instead')
data.aqr(data.set, frequency, sheet, force.download, last.col2extract)
}
data.aqr = function
(
data.set = 'betting-against-beta-equity-factors', #'time-series-momentum-factors'
frequency = c('monthly','daily'),
sheet = 1,
force.download = F,
last.col2extract = 'Global'
)
{
data.folder = paste(getwd(), 'aqr.data', sep='/')
url = paste0('http://www.aqr.com/library/data-sets/', data.set, '-', frequency[1], '/data')
filename = file.path(data.folder, paste0(data.set, '-', frequency[1],'.xlsx'))
if( !file.exists(filename) || force.download) {
dir.create(data.folder, F)
download.file(url, filename,  mode = 'wb')
}
require(readxl)
data = read_excel(filename, sheet=sheet)
skip = which(data[,1]=='DATE')
data = read_excel(filename, sheet=sheet,skip=skip)
if( is.character(last.col2extract) ) last.col2extract = which(colnames(data)==last.col2extract)-1
data = data[!is.na(data[,1]), 1:last.col2extract]
data = data[rowSums(!is.na(data[,-1,drop=F])) > 0,]
make.xts(data[,-1], as.Date(data[,1]))
}
load.csi.security.master = function(force.download = F) {
data.folder = paste(getwd(), 'csi.data', sep='/')
url = 'http://www.csidata.com/factsheets.php?type=commodity&format=csv'
filename = file.path(data.folder, 'commodityfactsheet.csv')
if( !file.exists(filename) || force.download) {
dir.create(data.folder, F)
download.file(url, filename,  mode = 'wb')
}
read.csv(filename)
}
fred.fx.symbol = function() {
url = 'https://research.stlouisfed.org/fred2/release/tables?rid=17&eid=23340'
txt = join(readLines(url))
temp = gsub(pattern = 'series', replacement = '<td>', txt, perl = TRUE)
temp = gsub(pattern = 'target', replacement = '</td><', temp, perl = TRUE)
temp = extract.table.from.webpage(temp, 'Country', has.header = F)
data = gsub('/','',gsub('"','',trim(temp[,c(2,3,7)])))
colnames(data) = spl('symbol,name,description')
data[,'description']
keep.index = !is.na(data[,'description']) & nchar(data[,'description']) > 0
data = data.frame(data[keep.index,])
index = grep('index',data[,'description'],T)
list(fx = data[-index,], index = data[index,])
}
fxhistoricaldata.fx.symbol = function() {
url = 'http://www.fxhistoricaldata.com/'
txt = join(readLines(url))
temp = gsub(pattern = '<ul>', replacement = '<table>', txt, perl = TRUE)
temp = gsub(pattern = '</ul>', replacement = '</table>', temp, perl = TRUE)
temp = gsub(pattern = '<li>', replacement = '<td>', temp, perl = TRUE)
temp = gsub(pattern = '</li>', replacement = '</td>', temp, perl = TRUE)
temp = gsub(pattern = '<!--.*?-->', replacement = '', temp, perl = TRUE)
temp = extract.table.from.webpage(temp, 'EURUSD', has.header = F)
as.character(temp)
}
consecutive.changes <- function
(
data,
positive=T
)
{
if(positive) dir = diff(data) > 0 else dir = diff(data) < 0
temp = cumsum(iif(dir, 1, 0))
temp - ifna.prev(iif(dir, NA, coredata(temp)))
}
factor.avgcor <- function(data, next.month.ret, name) {
load.packages('abind')
temp = abind(data, along = 3)
temp = abind(next.month.ret, temp, along = 3)
dimnames(temp)[[3]][1] = 'Ret'
temp = t(compute.avgcor(temp, 'spearman')[,-1])
temp[] = plota.format(100 * temp, 0, '', '%')
plot.table(temp, smain=paste(name,'Correlation',sep=' \n '))
}
compute.avgcor <- function
(
data,
method = c('pearson', 'kendall', 'spearman')
)
{
nr = dim(data)[1]
nc = dim(data)[3]
corm = matrix(NA,nc,nc)
colnames(corm) = rownames(corm) = dimnames(data)[[3]]
for( i in 1:(nc-1) ) {
for( j in (i+1):nc ) {
corm[i,j] = mean( as.double( sapply(1:nr, function(t)
try(cor(data[t,,i], data[t,,j], use = 'complete.obs', method[1]),TRUE)
)), na.rm=T)
}
}
return(corm)
}
cap.weighted.mean <- function
(
data,
capitalization
)
{
capitalization = capitalization * (!is.na(data))
weight = capitalization / rowSums(capitalization,na.rm=T)
rowSums(data * weight,na.rm=T)
}
sector.mean <- function
(
data,
sectors
)
{
out = data * NA
for(sector in levels(sectors)) {
index = (sector == sectors)
out[,index] = ifna(apply(data[,index, drop=F], 1, mean, na.rm=T),NA)
}
return(out)
}
compute.quantiles <- function
(
data,
next.month.ret,
smain='',
n.quantiles=5,
plot=T
)
{
n = ncol(data)
nperiods = nrow(data)
data = coredata(ifna(data,NA))
next.month.ret = coredata(ifna(next.month.ret,NA))
temp = matrix(NA, nperiods, n.quantiles)
hist.factor.quantiles = hist.ret.quantiles = temp
temp = matrix(NA, nperiods, n)
quantiles = weights = ranking = temp
index = which(rowSums(!is.na(data)) >= n.quantiles)
for(t in index) {
factor = data[t,]
ret = next.month.ret[t,]
ranking[t,] = rank(factor, na.last = 'keep','first')
t.ranking = ceiling(n.quantiles * ranking[t,] / count(factor))
quantiles[t,] = t.ranking
weights[t,] = 1/tapply(rep(1,n), t.ranking, sum)[t.ranking]
hist.factor.quantiles[t,] = tapply(factor, t.ranking, mean)
hist.ret.quantiles[t,] = tapply(ret, t.ranking, mean)
}
if(plot) {
par(mar=c(4,4,2,1))
temp = 100*apply(hist.ret.quantiles,2,mean,na.rm=T)
barplot(temp, names.arg=paste(1:n.quantiles), ylab='%',
main=paste(smain, ', spread =',round(temp[n.quantiles]-temp[1],2), '%'))
}
return(list(quantiles=quantiles, weights=weights, ranking=ranking,
hist.factor.quantiles = hist.factor.quantiles, hist.ret.quantiles = hist.ret.quantiles))
}
add.avg.factor <- function
(
data
)
{
temp = abind(data, along = 3)
data$AVG = data[[1]]
data$AVG[] = ifna(apply(temp, c(1,2), mean, na.rm=T),NA)
return(data)
}
normalize.mkval <- function
(
data,
MKVAL
)
{
for(i in names(data)) {
data[[i]] = (data[[i]] - cap.weighted.mean(data[[i]], MKVAL)) /
apply(data[[i]], 1, sd, na.rm=T)
}
return(data)
}
normal.transform <- function(data)
{
rk=rank(data, na.last='keep', ties.method = 'first')
n = count(data)
x = qnorm((1:n) / (n+1))
return(x[rk])
}
normalize.normal <- function
(
data
)
{
for(i in names(data)) {
data[[i]][] = t(apply(data[[i]], 1, normal.transform))
}
return(data)
}
plot.quantiles <- function
(
data,
next.month.ret,
smain=''
)
{
layout(matrix(1:(2*ceiling(len(data)/2)), nc=2))
sapply(1:len(data), function(i)
compute.quantiles(data[[i]], next.month.ret, paste(names(data)[i],smain))
)
}
plot.bt.quantiles <- function
(
factors,
next.month.ret,
smain='',
data
)
{
out = compute.quantiles(factors, next.month.ret, plot=F)
prices = data$prices
prices = bt.apply.matrix(prices, function(x) ifna.prev(x))
month.ends = endpoints(prices, 'months')
models = list()
for(i in 1:5) {
data$weight[] = NA
data$weight[month.ends,] = iif(out$quantiles == i, out$weights, 0)
capital = 100000
data$weight[] = (capital / prices) * (data$weight)
models[[paste('Q',i,sep='')]] = bt.run(data, type='share', capital=capital)
}
data$weight[] = NA
data$weight[month.ends,] = iif(out$quantiles == 5, out$weights,
iif(out$quantiles == 1, -out$weights, 0))
capital = 100000
data$weight[] = (capital / prices) * (data$weight)
models$Q5_Q1 = bt.run(data, type='share', capital=capital)
plotbt(models, plotX = T, log = 'y', LeftMargin = 3, main=smain)
mtext('Cumulative Performance', side = 2, line = 1)
}
plot.factors <- function
(
data,
name,
next.month.ret
)
{
x = as.vector(t(data))
y = as.vector(t(next.month.ret))
x = ifna(x,NA)
y = ifna(y,NA)
index = !is.na(x) & !is.na(y)
x = x[index]
y = y[index]
cor.p = round(100*cor(x, y, use = 'complete.obs', method = 'pearson'),1)
cor.s = round(100*cor(x, y, use = 'complete.obs', method = 'spearman'),1)
layout(1:2)
plot(x, pch=20)
par(mar=c(4,4,2,1))
plot(x, y, pch=20, xlab=name, ylab='Next Month Return')
abline(lm(y ~ x), col='blue', lwd=2)
plota.legend(paste('Pearson =',cor.p,',Spearman =', cor.s))
}
fund.data <- function
(
Symbol,
n=10,
mode=c('quarterly','annual'),
max.attempts=5
)
{
all.data = c()
option.value = -1
start_date = spl('istart_date,start_date')
names(start_date) = spl('quarterly,annual')
repeat {
if(option.value >= 0) {
url = paste('http://uk.advfn.com/p.php?pid=financials&symbol=', Symbol, '&btn=', mode[1], '_reports&', start_date[mode[1]], '=', option.value, sep = '')
} else {
url = paste('http://uk.advfn.com/p.php?pid=financials&symbol=', Symbol, '&btn=', mode[1], '_reports', sep = '')
}
cat('Downloading', url, '\n')
for(iattempt in 1:max.attempts) {
flag = T
tryCatch({
txt = join(readLines(url))
}, interrupt = function(ex) {
flag <<-  F
Sys.sleep(0.1)
}, error = function(ex) {
flag <<-  F
Sys.sleep(0.1)
}, finally = {
if(flag) break
})
}
if( len(grep('INDICATORS', txt, ignore.case = T)) == 0 ) {
cat('No Data Found for', Symbol, '\n')
return(all.data)
}
pos = regexpr(pattern = '<title>(.*?)</title>', txt, ignore.case = TRUE, perl = TRUE)
if(len(pos) == 1)
title = substr(txt, attr(pos, 'capture.start'), attr(pos, 'capture.start') + attr(pos, 'capture.length') - 1)
data = extract.table.from.webpage(txt, 'INDICATORS', has.header = T)
colnames(data) = data[1,]
rownames(data) = data[,1]
data = data[,-1,drop=F]
add.index = which( is.na(match( colnames(data), colnames(all.data) )) )
all.data = cbind(data[,add.index,drop=F], all.data)
if(ncol(all.data) >= n) break
if(option.value == 0)  break
temp = gsub(pattern = '<option', replacement = '<tr>', txt, perl = TRUE)
temp = gsub(pattern = '</option>', replacement = '</tr>', temp, perl = TRUE)
temp = extract.table.from.webpage(temp, 'All amounts', has.header = T)
temp = apply(temp,1,join)
index.selected = grep('selected', temp)
option.value = 0
if(	len(index.selected) )
option.value = as.double( gsub('.*value=\'([0-9]*).*', '\\1', temp[index.selected]) )
if(option.value > 0) {
option.value = option.value - 5
option.value = max(0, option.value)
} else {
break
}
}
all.data = all.data[, colSums(nchar(trim(all.data))) > 0, drop=F]
all.data = rbind(all.data, title)
rownames(all.data)[nrow(all.data)] = 'HTMLTITLEtext'
if( ncol(all.data) > n ) {
return(all.data[,(ncol(all.data)-n+1):ncol(all.data), drop=F])
} else {
return(all.data)
}
}
date.fund.data <- function(data)
{
quarter.end.date = as.Date(paste(data[1,], '/1', sep=''), '%Y/%m/%d')
quarterly.indicator = data['quarterly indicator',]
date.preliminary.data.loaded = as.Date(data['date preliminary data loaded',], '%Y-%m-%d') + 1
months = seq(quarter.end.date[1], tail(quarter.end.date,1)+365, by='1 month')
index = match(quarter.end.date, months)
quarter.end.date = months[ iif(quarterly.indicator == '4', index+3, index+2) + 1 ] - 1
fund.date = date.preliminary.data.loaded
fund.date[is.na(fund.date)] = quarter.end.date[is.na(fund.date)]
return(fund.date)
}
get.fund.data.index <- function
(
label,
fund,
silent = T
)
{
names = rownames(fund)
index = grep(label, names, ignore.case = T)
if( len(index) == 0 ) {
labels = spl(label,' ')
n = len(labels)
temp.count = rep(0,nrow(fund))
for(ilabel in labels) {
index = grep(ilabel, rownames(fund), ignore.case = T)
if(len(index)>0) temp.count[index] = temp.count[index]+1
}
index = which(temp.count == n)
if( !silent ) cat('Exact label not found, trying partial match\n')
}
if( len(index) > 0 ) {
if( len(index) > 1 ) {
if( !silent ) cat('Possible Matches', rownames(fund)[index], '\n', sep=' | ')
index = index[ which.min(nchar(names[index]) - nchar(label)) ]
}
if( !silent ) cat('Match =', rownames(fund)[index], '\n')
index[1]
} else {
if( !silent ) cat('No Match Found for', label, '\n')
c()
}
}
get.fund.data <- function
(
label,
fund,
fund.date,
is.12m.rolling=F,
cash.flow=F
)
{
index = get.fund.data.index(label, fund)
if( len(index) == 0 ) return(as.xts(rep(NA,len(fund.date)), fund.date))
temp.q = as.double(gsub(',', '', fund[index,]))
temp.q = ifna(temp.q, 0)
if(cash.flow) {
quarterly.indicator = fund['quarterly indicator',]
temp.q = iif(quarterly.indicator == '1', temp.q, temp.q - mlag(temp.q))
}
temp.q = as.xts(temp.q, fund.date)
iif(is.12m.rolling, runSum(temp.q, 4), temp.q)
}
iline = function(type=c('v','h','cross'), col='red', remove.col='white',stop.key = 'q') {
type = tolower(substr(type[1],1,1))
keydown <- function(key) {
if (key == stop.key) return(invisible(1))
NULL
}
prev.x = NULL
prev.y = NULL
mouseup <- function(buttons, x, y) {
par(xpd=NA)
if(type == 'v' || type == 'c')
if(!is.null(prev.x)) abline(v=prev.x, col=remove.col)
if(type == 'h' || type == 'c')
if(!is.null(prev.y)) abline(h=prev.y, col=remove.col)
prev.x <<- grconvertX(x, "ndc", "user")
prev.y <<- grconvertY(y, "ndc", "user")
if(type == 'v' || type == 'c')
abline(v=prev.x, col=col)
if(type == 'h' || type == 'c')
abline(h=prev.y, col=col)
par(xpd=FALSE)
NULL
}
getGraphicsEvent(prompt = "Click to plot v/h/cross line, hit q to quit",
onMouseDown = mouseup,onKeybd = keydown)
}
visualize.system.parameter.optimization = function(result) {
load.packages('rpanel,tkrplot,MASS')
draw = function(panel) {
par(bg='white')
col = cols
index = colSums(result.t >= panel$min & result.t <= panel$max) == n
col[ index ] = 2
parcoord(result, var.label =T, col=col, lwd=col)
panel
}
redraw = function(panel) {
rp.tkrreplot(panel, tkrp)
panel
}
max.val = apply(result,2,max)
min.val = apply(result,2,min)
max.val = max.val + 0.1 * abs(max.val)
min.val = min.val - 0.1 * abs(min.val)
n = ncol(result)
cols = rep(1, nrow(result))
result.t = t(result)
colnames(result) = paste0(colnames(result), '\n(', 1:n,')')
panel  = rp.control(title = 'Parallel Coordinates', max=max.val, min=min.val)
rp.slider(panel, max, max.val, min.val, horizontal=F, showvalue = T, action = redraw,initval=max.val)
rp.slider(panel, min, max.val, min.val, horizontal=F, showvalue = T, action = redraw,initval=min.val)
rp.tkrplot(panel, tkrp, draw)
}
min.corr.paper.examples <- function()
{
load.packages('quantmod')
data <- new.env()
getSymbols.TB(env = data, auto.assign = T, download = T)
bt.prep(data, align='remove.na', dates='1990::')
save(data,file='FuturesForex.Rdata')
tickers = spl('SPY,QQQ,EEM,IWM,EFA,TLT,IYR,GLD')
data <- new.env()
getSymbols(tickers, src = 'yahoo', from = '1980-01-01', env = data, auto.assign = T)
for(i in ls(data)) data[[i]] = adjustOHLC(data[[i]], use.Adjusted=T)
bt.prep(data, align='keep.all', dates='2002:08::')
save(data,file='ETF.Rdata')
load.packages('quantmod,quadprog')
tickers = spl('AA,AXP,BA,CAT,DD,DIS,GE,IBM,IP,JNJ,JPM,KO,MCD,MMM,MO,MRK,MSFT')
data <- new.env()
getSymbols(tickers, src = 'yahoo', from = '1980-01-01', env = data, auto.assign = T)
for(i in ls(data)) data[[i]] = adjustOHLC(data[[i]], use.Adjusted=T)
bt.prep(data, align='keep.all', dates='1980::')
prices = coredata(data$prices)
prices[is.na(prices)] = mlag(prices)[is.na(prices)]
prices[is.na(prices)] = mlag(prices)[is.na(prices)]
data$prices[] = prices
save(data,file='Dow.Engle.Rdata')
load.packages('quantmod,quadprog')
tickers = spl('VTI,IEV,EEM,EWJ,AGG,GSG,GLD,ICF')
data <- new.env()
getSymbols(tickers, src = 'yahoo', from = '1980-01-01', env = data, auto.assign = T)
for(i in ls(data)) data[[i]] = adjustOHLC(data[[i]], use.Adjusted=T)
bt.prep(data, align='keep.all', dates='2003:10::')
save(data,file='ETF2.Rdata')
load.packages('quantmod,quadprog')
tickers = spl('ATVI,ADBE,ALTR,AMZN,AMGN,APOL,AAPL,AMAT,ADSK,ADP,BBBY,BIIB,BMC,BRCM,CHRW,CA,CELG,CERN,CHKP,CTAS,CSCO,CTXS,CTSH,CMCSA,COST,DELL,XRAY,DISH,EBAY,EA,EXPD,ESRX,FAST,FISV,FLEX,FLIR,FWLT,GILD,HSIC,HOLX,INFY,INTC,INTU,JBHT,KLAC,LRCX,LIFE,LLTC,LOGI,MAT,MXIM,MCHP,MSFT,MYL,NTAP,NWSA,NVDA,ORLY,ORCL,PCAR,PDCO,PAYX,PCLN,QGEN,QCOM,BBRY,ROST,SNDK,SIAL,SPLS,SBUX,SRCL,SYMC,TEVA,URBN,VRSN,VRTX,VOD,XLNX,YHOO')
data <- new.env()
for(i in tickers) {
try(getSymbols(i, src = 'yahoo', from = '1980-01-01', env = data, auto.assign = T), TRUE)
data[[i]] = adjustOHLC(data[[i]], use.Adjusted=T)
}
bt.prep(data, align='keep.all', dates='1995::')
prices = coredata(data$prices)
prices[is.na(prices)] = mlag(prices)[is.na(prices)]
prices[is.na(prices)] = mlag(prices)[is.na(prices)]
data$prices[] = prices
save(data,file='nasdaq.100.Rdata')
names = spl('ETF,FuturesForex,Dow.Engle,ETF2,nasdaq.100')
lookback.len = 60
periodicitys = spl('weeks,months')
periodicity = periodicitys[1]
prefix = paste(substr(periodicity,1,1), '.', sep='')
for(name in names) {
load(file = paste(name, '.Rdata', sep=''))
obj = portfolio.allocation.helper(data$prices, periodicity, lookback.len = lookback.len, prefix = prefix,
min.risk.fns = 'min.corr.portfolio,min.corr2.portfolio,max.div.portfolio,min.var.portfolio,risk.parity.portfolio(),equal.weight.portfolio',
custom.stats.fn = 'portfolio.allocation.custom.stats')
save(obj, file=paste(name, lookback.len, periodicity, '.bt', '.Rdata', sep=''))
}
for(name in names) {
load(file=paste(name, '.Rdata', sep=''))
custom.input.report.helper(paste('report.', name, sep=''), data)
load(file=paste(name, lookback.len, periodicity, '.bt', '.Rdata', sep=''))
custom.report.helper(paste('report.', name, lookback.len, periodicity, sep=''),
create.strategies(obj, data))
}
names = spl('FuturesForex')
for(name in names) {
load(file=paste(name, '.Rdata', sep=''))
load(file=paste(name, lookback.len, periodicity, '.bt', '.Rdata', sep=''))
leverage = c(5, 4, 15, 20, 3, 1)
custom.report.helper(paste('report.leverage.', name, lookback.len, periodicity, sep=''),
create.strategies(obj, data, leverage))
}
}
custom.input.report.helper <- function(filename, data) {
filename.pdf = paste(filename, '.pdf', sep='')
filename.csv = paste(filename, '.csv', sep='')
pdf(file = filename.pdf, width=8.5, height=11)
layout(1:2)
asset.models = list()
for(i in data$symbolnames) {
data$weight[] = NA
data$weight[,i] = 1
asset.models[[i]] = bt.run(data, silent=T)
}
asset.summary = plotbt.strategy.sidebyside(asset.models, return.table=T)
ret.log = bt.apply.matrix(data$prices, ROC, type='continuous')
temp = cor(ret.log, use='complete.obs', method='pearson')
temp[] = plota.format(100 * temp, 0, '', '%')
plot.table(temp, smain='Correlation', highlight = TRUE, colorbar = TRUE)
layout(matrix(1:4,2,2))
if( is.null(data$symbol.groups) ) {
index = order(data$symbolnames)
for(i in data$symbolnames[index])
plota(data[[i]], type='l', cex.main=0.7,main= i)
} else {
index = order(data$symbol.groups)
for(i in data$symbolnames[index])
plota(data[[i]], type='l', cex.main=0.7, main= paste(i, data$symbol.groups[i], data$symbol.descriptions.print[i], sep=' / ') )
asset.summary = rbind(data$symbol.groups, data$symbol.descriptions.print, asset.summary)
}
dev.off()
load.packages('abind')
write.csv(asset.summary, filename.csv)
cat('\n\n', file=filename.csv, append=TRUE)
write.table(temp, sep=',',  row.names = , col.names = NA,
file=filename.csv, append=TRUE)
}
custom.report.helper <- function(filename, obj) {
filename.pdf = paste(filename, '.pdf', sep='')
filename.csv = paste(filename, '.csv', sep='')
models = obj$models
pdf(file = filename.pdf, width=8.5, height=11)
plotbt(models, plotX = T, log = 'y', LeftMargin = 3)
mtext('Cumulative Performance', side = 2, line = 1)
out = plotbt.strategy.sidebyside(models, perfromance.fn = 'custom.returns.kpi', return.table=T)
cdi = custom.composite.diversification.indicator(obj)
out = rbind(colMeans(cdi, na.rm=T), out)
rownames(out)[1] = 'Composite Diversification Indicator(CDI)'
y = 100 * sapply(models, compute.turnover, data)
out = rbind(y, out)
rownames(out)[1] = 'Portfolio Turnover'
performance.barchart.helper(out, 'Sharpe,Cagr,RC Gini,RC Herfindahl,Volatility,Portfolio Turnover,Composite Diversification Indicator(CDI)', c(T,T,F,F,F,F,T))
custom.summary.positions(obj$weights)
custom.period.chart(models)
layout(1:len(models))
for(m in names(models)) {
plotbt.transition.map(models[[m]]$weight, name=m)
legend('topright', legend = m, bty = 'n')
}
dates = index(models[[1]]$weight)[obj$period.ends]
layout(1:len(models))
for(m in names(models)) {
plotbt.transition.map(make.xts(obj$risk.contributions[[m]], dates),
name=paste('Risk Contributions',m))
legend('topright', legend = m, bty = 'n')
}
dev.off()
load.packages('abind')
write.csv(out, filename.csv)
cat('\n\n', file=filename.csv, append=TRUE)
out = abind(lapply(models, function(m) m$equity))
colnames(out) = names(models)
write.xts(make.xts(out, index(models[[1]]$equity)), filename.csv, append=TRUE)
}
custom.composite.diversification.indicator <- function
(
obj,
avg.flag = T,
avg.len = 10,
plot.main = T,
plot.table = T
)
{
cdi = 0.5 * obj$risk.gini + 0.5 * obj$degree.diversification
if(avg.flag) cdi = bt.apply.matrix(cdi, EMA, avg.len)
if(plot.main) {
avg.name = iif(avg.flag, paste(avg.len,'period EMA') , '')
layout(1:3)
out = obj$degree.diversification
if(avg.flag) out = bt.apply.matrix(out, EMA, avg.len)
plota.matplot(out, cex.main = 1,
main=paste('D = 1 - Portfolio Risk/Weighted Average of asset vols in the portfolio', avg.name))
out = obj$risk.gini
if(avg.flag) out = bt.apply.matrix(out, EMA, avg.len)
plota.matplot(out, cex.main = 1,
main=paste('1 - Gini(Risk Contributions)', avg.name))
plota.matplot(cdi, cex.main = 1,
main=paste('Composite Diversification Indicator (CDI) = 50/50 Gini/D', avg.name))
}
if(plot.table) {
weights = seq(0,1,0.1)
temp = matrix(NA, nc=len(weights), nr=ncol(cdi))
colnames(temp) = paste(weights)
rownames(temp) = colnames(cdi)
for(j in 1:len(weights)) {
i = weights[j]
temp[,j] = rank(-colMeans((1-i) * obj$risk.gini + i * obj$degree.diversification, na.rm=T))
}
temp = cbind(temp, round(rowMeans(temp),1))
colnames(temp)[ncol(temp)] = 'AVG'
highlight = apply(temp,2, function(x) plot.table.helper.color(t(x)) )
layout(1)
plot.table(temp, smain = 'CDI Rank\nAlgo vs %D',highlight = highlight, colorbar = TRUE)
}
return(cdi)
}
custom.summary.positions <- function(weights) {
layout(1:len(weights))
for(w in names(weights)) {
tickers = colnames(weights[[w]])
n = len(tickers)
temp = matrix(NA, nr=4, nc=n)
colnames(temp) = tickers
rownames(temp) = spl('Avg Pos,Max Pos,Min Pos,# Periods')
temp['Avg Pos',] = 100 * apply(weights[[w]],2,mean,na.rm=T)
temp['Max Pos',] = 100 * apply(weights[[w]],2,max,na.rm=T)
temp['Min Pos',] = 100 * apply(weights[[w]],2,min,na.rm=T)
temp['# Periods',] = apply(weights[[w]] > 1/1000,2,sum,na.rm=T)
temp[] = plota.format(temp, 0, '', '')
plot.table(temp, smain=w)
}
}
custom.profit.chart <- function(data, main, cols) {
par(mar=c(4, 3, 2, 2),cex.main=2,cex.sub=1, cex.axis=1.5,cex.lab=1.5)
barplot(data, names.arg = names(data),
col=iif(data > 0, cols[1], cols[2]),
main=main,
cex.names = 1.5, border = 'darkgray',las=2)
grid(NA,NULL)
abline(h=0,col='black')
abline(h=mean(data),col='gray',lty='dashed',lwd=3)
}
custom.period.chart <- function(models) {
for(imodel in 1:len(models)) {
equity = models[[imodel]]$equity
period.ends = endpoints(equity, 'months')
period.ends = unique(c(1, period.ends[period.ends > 0]))
ret = equity[period.ends,] / mlag(equity[period.ends,]) - 1
ret = ret[-1]
ret.by.month = create.monthly.table(ret)
ret.by.month = 100 * apply(ret.by.month, 2, mean, na.rm=T)
period.ends = endpoints(equity, 'years')
period.ends = unique(c(1, period.ends[period.ends > 0]))
ret = equity[period.ends,] / mlag(equity[period.ends,]) - 1
ret.by.year = ret[-1]
ilayout =
'1,1
2,2
2,2
2,2
2,2
2,2
3,4
3,4
3,4
3,4
3,4
3,4'
plota.layout(ilayout)
make.table(1,1)
a = matrix(names(models)[imodel],1,1)
cex = plot.table.helper.auto.adjust.cex(a)
draw.cell(a[1],1,1, text.cex=cex,frame.cell=F)
temp = plotbt.monthly.table(equity)
cols = spl('green,red')
custom.profit.chart(ret.by.month, 'Average Monthly Returns', cols)
ret = 100*as.vector(ret.by.year)
names(ret) = date.year(index(ret.by.year))
custom.profit.chart(ret, 'Annual Returns', cols)
}
}
custom.returns.kpi <- function
(
bt,
trade.summary = NULL
)
{
out = list()
w = bt$period.weight
rc = bt$risk.contribution
out[[ 'Avg #' ]] =  mean(rowSums(w > 1/1000)) / 100
out[[ 'W Gini' ]] = mean(portfolio.concentration.gini.coefficient(w), na.rm=T)
out[[ 'W Herfindahl' ]] = mean(portfolio.concentration.herfindahl.index(w), na.rm=T)
out[[ 'RC Gini' ]] = mean(portfolio.concentration.gini.coefficient(rc), na.rm=T)
out[[ 'RC Herfindahl' ]] = mean(portfolio.concentration.herfindahl.index(rc), na.rm=T)
out = lapply(out, function(x) if(is.double(x)) round(100*x,1) else x)
out = c(bt.detail.summary(bt)$System, out)
return( list(System=out))
}
min.corr.paper.numerical.examples <- function()
{
n = 3
ia = list()
ia$n = 3
ia$risk = c(14, 18, 22) / 100;
ia$correlation = matrix(
c(1, 0.90, 0.85,
0.90, 1, 0.70,
0.85, 0.70, 1), nr=3, byrow=T)
ia$cov = ia$correlation * (ia$risk %*% t(ia$risk))
constraints = new.constraints(n)
constraints = new.constraints(n, lb = 0, ub = 1)
constraints = add.constraints(diag(n), type='>=', b=0, constraints)
constraints = add.constraints(diag(n), type='<=', b=1, constraints)
constraints = add.constraints(rep(1, n), 1, type = '=', constraints)
x = min.var.portfolio(ia, constraints)
sol = solve.QP(Dmat=ia$cov, dvec=rep(0, ia$n),
Amat=constraints$A, bvec=constraints$b, meq=constraints$meq)
x = sol$solution
round(x,4)
sqrt(x %*% ia$cov %*% x)
x %*% ia$cov
sol = solve.QP(Dmat=ia$correlation, dvec=rep(0, ia$n),
Amat=constraints$A, bvec=constraints$b, meq=constraints$meq)
x = sol$solution
round(x,4)
x %*% ia$correlation
x = x / sqrt( diag(ia$cov) )
x = x / sum(x)
round(x,4)
sqrt(x %*% ia$cov %*% x)
upper.index = upper.tri(ia$correlation)
cor.m = ia$correlation[upper.index]
cor.mu = mean(cor.m)
cor.sd = sd(cor.m)
norm.dist.m = 0 * ia$correlation
diag(norm.dist.m) = NA
norm.dist.m[upper.index] = sapply(cor.m, function(x) 1-pnorm(x, cor.mu, cor.sd))
norm.dist.m = (norm.dist.m + t(norm.dist.m))
norm.dist.avg = apply(norm.dist.m, 1, mean, na.rm=T)
norm.dist.rank = rank(-norm.dist.avg)
adjust.factor = 1
adjusted.norm.dist.rank = norm.dist.rank ^ adjust.factor
norm.dist.weight = adjusted.norm.dist.rank / sum(adjusted.norm.dist.rank)
weighted.norm.dist.average = norm.dist.weight %*% ifna(norm.dist.m,0)
final.weight = weighted.norm.dist.average / sum(weighted.norm.dist.average)
x = final.weight
x = x / sqrt( diag(ia$cov) )
x = x / sum(x)
round(x,4)
x = as.vector(x)
sqrt(x %*% ia$cov %*% x)
cor.m = ia$correlation
diag(cor.m) = 0
avg = rowMeans(cor.m)
cor.mu = mean(avg)
cor.sd = sd(avg)
norm.dist.avg = 1-pnorm(avg, cor.mu, cor.sd)
norm.dist.rank = rank(-norm.dist.avg)
adjust.factor = 1
adjusted.norm.dist.rank = norm.dist.rank ^ adjust.factor
norm.dist.weight = adjusted.norm.dist.rank / sum(adjusted.norm.dist.rank)
weighted.norm.dist.average = norm.dist.weight %*% (1-cor.m)
final.weight = weighted.norm.dist.average / sum(weighted.norm.dist.average)
x = final.weight
x = x / sqrt( diag(ia$cov) )
x = x / sum(x)
round(x,4)
x = as.vector(x)
sqrt(x %*% ia$cov %*% x)
}
solve.LP.bounds <- function
(
direction,
objective.in,
const.mat,
const.dir,
const.rhs,
binary.vec = 0,
lb = 0,
ub = +Inf,
default.lb = -100
)
{
n = len(objective.in)
if( len(lb) == 1 ) lb = rep(lb, n)
if( len(ub) == 1 ) ub = rep(ub, n)
lb = ifna(lb, default.lb)
ub = ifna(ub, +Inf)
lb[ lb < default.lb ] = default.lb
dvec = lb
index = which( ub < +Inf )
if( len(index) > 0 ) {
const.rhs = c(const.rhs, ub[index])
const.dir = c(const.dir, rep('<=', len(index)))
const.mat = rbind(const.mat, diag(n)[index, ])
}
if ( binary.vec[1] == 0 ) {
sol = lp( direction, objective.in, const.mat, const.dir,
const.rhs - const.mat %*% dvec )
} else {
dvec[binary.vec] = 0
sol = lp( direction, objective.in, const.mat, const.dir,
const.rhs - const.mat %*% dvec, binary.vec = binary.vec )
}
sol$solution = sol$solution + dvec
sol$value = objective.in %*% sol$solution
return( sol )
}
solve.QP.bounds <- function
(
Dmat,
dvec,
Amat,
bvec,
meq=0,
factorized=FALSE,
binary.vec = 0,
lb = -Inf,
ub = +Inf
)
{
Amat1 = Amat
bvec1 = bvec
n = len(dvec)
if( len(lb) == 1 ) lb = rep(lb, n)
if( len(ub) == 1 ) ub = rep(ub, n)
lb = ifna(lb, -Inf)
ub = ifna(ub, +Inf)
index = which( ub < +Inf )
if( len(index) > 0 ) {
bvec = c(bvec, -ub[index])
Amat = cbind(Amat, -diag(n)[, index])
}
index = which( lb > -Inf )
if( len(index) > 0 ) {
bvec = c(bvec, lb[index])
Amat = cbind(Amat, diag(n)[, index])
}
if ( binary.vec[1] == 0 ) {
qp.data.final = solve.QP.remove.equality.constraints(Dmat, dvec, Amat, bvec, meq)
Dmat = qp.data.final$Dmat
dvec = qp.data.final$dvec
Amat = qp.data.final$Amat
bvec = qp.data.final$bvec
meq = qp.data.final$meq
sol = try(solve.QP(Dmat, dvec, Amat, bvec, meq, factorized),TRUE)
if(inherits(sol, 'try-error')) {
ok = F
sol = list()
} else {
tol = 1e-3
ok = T
check = sol$solution %*% Amat - bvec
if(meq > 0) ok = ok & all(abs(check[1:meq]) <= tol)
ok = ok & all(check[-c(1:meq)] > -tol)
}
if(!ok) {
require(kernlab)
index.constant.variables = which(!is.na(qp.data.final$solution))
if( len(index.constant.variables) > 0 ) {
Amat1 = Amat[,1:ncol(Amat1)]
bvec1 = bvec[1:ncol(Amat1)]
lb = lb[-index.constant.variables]
ub = ub[-index.constant.variables]
}
sv = ipop(c = matrix(-dvec), H = Dmat, A = t(Amat1),
b = bvec1, l = ifna(lb,-100), u = ifna(ub,100),
r = c(rep(0,meq), rep(100, len(bvec1) - meq))
)
sol$solution = primal(sv)
}
x = qp.data.final$solution
x[qp.data.final$var.index] = sol$solution
sol$solution = x
} else {
qp_data = qp_new(binary.vec, Dmat = Dmat, dvec = dvec,
Amat=Amat, bvec=bvec, meq=meq, factorized=factorized)
sol = binary_branch_bound(binary.vec, qp_data, qp_solve,
control = bbb_control(silent=T, branchvar='max', searchdir='best' ))
qp_delete(qp_data)
sol$value = sol$fmin
sol$solution = sol$xmin
}
return(sol)
}
bbb_control <- function
(
itermax = 200,
depthmax = Inf,
bineps = 1e-4,
precisioneps = 0,
silent = T,
branchvar = c('first', 'max','min'),
proborder = c('0', '1', 'mindiff'),
searchdir = c('depth', 'breadth', 'best', 'normbest')
)
{
branchvar = switch(branchvar[1],
'first' = 0,
'max' = 1,
'min' = 2,
0)
branchvar = iif( is.null(branchvar),0, branchvar)
proborder = switch(proborder[1],
'0' = 0,
'1' = 1,
'mindiff' = 2,
0)
proborder = iif( is.null(proborder),0, proborder)
searchdir = switch(searchdir[1],
'depth' = 0,
'breadth' = 1,
'best' = 2,
'normbest' = 2,
0)
searchdir = iif( is.null(searchdir),0, searchdir)
control = list(itermax = itermax, depthmax = depthmax, bineps = bineps, precisioneps = precisioneps, silent = silent,
branchvar = branchvar,
proborder = proborder,
searchdir = searchdir)
return(control)
}
qp_new <- function
(
index_binvar,
Dmat,
dvec,
Amat,
bvec,
meq = 0,
factorized = FALSE
)
{
nbinvar = length(index_binvar)
nx = nrow(Dmat)
nbvec = length(bvec)
Amat = cbind( Amat, diag(nx)[,index_binvar], -diag(nx)[,index_binvar] )
bvec = c(bvec, rep(0,nx)[index_binvar], rep(1,nx)[index_binvar] )
lb_bin_index = (1:nbinvar) + nbvec
ub_bin_index = (1:nbinvar) + nbvec + nbinvar
qp_data = new.env()
qp_data$Dmat = Dmat
qp_data$dvec = dvec
qp_data$Amat = Amat
qp_data$bvec = bvec
qp_data$meq = meq
qp_data$factorized = factorized
qp_data$x0 = rep(0,nx)
qp_data$lb_bin_index = lb_bin_index
qp_data$ub_bin_index = ub_bin_index
qp_data$lb = bvec[lb_bin_index]
qp_data$ub = bvec[ub_bin_index]
return(qp_data)
}
qp_delete <- function(qp_data)
{
rm(list = ls(qp_data,all=TRUE), envir = qp_data)
}
qp_solve <- function
(
qp_data,
lb,
ub
)
{
bvec = qp_data$bvec
bvec[qp_data$lb_bin_index] = lb
bvec[qp_data$ub_bin_index] = -ub
qp.data.final = solve.QP.remove.equality.constraints(qp_data$Dmat, qp_data$dvec,
qp_data$Amat, bvec, qp_data$meq)
sol = tryCatch( solve.QP(Dmat=qp.data.final$Dmat, dvec=qp.data.final$dvec, Amat=qp.data.final$Amat,
bvec=qp.data.final$bvec, meq=qp.data.final$meq, factorized=qp_data$factorized),
error=function( err ) FALSE,
warning=function( warn ) FALSE )
if( !is.logical( sol ) ) {
x = qp.data.final$solution
x[qp.data.final$var.index] = sol$solution
return(list( ok = TRUE, x = x, fval = sol$value ))
} else {
return(list( ok = FALSE ))
}
}
mbqp.test <- function()
{
load.packages('quadprog')
Q = diag(4)
b	= c( 2, -3, -2, -3)
C	= matrix(c(-1,  -1,  -1,  -1,
10,	5,   3,	4,
-1,	0,   0,	0),
3,4,byrow = TRUE)
d = c(-2, 10,  0)
vlb  = c(-1e10, 0, 0, 0);
vub  = c( 1e10, 1, 1, 1);
index_binvar = c(2, 3, 4);
Dmat = Q
dvec = -b
Amat = -t(C)
bvec = -d
n = nrow(Dmat)
Amat = cbind( Amat, diag(n), -diag(n) )
bvec = c( bvec, vlb, -vub )
sol = solve.QP(Dmat=Dmat, dvec=dvec, Amat=Amat, bvec=bvec, meq=0)
xsol = sol$solution
fsol = sol$value
cat('QP.solve fsol =', fsol, 'xsol =', xsol, '\n')
sol = solve.QP.bounds(Dmat=Dmat, dvec=dvec, Amat=Amat, bvec=bvec, meq=0, binary.vec = index_binvar)
xsol = sol$solution
fsol = sol$value
cat('QP.solve Binary Branch and Bound fsol =', fsol, 'xsol =', xsol, '\n')
}
solve.QP.remove.equality.constraints <- function
(
Dmat,
dvec,
Amat,
bvec,
meq=0
)
{
qp.data = list()
qp.data$Amat = Amat
qp.data$bvec = bvec
qp.data$Dmat = Dmat
qp.data$dvec = dvec
qp.data$meq = meq
Amat1 = t(qp.data$Amat)
bvec1 = qp.data$bvec
Dmat1 = qp.data$Dmat
dvec1 = qp.data$dvec
meq1 = qp.data$meq
qp.data$solution = rep(NA, ncol(Amat1))
qp.data$var.index = 1:ncol(Amat1)
while(T) {
one.non.zero.index = which( rowSums(Amat1!=0) == 1 )
if( len(one.non.zero.index) == 0 ) break
temp0 = rowSums(Amat1[one.non.zero.index,])
temp = abs( temp0 )
bvec1[one.non.zero.index] = bvec1[one.non.zero.index] / temp
Amat1[one.non.zero.index,] = Amat1[one.non.zero.index,] / temp
temp0.index = matrix(1:ncol(Amat1), nr=ncol(Amat1), nc=len(one.non.zero.index))[t(Amat1[one.non.zero.index,]!=0)]
equality.constraints = rep(NA, ncol(Amat1))
lb = ub = rep(NA, ncol(Amat1))
index = temp0 > 0
temp = order(bvec1[one.non.zero.index[index]], decreasing = FALSE)
lb[temp0.index[index][temp]] = bvec1[one.non.zero.index[index]][temp]
index = temp0 < 0
temp = order(-bvec1[one.non.zero.index[index]], decreasing = TRUE)
ub[temp0.index[index][temp]] = -bvec1[one.non.zero.index[index]][temp]
remove.index = which(lb == ub)
if( len(remove.index) > 0 ) {
equality.constraints[remove.index] = lb[remove.index]
Dmat1 = Dmat1[-remove.index, -remove.index,drop=F]
dvec1 = dvec1[-remove.index]
bvec1 = bvec1 - Amat1[,remove.index,drop=F] %*% equality.constraints[remove.index]
Amat1 = Amat1[,-remove.index,drop=F]
qp.data$solution[ qp.data$var.index[remove.index] ] = lb[remove.index]
qp.data$var.index = which(is.na(qp.data$solution))
if( ncol(Amat1) > 0 ) {
remove.index = which( rowSums(Amat1!=0) == 0 & bvec1 == 0 )
if(len(remove.index)>0) {
bvec1 = bvec1[-remove.index]
Amat1 = Amat1[-remove.index,,drop=F]
if( meq1 > 0 ) meq1 = meq1 - len(intersect((1:meq1), remove.index))
}
} else break
} else break
}
qp.data$Amat = t(Amat1)
qp.data$bvec = bvec1
qp.data$Dmat = Dmat1
qp.data$dvec = dvec1
qp.data$meq = meq1
return(qp.data)
}
solve.QP.remove.equality.constraints12 <- function
(
Dmat,
dvec,
Amat,
bvec,
meq=0
)
{
qp.data.temp = list()
qp.data.temp$Amat = Amat
qp.data.temp$bvec = bvec
qp.data.temp$Dmat = Dmat
qp.data.temp$dvec = dvec
qp.data.temp$meq = meq
qp.data.temp = remove.equality.constraints.old(qp.data.temp)
if( len(qp.data.temp$var.index) == len(qp.data.temp$solution) ) {
qp.data.final = qp.data.temp
} else {
qp.data.final = remove.equality.constraints.old(qp.data.temp)
qp.data.temp$solution[qp.data.temp$var.index] = qp.data.final$solution
qp.data.final$solution = qp.data.temp$solution
qp.data.final$var.index = qp.data.temp$var.index[qp.data.final$var.index]
}
return(qp.data.final)
}
remove.equality.constraints.old <- function(qp.data)
{
Amat1 = qp.data$Amat
bvec1 = qp.data$bvec
Dmat1 = qp.data$Dmat
dvec1 = qp.data$dvec
meq1 = qp.data$meq
one.non.zero.index = which( colSums(Amat1!=0) == 1 )
if( len(one.non.zero.index) == 0 ) {
equality.constraints = rep(NA, nrow(Amat1))
qp.data$solution = equality.constraints
qp.data$var.index = which(is.na(equality.constraints))
return(qp.data)
}
bvec1[one.non.zero.index] = bvec1[one.non.zero.index] / abs( colSums(Amat1[,one.non.zero.index]) )
Amat1[,one.non.zero.index] = Amat1[,one.non.zero.index] / abs( Amat1[,one.non.zero.index] )
Amat1[is.na(Amat1)] = 0
equality.constraints = rep(NA, nrow(Amat1))
remove.constraints = rep(NA, ncol(Amat1))
for( r in which(rowSums(Amat1[,one.non.zero.index]!=0) > 1) ) {
temp.index = which( Amat1[,one.non.zero.index][r,] != 0 )
for( r1 in temp.index[-1] ) {
temp.index1 = colSums(abs(
rbind(Amat1,bvec1)[,one.non.zero.index[temp.index]] +
rbind(Amat1,bvec1)[,one.non.zero.index[r1]])
)
if( any(temp.index1 == 0 ) ) {
equality.constraints[r] =
bvec1[one.non.zero.index[r1]] / Amat1[r,one.non.zero.index[r1]]
remove.constraints[ one.non.zero.index[r1] ] = 1
remove.constraints[ one.non.zero.index[ temp.index[ temp.index1 == 0]]] = 1
break;
}
}
}
remove.index = which(!is.na(equality.constraints))
if(len(remove.index)>0) {
Dmat1 = Dmat1[-remove.index, -remove.index,drop=F]
dvec1 = dvec1[-remove.index]
bvec1 = bvec1 - equality.constraints[remove.index] %*% Amat1[remove.index,,drop=F]
Amat1 = Amat1[-remove.index,,drop=F]
remove.index1 = which(remove.constraints==1)
if(len(remove.index1)>0) {
bvec1 = bvec1[-remove.index1]
Amat1 = Amat1[,-remove.index1,drop=F]
if( meq1 > 0 ) meq1 = meq1 - len(intersect((1:meq1), remove.index1))
}
}
qp.data$Amat = Amat1
qp.data$bvec = bvec1
qp.data$Dmat = Dmat1
qp.data$dvec = dvec1
qp.data$meq = meq1
qp.data$solution = equality.constraints
qp.data$var.index = which(is.na(equality.constraints))
return(qp.data)
}
lm.constraint <- function
(
x,
y,
constraints = NULL
)
{
if( is.null(constraints) ) {
fit = lm.fit(x, y)
return( ols.summary(x, y, fit$coefficients) )
} else {
temp = cov(cbind(y, x))
Dmat = temp[-1,-1]
dvec = temp[-1,1]
sol = solve.QP.bounds(Dmat = Dmat, dvec = dvec ,
Amat=constraints$A, bvec=constraints$b, constraints$meq,
lb = constraints$lb, ub = constraints$ub)
return( ols.summary(x, y, sol$solution) )
}
}
ols <- function
(
x,
y,
computeSummary=F
)
{
xx = t(x) %*% x
if(is.null(ncol(xx))) { xinv = inv1(xx)
} else if(ncol(xx) == 1) { xinv = inv1(xx)
} else if(ncol(xx) == 2) { xinv = inv2(xx)
} else if(ncol(xx) == 3) { xinv = inv3(xx)
} else { xinv = inv(xx) }
coefficients = xinv %*% t(x) %*% y
if(computeSummary) {
return( ols.summary(x, y, coefficients, xinv) )
} else {
return(list(coefficients = coefficients))
}
}
ols.summary <- function
(
x,
y,
coefficients,
xinv = NULL
)
{
n = length(y)
p = length(coefficients)
rdf = n-p
e = y - x %*% coefficients
ess=sum(e^2)
mss = sum((y - sum(y)/n)^2)
r.squared = 1 - ess/mss
if( !is.null(xinv) ) {
s2=ess/(rdf)
seb=sqrt(diag(s2*xinv))
tratio=coefficients/seb
return(list(coefficients = coefficients, seb = seb,tratio = tratio, r.squared = r.squared))
} else {
return(list(coefficients = coefficients, r.squared = r.squared))
}
}
ols.test <- function() {
x = matrix( rnorm(4*10), ncol=4)
y = rnorm(10)
summary(lm(y ~ x+0))
ols(x, y, T)
}
inv <- function(x) { solve(x) }
inv1 <- function(x) { 1/x }
inv2 <- function(x)
{
matrix(c(x[2,2],-x[1,2],-x[2,1],x[1,1]),nrow=2,byrow=T) / (x[1,1]*x[2,2] - x[1,2]*x[2,1])
}
inv3 <- function(x)
{
matrix(c(x[3,3]*x[2,2]-x[3,2]*x[2,3],-(x[3,3]*x[1,2]-x[3,2]*x[1,3]),x[2,3]*x[1,2]-x[2,2]*x[1,3],
-(x[3,3]*x[2,1]-x[3,1]*x[2,3]),x[3,3]*x[1,1]-x[3,1]*x[1,3],-(x[2,3]*x[1,1]-x[2,1]*x[1,3]),
x[3,2]*x[2,1]-x[3,1]*x[2,2],-(x[3,2]*x[1,1]-x[3,1]*x[1,2]),x[2,2]*x[1,1]-x[2,1]*x[1,2]),nrow=3,byrow=T) /
(x[1,1]*(x[3,3]*x[2,2]-x[3,2]*x[2,3])-x[2,1]*(x[3,3]*x[1,2]-x[3,2]*x[1,3])+x[3,1]*(x[2,3]*x[1,2]-x[2,2]*x[1,3]))
}
inv.test <- function() {
m=matrix(c(4,3,3,2),nrow=2,byrow=T)
inv2(m) %*% m
inv(m) %*% m
m = matrix(c(1,2,3,4,5,6,7,8,8),ncol=3,byrow=T)
inv3(m) %*% m
m %*% inv3(m)
inv(m) %*% m
}
find.maximum.distance.point <- function
(
y,
x=1:len(y)
)
{
allCoord = rbind(vec(y), vec(x))
firstPoint = allCoord[,1]
lineVec = allCoord[,len(y)] - firstPoint
lineVecN = lineVec / sqrt(sum(lineVec^2))
vecFromFirst = allCoord - firstPoint
scalarProduct = lineVecN %*% vecFromFirst
vecFromFirstParallel = t(scalarProduct) %*% lineVecN
vecToLine = t(vecFromFirst) - vecFromFirstParallel
distToLine = sqrt(rowSums(vecToLine^2,2))
which.max(distToLine)
}
make.table <- function
(
nr,
nc
)
{
savepar = par(mar = rep(1, 4))
plot(c(0.5, nc*2 + 0.5), c(-0.5, -(nr + 0.5)), xaxs = 'i', yaxs = 'i',
type = 'n', xlab = '', ylab = '', axes = FALSE)
savepar
}
draw.cell <- function
(
title,
r,
c,
text.cex = 1,
bg.col = 'white',
frame.cell = T
)
{
if(!frame.cell) bcol = bg.col else bcol = 'black'
rect((2*(c - 1) + .5), -(r - .5), (2*c + .5), -(r + .5), col = bg.col, border = bcol)
if( c == 1) {
text((2*(c - 1) + .5), -r, title, adj = 0, cex = text.cex)
} else if( r == 1 ) {
text((2*(c - 1) + .5), -r, title, adj = 0, cex = text.cex)
} else {
text((2*c + .5), -r, title, adj = 1, cex = text.cex)
}
}
plot.table.helper.auto.adjust.cex <- function
(
temp.table,
keep.all.same.cex = FALSE
)
{
nr = nrow(temp.table)
nc = ncol(temp.table)
all.xrange = diff(par()$usr[1:2]) / nc
xrange = matrix( strwidth(paste('  ', temp.table), units = 'user', cex = 1), nc = nc)
all.yrange = diff(par()$usr[3:4]) / nr
yrange = matrix( 5/3 * strheight(temp.table, units = 'user', cex = 1), nc = nc)
plot.matrix.cex = pmin( round(all.yrange / yrange, 2) , round(all.xrange / xrange, 2) )
header.col.cex = min(plot.matrix.cex[1,-1])
header.row.cex = min(plot.matrix.cex[-1,1])
title.cex = plot.matrix.cex[1, 1]
data.cex = min(plot.matrix.cex[-1, -1])
if ( keep.all.same.cex ) {
plot.matrix.cex[] = min(plot.matrix.cex)
} else {
plot.matrix.cex[1,-1] = min(c(header.col.cex, header.row.cex))
plot.matrix.cex[-1,1] = min(c(header.col.cex, header.row.cex))
plot.matrix.cex[-1,-1]= min(c(header.col.cex, header.row.cex, data.cex))
plot.matrix.cex[1,1]= min(c(header.col.cex, header.row.cex, data.cex, title.cex))
plot.matrix.cex[1,-1] = min(c(header.col.cex))
plot.matrix.cex[-1,1] = min(c(header.row.cex))
plot.matrix.cex[-1,-1]= min(c(data.cex))
plot.matrix.cex[1,1]= min(c(title.cex))
}
return(plot.matrix.cex)
}
plot.table.param <- function
(
plot.matrix,
smain = '',
plot.matrix.cex,
plot.matrix_bg.col,
frame.cell = T,
keep.all.same.cex = FALSE
)
{
n = nrow(plot.matrix)
pages = unique(c(seq(0, n, by = 120), n))
for(p in 1:(len(pages)-1)) {
rindex = (pages[p]+1) : pages[p+1]
temp.table = matrix('', nr = len(rindex)+1, nc = ncol(plot.matrix)+1)
temp.table[-1, -1] = plot.matrix[rindex,]
temp.table[1, -1] = colnames(plot.matrix)
temp.table[-1, 1] = rownames(plot.matrix)[rindex]
temp.table[1, 1] = smain
nr = nrow(temp.table)
nc = ncol(temp.table)
par(mar = c(0, 0, 0, 0), cex = 0.5)
oldpar = make.table(nr, nc)
text.cex = plot.matrix.cex[c(1, 1 + rindex), ]
text.cex = plot.table.helper.auto.adjust.cex(temp.table, keep.all.same.cex)
bg.col = plot.matrix_bg.col[c(1, 1 + rindex), ]
for(r in 1:nr) {
for(c in 1:nc) {
draw.cell( paste('', temp.table[r,c], '', sep=' '), r, c,
text.cex = text.cex[r,c], bg.col = bg.col[r,c], frame.cell = frame.cell)
}
}
}
}
plot.table.helper.color <- function
(
temp
){
temp = matrix(as.double(gsub('[%,$]', '', temp)), nrow(temp), ncol(temp))
highlight = as.vector(temp)
cols = rep(NA, len(highlight))
ncols = len(highlight[!is.na(highlight)])
cols[1:ncols] = rainbow(ncols, start = 0, end = 0.3)
o = sort.list(highlight, na.last = TRUE, decreasing = FALSE)
o1 = sort.list(o, na.last = TRUE, decreasing = FALSE)
highlight = matrix(cols[o1], nrow = nrow(temp))
highlight[is.na(temp)] = NA
return(highlight)
}
plot.table.helper.colorbar <- function
(
plot.matrix
)
{
nr = nrow(plot.matrix) + 1
nc = ncol(plot.matrix) + 1
c = nc
r1 = 1
r2 = nr
rect((2*(c - 1) + .5), -(r1 - .5), (2*c + .5), -(r2 + .5), col='white', border='white')
rect((2*(c - 1) + .5), -(r1 - .5), (2*(c - 1) + .5), -(r2 + .5), col='black', border='black')
y1= c( -(r2) : -(r1) )
graphics::image(x = c(  (2*(c - 1) + 1.5) : (2*c + 0.5) ),
y   = y1,
z   = t(matrix(  y1  , ncol = 1)),
col = t(matrix( rainbow(len( y1  ), start = 0, end = 0.3) , ncol = 1)),
add = T)
}
plot.table <- function
(
plot.matrix,
smain = NULL,
text.cex = 1,
frame.cell = T,
highlight = F,
colorbar = FALSE,
keep_all.same.cex = FALSE
)
{
if( is.null(rownames(plot.matrix)) & is.null(colnames(plot.matrix)) ) {
temp.matrix = plot.matrix
if( nrow(temp.matrix) == 1 ) temp.matrix = rbind('', temp.matrix)
if( ncol(temp.matrix) == 1 ) temp.matrix = cbind('', temp.matrix)
plot.matrix = temp.matrix[-1, -1, drop = FALSE]
colnames(plot.matrix) = temp.matrix[1, -1]
rownames(plot.matrix) = temp.matrix[-1, 1]
smain = iif(is.null(smain), temp.matrix[1, 1], smain)
} else if( is.null(rownames(plot.matrix)) ) {
temp.matrix = plot.matrix
if( ncol(plot.matrix) == 1 ) temp.matrix = cbind('', temp.matrix)
plot.matrix = temp.matrix[, -1, drop = FALSE]
colnames(plot.matrix) = colnames(temp.matrix)[-1]
rownames(plot.matrix) = temp.matrix[,1]
smain = iif(is.null(smain), colnames(temp.matrix)[1], smain)
} else if( is.null(colnames(plot.matrix)) ) {
temp.matrix = plot.matrix
if( nrow(temp.matrix) == 1 ) temp.matrix = rbind('', temp.matrix)
plot.matrix = temp.matrix[-1, , drop = FALSE]
rownames(plot.matrix) = rownames(temp.matrix)[-1]
colnames(plot.matrix) = temp.matrix[1, ]
smain = iif(is.null(smain), rownames(temp.matrix)[1], smain)
}
smain = iif(is.null(smain), '', smain)
plot.matrix[which(trim(plot.matrix) == 'NA')] = ''
plot.matrix[which(trim(plot.matrix) == 'NA%')] = ''
plot.matrix[which(is.na(plot.matrix))] = ''
if(colorbar) {
plot.matrix = cbind(plot.matrix, '')
if(!is.null(highlight)) if(!is.logical(highlight)) { highlight = cbind(highlight, NA) }
}
nr = nrow(plot.matrix) + 1
nc = ncol(plot.matrix) + 1
is_highlight = T
if(is.logical(highlight)) {
is_highlight = highlight
if(highlight) highlight = plot.table.helper.color(plot.matrix)
}
if(!is_highlight) {
plot.matrix.cex = matrix(1, nr = nr, nc = nc )
plot.matrix_bg.col = matrix('white', nr = nr, nc = nc )
plot.matrix_bg.col[seq(1, nr, 2), ] = 'yellow'
plot.matrix_bg.col[1,] = 'gray';
plot.table.param( plot.matrix, smain, plot.matrix.cex, plot.matrix_bg.col,
frame.cell, keep_all.same.cex)
} else {
plot.matrix.cex = matrix(1, nr = nr, nc = nc )
plot.matrix_bg.col = matrix('white', nr = nr, nc = nc )
plot.matrix_bg.col[1,] = 'gray'
plot.matrix_bg.col[2:nr,2:nc] = highlight
plot.table.param(plot.matrix, smain, plot.matrix.cex, plot.matrix_bg.col,
frame.cell, keep_all.same.cex)
}
if(colorbar) plot.table.helper.colorbar(plot.matrix);
}
plot.table.test <- function()
{
mrownames = spl('row one,row two,row 3')
mcolnames = spl('col 1,col 2,col 3,col 4')
temp = matrix(NA, len(mrownames), len(mcolnames))
rownames(temp) = mrownames
colnames(temp) = mcolnames
temp[,] = matrix(1:12,3,4)
png(filename = 'plot1.png', width = 500, height = 500, units = 'px', pointsize = 12, bg = 'white')
plot.table(temp, format(as.Date(Sys.time()), '%d %b %Y'))
dev.off()
data =  matrix(rnorm(1000), nc=10)
colnames(data) = paste('data', 1:10, sep='')
temp = cor(data, use='complete.obs', method='pearson')
temp[] = plota.format(100 * temp, 0, '', '%')
png(filename = 'plot2.png', width = 500, height = 500, units = 'px', pointsize = 12, bg = 'white')
plot.table(temp, smain='Correlation', highlight = TRUE, colorbar = TRUE)
dev.off()
}
plot.periodic.table1 <- function(hist.returns)
{
n = ncol(hist.returns)
temp = t(coredata(hist.returns))
colnames(temp) = format(index.xts(hist.returns), '%Y')
rownames(temp) = 1:n
rownames(temp)[1] = ' Best '
rownames(temp)[n] = ' Worst '
col = plota.colors(n)
highlight = apply(temp,2, function(x) col[order(x, decreasing = T)] )
temp[] = apply(temp,2, sort, decreasing = T)
temp[] = plota.format(100 * temp, 0, '', '%')
plot.table(temp, highlight = highlight)
plota.legend(colnames(hist.returns), col)
}
plot.periodic.table2 <- function(hist.returns)
{
temp = t(coredata(hist.returns))
colnames(temp) = format(index.xts(hist.returns), '%Y')
temp[] = plota.format(100 * temp, 0, '', '%')
highlight = apply(temp,2, function(x) plot.table.helper.color(t(x)) )
plot.table(temp, highlight = highlight, colorbar = TRUE)
}
plota.theme <- function
(
col.border = 'black',
col.up = 'green',
col.dn = 'red',
col.x.highlight = 'orange',
col.y.highlight = 'orange',
alpha=NA
)
{
col = c(col.border, col.up, col.dn, col.x.highlight, col.y.highlight)
if(!is.na(alpha)) col = col.add.alpha(col, alpha)
plota.control$col.border = col[1]
plota.control$col.up = col[2]
plota.control$col.dn = col[3]
plota.control$col.x.highlight = col[4]
plota.control$col.y.highlight = col[5]
}
plota.theme.blue.red <- function(alpha=NA)
{
plota.theme(
col.border = 'black',
col.up = 'blue',
col.dn = 'red',
alpha = alpha
)
}
plota.theme.green.orange <- function(alpha=NA)
{
plota.theme(
col.border = rgb(68,68,68, maxColorValue=255),
col.up = rgb(0,204,0, maxColorValue=255),
col.dn = rgb(255,119,0, maxColorValue=255),
alpha = alpha
)
}
plota.theme.gray.orange <- function(alpha=NA)
{
plota.theme(
col.border = '#444444',
col.up = '#BEBEBE',
col.dn = '#FF7700',
alpha = alpha
)
}
plota.control = new.env()
plota.control$col.border = 'black'
plota.control$col.up = 'green'
plota.control$col.dn = 'red'
plota.control$col.x.highlight = 'orange'
plota.control$col.y.highlight = 'orange'
plota.control$xaxis.ticks = c()
plota.theme.green.orange();
col.add.alpha <- function
(
col,
alpha=150
)
{
rgb(t(col2rgb(col)), alpha=alpha, maxColorValue = 255)
}
plota <- function
(
y,
main = NULL,
plotX = TRUE,
LeftMargin = 0,
x.highlight = NULL,
y.highlight = NULL,
las = 1,
type = 'n',
xlab = '',
ylab = '',
ylim = NULL,
log = '',
...
)
{
hasTitle = !is.null(main);
par( mar = c(iif(plotX,2,0), LeftMargin , iif(hasTitle,2,0), 3) )
if(has.Cl(y)) y1 = Cl(y) else y1 = y[,1]
if( is.null(ylim) ) {
ylim = range(y1, na.rm = T)
switch(type,
'ohlc' = ,
'hl' = ,
'candle' = { ylim = range(OHLC(y), na.rm = T) },
'volume' = { y1 = Vo(y); ylim = range(Vo(y), na.rm = T) }
)
}
temp.x = attr(y, 'index')
plot( temp.x, y1, xlab = xlab, ylab = ylab, main = main,
type = 'n', yaxt = 'n', xaxt = 'n', ylim = ylim, log = log, ... )
axis(4, las = las)
class(temp.x) = c('POSIXct', 'POSIXt')
plota.control$xaxis.ticks = axis.POSIXct(1, temp.x,labels = plotX, tick = plotX)
if( !is.null(x.highlight) ) plota.x.highlight(y, x.highlight);
if( !is.null(y.highlight) ) plota.y.highlight(y, y.highlight);
plota.grid()
switch(type,
'candle' = plota.candle(y, ...),
'hl' = plota.hl(y, ...),
'ohlc' = plota.ohlc(y, ...),
'volume' = plota.volume(y, ...),
{  lines(temp.x, y1, type=type, ...) }
)
box();
}
plota2Y <- function(
y,
las = 1,
type = 'n',
...
)
{
xlim = par('usr')[1:2]
class(xlim) = c('POSIXct', 'POSIXt')
y1 = y[paste(format(xlim, '%Y:%m:%d %H:%M:%S'), sep = '', collapse = '::')]
par(new = TRUE)
xlim = par('usr')[1:2]
plot( attr(y1, 'index') , y1[,1], xlim = xlim, xaxs = 'i', type = type,
yaxt = 'n', xaxt = 'n', xlab = '', ylab = '', axes = F, ... )
axis(2, las = las, ...)
}
plota.grid <- function()
{
abline( h = axTicks(2), col = 'lightgray', lty = 'dotted')
abline( v = plota.control$xaxis.ticks, col = 'lightgray', lty = 'dotted')
}
plota.lines <- function(
y,
type = 'l',
col = par('col'),
...
)
{
if(has.Cl(y)) y1 = Cl(y) else y1 = y[,1]
temp.x = attr(y, 'index')
if( type == 'l' & len(col) > 1 ) {
for( icol in unique(col) ) {
lines(temp.x, iif(col == icol, y1, NA), type = type, col = icol, ...)
}
} else {
lines(temp.x, y1, type = type, col = col, ...)
}
}
plota.text <- function(
y,
...
)
{
if(has.Cl(y)) y1 = Cl(y) else y1 = y[,1]
text(index4xts(y1), y1, ...)
}
plota.format <- function(
temp,
nround = 2,
sprefix = '',
eprefix = ''
)
{
return( paste(sprefix,
format(round(as.numeric(temp), nround), big.mark = ',', scientific=FALSE),
eprefix ,sep='') )
}
plota.legend <- function
(
labels,
fill = NULL,
lastobs = NULL,
x = 'topleft',
merge = F,
bty = 'n',
yformat = plota.format,
...
)
{
if( !is.null(fill) ) fill = spl( as.character(fill) )
labels = spl( as.character(labels) )
if( !is.null(lastobs) ) {
if( is.list(lastobs) ) {
labels1 = sapply(lastobs, function(x) unclass(last(x))[1])
} else {
labels1 = unclass(last(lastobs))[1];
}
labels = paste(labels, match.fun(yformat)( labels1 ))
}
legend(x, legend = labels, fill = fill, merge = merge, bty = bty, ...)
}
plota.layout <- function(
ilayout,
delim = ','
)
{
ilayout = matrix( as.double(spl( gsub('\n', delim, ilayout), delim)),
nrow = len(spl(ilayout, '\n')), byrow=TRUE)
layout(mat = ilayout)
}
plota.dx <- function
(
y
)
{
xlim = par('usr')[1:2]
class(xlim) = c('POSIXct', 'POSIXt')
y1 = y[paste(format(xlim, '%Y:%m:%d %H:%M:%S'), sep = '', collapse = '::')]
xlim = par('usr')[1:2]
xportion = min(1, diff(unclass(range(attr(y1, 'index'))))*1.08 / diff(xlim) )
return( xportion * diff(xlim) / ( 2* nrow(y1)  ) )
}
plota.x.highlight <- function
(
y,
highlight,
col = plota.control$col.x.highlight
)
{
if(len(col)==1) {
plota.x.highlight.helper(y, highlight, col = col)
} else {
for( icol in unique(col[highlight]) ) {
plota.x.highlight.helper(y, iif(col == icol, highlight, FALSE), col = icol)
}
}
}
plota.x.highlight.helper <- function
(
y,
highlight,
col = plota.control$col.x.highlight
)
{
dx = plota.dx(y);
hl_index = highlight;
if( is.logical(highlight) ) hl_index = which(highlight);
if( identical(unique(highlight) , c(0, 1)) ) hl_index = which(as.logical(highlight));
hl_index1 = which(diff(hl_index) > 1 )
hl_index = hl_index[ sort(c(1, len(hl_index), hl_index1, (hl_index1+1))) ]
temp.y = par('usr')[3:4]
if(par('ylog')) temp.y = 10^temp.y
temp.x = attr(y, 'index')
for( i in seq(1,len(hl_index),2) ) {
rect(temp.x[hl_index[i]] - dx/2, temp.y[1],
temp.x[hl_index[(i + 1)]] + dx/2, temp.y[2],
col = col, border = col )
}
box();
}
plota.y.highlight <- function
(
y,
highlight,
col = plota.control$col.y.highlight
)
{
temp.y = par('usr')[3:4]
if(par('ylog')) temp.y = 10^temp.y
temp.x = par('usr')[1:2]
if(par('xlog')) temp.x = 10^temp.x
highlight[highlight == Inf] = temp.y[2]
highlight[highlight == -Inf] = temp.y[1]
for( i in seq(1,len(highlight),by=2) ) {
rect(temp.x[1], highlight[i],
temp.x[2], highlight[(i + 1)],
col = col, border = col )
}
box();
}
plota.candle.col <- function(	y ) {
return( iif( Cl(y)>Op(y), plota.control$col.up, plota.control$col.dn) )
}
plota.volume.col <- function( y ) {
return( iif( Cl(y)>mlag(Cl(y)), plota.control$col.up, plota.control$col.dn) )
}
plota.candle <- function
(
y,
col = plota.candle.col(y)
)
{
dx = plota.dx(y)
dxi0 = ( dx / xinch() ) * 96
if( dxi0 < 1 ) {
plota.hl.lwd(y, col = col, lwd = 1)
} else if ( dxi0 < 1.75 ) {
plota.ohlc.lwd(y, col = col, lwd = 1)
} else {
temp.x = attr(y, 'index')
rect(temp.x - dx/10, Lo(y), temp.x + dx/10, Hi(y),
col = plota.control$col.border, border = plota.control$col.border)
rect(temp.x - dx/2, Op(y), temp.x + dx/2, Cl(y),
col = col, border = plota.control$col.border)
}
}
plota.ohlc <- function
(
y,
col = plota.control$col.border
)
{
dx = plota.dx(y)
dxi0 = ( dx / xinch() ) * 96
if( dxi0 < 1 ) {
plota.hl.lwd(y, col = col, lwd = 1)
} else if ( dxi0 < 1.75 ) {
plota.ohlc.lwd(y, col = col, lwd = 1)
} else {
temp.x = attr(y, 'index')
rect(temp.x - dx/8, Lo(y), temp.x + dx/8, Hi(y), col = col, border = col)
segments(temp.x - dx/2, Op(y), temp.x, Op(y), col = col)
segments(temp.x + dx/2, Cl(y), temp.x, Cl(y), col = col)
}
}
plota.hl <- function
(
y,
col = plota.volume.col(y),
border = plota.control$col.border
)
{
dx = plota.dx(y)
dxi0 = ( dx / xinch() ) * 96
if( dxi0 < 1.75 ) {
plota.hl.lwd(y, col = col, lwd = 1)
} else {
temp.x = attr(y, 'index')
rect(temp.x - dx/2, Lo(y), temp.x + dx/2, Hi(y),
col = col, border = border)
}
}
plota.ohlc.lwd <- function
(
y,
lwd=1,
...
)
{
dx = plota.dx(y)
temp.x = attr(y, 'index')
segments(temp.x, Lo(y), temp.x, Hi(y), lwd = lwd, lend = 2,  ...)
segments(temp.x - dx/2, Op(y), temp.x, Op(y), lwd = lwd, lend = 2, ...)
segments(temp.x + dx/2, Cl(y), temp.x, Cl(y), lwd = lwd, lend = 2, ...)
}
plota.hl.lwd <- function
(
y,
lwd=1,
...
)
{
temp.x = attr(y, 'index')
segments(temp.x, Lo(y), temp.x, Hi(y), lwd = lwd, lend = 2, ...)
}
plota.volume <- function
(
y,
col = plota.volume.col(y),
border = plota.control$col.border
)
{
dx = plota.dx(y)
dxi0 = ( dx / xinch() ) * 96
temp.x = attr(y, 'index')
if( dxi0 < 1.75 ) {
segments(temp.x, 0, temp.x, Vo(y), col = col, lwd = 1, lend = 2)
} else {
rect(temp.x - dx/2, 0, temp.x + dx/2, Vo(y),
col = col, border = border)
}
idv = grep('Volume', colnames(y))
temp = spl(colnames(y)[idv], ';')
if( len(temp) > 1 ) legend('topright',legend = temp[len(temp)], bty='n');
}
plota.scale.volume <- function(y)
{
Volumes = Vo(y)
max.vol = max(Volumes, na.rm = T)
vol.scale = list(100, '100s')
if (max.vol > 10000)
vol.scale = list(1000, '1000s')
if (max.vol > 1e+05)
vol.scale = list(10000, '10,000s')
if (max.vol > 1e+06)
vol.scale = list(1e+05, '100,000s')
if (max.vol > 1e+07)
vol.scale = list(1e+06, 'millions')
idv = grep('Volume', colnames(y))
y[, idv] = Volumes/vol.scale[[1]]
colnames(y)[idv] = paste( colnames(y)[idv], vol.scale[[2]], sep=';' )
return(y)
}
plota.test <- function() {
load.packages('quantmod')
data.spy = getSymbols('SPY', from = '1980-01-01', auto.assign = FALSE)
data.ibm = getSymbols('IBM', from = '1980-01-01', auto.assign = FALSE)
y = data.spy['2011:01:01::2011:02:01']
highlight = which(Cl(y) < 127)
png(filename = 'plot1.png', width = 500, height = 500, units = 'px', pointsize = 12, bg = 'white')
layout(c(1,1,2))
plota(y, type = 'candle', main = 'SPY', plotX = F, x.highlight = highlight)
y = plota.scale.volume(y)
plota(y, type = 'volume', x.highlight = highlight)
dev.off()
y = data.spy['2010:01:01::2011:02:01']
png(filename = 'plot2.png', width = 500, height = 500, units = 'px', pointsize = 12, bg = 'white')
layout(c(1,1,2,3))
plota(y, type = 'candle', plotX = F)
plota.legend('SPY', 'blue', y)
y = plota.scale.volume(y)
plota(y, type = 'volume', plotX = F)
plota.legend('Volume', 'blue', Vo(y))
rsi = RSI(Cl(y),2)
plota(rsi, type = 'l', y.highlight = c(c(Inf,80),c(20,-Inf)))
abline(h = 20, col = 'red')
abline(h = 80, col = 'red')
plota.legend('RSI(2)', 'black', rsi)
dev.off()
y = data.spy['2010:01:01::2011:02:01']
png(filename = 'plot3.png', width = 500, height = 500, units = 'px', pointsize = 12, bg = 'white')
plota(y, type = 'ohlc', LeftMargin=3)
y0 = y;
y = data.ibm['2010:10:15::2011:02:01']
plota2Y(y, ylim = range(OHLC(y)),las=1, col='red', col.axis = 'red')
plota.ohlc(y, col = 'red')
plota.legend('SPY(rhs),IBM(lhs)', 'blue,red', list(y0,y))
dev.off()
y = data.spy['2010:01:01::2011:02:01']
png(filename = 'plot4.png', width = 500, height = 500, units = 'px', pointsize = 12, bg = 'white')
plota(y, type = 'candle')
y1 = to.monthly(y)
index(y1) = as.Date(index(y1))
plota.ohlc(y1, col = 'pink')
plota.candle(y)
plota.legend('Daily,Monthly', 'red,pink')
dev.off()
y = data.spy['2010:01:01::2011']
png(filename = 'plot5.png', width = 500, height = 500, units = 'px', pointsize = 12, bg = 'white')
layout(c(1,2,3))
plota(y, type = 'candle', plotX = F)
plota.legend('Daily', 'blue', y)
plota(y, ylim = range(OHLC(y)), plotX = F)
y1 = to.weekly(y)
index(y1) = as.Date(index(y1))
plota.candle(y1)
plota.legend('Weekly', 'blue', y1)
plota(y, ylim = range(OHLC(y)))
y1 = to.monthly(y)
index(y1) = as.Date(index(y1))
plota.candle(y1)
plota.legend('Monthly', 'blue', y1)
dev.off()
}
plota.colors <- function(N) {
if( is.list(N) ) N = len(N)
col = rev(c('yellow','cyan','magenta','red','gray','green','blue'))
temp = list()
for(j in 1:length(col)) {
temp[[j]] = colors()[grep(col[j],colors())]
temp[[j]] = temp[[j]][grep('^[^0-9]*$',temp[[j]])]
temp[[j]] = temp[[j]][order(nchar(temp[[j]]))]
index = which( colSums(col2rgb(temp[[j]])) < 100 )
if( length(index) > 0 ) temp[[j]] = temp[[j]][-index]
index = which( colSums(255 - col2rgb(temp[[j]])) < 100 )
if( length(index) > 0 ) temp[[j]] = temp[[j]][-index]
}
index = 1
col = rep('', N)
for(i in 1:10) {
for(j in 1:length(temp)) {
if(length(temp[[j]]) >= i) {
col[index] = temp[[j]][i]
index = index + 1
if(index > N) break
}
}
if(index > N) break
}
col
}
plota.colors.experimental = function(N) {
if( is.list(N) ) N = len(N)
library(RColorBrewer)
colorRampPalette(spl('red,green'))(N)
}
plota.stacked <- function
(
x,
y,
xlab='',
col = plota.colors(ncol(y)),
type=c('l','s'),
flip.legend = F,
...
)
{
y = 100 * y
y1 = list()
y1$positive = y
y1$positive[ y1$positive < 0 ] = 0
y1$negative = y
y1$negative[ y1$negative > 0 ] = 0
ylim = c(min(rowSums(y1$negative, na.rm = T)), max(1, rowSums(y1$positive, na.rm = T)))
if( class(x)[1] != 'Date' & class(x)[1] != 'POSIXct') {
plot(x, rep(0, len(x)), ylim = ylim, t = 'n', xlab = '', ylab = '', cex = par('cex'), ...)
grid()
} else {
plota(make.xts(y[,1], x), ylim = ylim, cex = par('cex'), LeftMargin = 4, ...)
axis(2, las = 1)
x = unclass(as.POSIXct(x))
}
mtext('Allocation %', side = 2,line = 3, cex = par('cex'))
mtext(xlab, side = 1,line = 2, cex = par('cex'))
if( type[1] == 'l' ) {
prep.x = c(x[1], x, x[len(x)])
for( y in y1 ) {
for (i in ncol(y) : 1) {
prep.y = c(0, rowSums(y[, 1 : i, drop = FALSE]), 0)
polygon(prep.x, prep.y, col = col[i], border = NA, angle = 90)
}
}
} else {
dx = mean(diff(x))
prep.x = c(rep(x,each=2), x[len(x)] + dx, x[len(x)] + dx)
for( y in y1 ) {
for (i in ncol(y) : 1) {
prep.y = c(0, rep(rowSums(y[, 1 : i, drop = FALSE]),each=2), 0)
polygon(prep.x, prep.y, col = col[i], border = NA, angle = 90)
}
}
}
if(flip.legend)
plota.legend(rev(colnames(y)), rev(col), cex = par('cex'))
else
plota.legend(colnames(y), col, cex = par('cex'))
}
plota.matplot <- function
(
y,
dates = NULL,
ylim = NULL,
type = 'l',
...
)
{
if( is.list(y) ) {
if(!is.null(dates)) y[[1]] = y[[1]][dates]
if(is.null(ylim)) {
ylim = c()
n = len(y)
for( i in 1:n ) {
if(!is.null(dates)) y[[i]] = y[[i]][dates]
ylim = range(ylim, y[[i]], na.rm = T)
}
}
plota(y[[1]], ylim = ylim, col = 1, type = type, ...)
if( n > 1 ) {
for( i in 2:n ) plota.lines(y[[i]], col = i, type = type, ...)
}
plota.legend(names(y), paste(1:n), y)
} else {
n = ncol(y)
if(!is.null(dates)) y = y[dates]
if(is.null(ylim)) ylim = range(y, na.rm = T)
plota(y[,1], ylim = ylim, col = 1, type = type, ...)
if( n > 1 ) {
for( i in 2:n ) plota.lines(y[,i], col = i, type = type, ...)
}
plota.legend(names(y), paste(1:n), as.list(y))
}
}
plota.add.copyright <- function(copyright = 'Systematic Investor') {
mtext(paste(copyright,'\uA9'), side = 1,line = -1, outer = T, adj = 1, font = 1, cex = 0.7, col='blue')
}
plota.recession <- function
(
col = col.add.alpha('gray', 50),	# Set alpha = 50 so it's relatively transparent
highlight = NULL
)
{
if( is.null(highlight) )
highlight = quantmod:::getSymbols('USREC', src = 'FRED', from = '1900-01-01', auto.assign = F)
plota.x.highlight(highlight, highlight != 0, col)
}
randsphere <- function
(
m,
n,
r
)
{
x = matrix(rnorm( m * n ), nrow = m, ncol = n);
s2 = apply(x^2, 1, sum)
return( x * repmat(r*(runif(m)^(1/n))/sqrt(s2),1,n) )
}
randsphere.test <- function()
{
load.packages('car,scatterplot3d')
png(filename = 'plot1.png', width = 500, height = 500, units = 'px', pointsize = 12, bg = 'white')
x = randsphere(1000, 2, 1)
y = x[, 2]
x = x[, 1]
par(mar = c(5,4,1,1))
plot(x,y, pch = 20)
ellipse(c(0, 0), matrix(c(1, 0, 0, 1), nrow = 2), radius = 1)
dev.off()
png(filename = 'plot2.png', width = 500, height = 500, units = 'px', pointsize = 12, bg = 'white')
layout(matrix(1:4,nrow=2))
x = randsphere(10000, 3, 1)
z = x[, 3]
y = x[, 2]
x = x[, 1]
scatterplot3d(x, y, z, highlight.3d = TRUE, pch = 20)
par(mar = c(5,4,1,1))
plot(x, y, pch = 20)
ellipse(c(0, 0), matrix(c(1, 0, 0, 1), nrow = 2), radius = 1)
plot(x, z, pch = 20)
ellipse(c(0, 0), matrix(c(1, 0, 0, 1), nrow = 2), radius = 1)
plot(y, z, pch = 20)
ellipse(c(0, 0), matrix(c(1, 0, 0, 1), nrow = 2), radius = 1)
dev.off()
}
randfixedsum <- function
(
m,
n,
s,
a,
b
)
{
if( (s<n*a) | (s>n*b) | (a>=b) )
stop('Inequalities n*a <= s <= n*b and a < b must hold.\n')
s = (s - n * a) / (b - a)
k = max( min( floor(s), n - 1), 0)
s = max( min( s, k + 1), k)
s1 = s - (k : (k - n + 1))
s2 = ((k + n) : (k+1)) - s
w = matrix(0, n, (n + 1))
realmax = 10^300
w[1, 2] = realmax
t = matrix(0, (n-1), n)
tiny = 2^(-1074)
for( i in 2:n ) {
tmp1 = w[(i - 1), 2 : (i + 1)] * s1[1 : i] / i
tmp2 = w[(i - 1), 1 : i] * s2[(n - i + 1) : n] / i
w[i, 2 : (i + 1)] = tmp1 + tmp2
tmp3 = w[i, 2 : (i+1)] + tiny
tmp4 = (s2[(n - i + 1) : n] > s1[1:i])
t[(i - 1), 1 : i] = (tmp2 / tmp3) * tmp4 + (1 - tmp1 / tmp3) * (!tmp4)
}
v = n^(3/2) * (w[n, (k + 2)] / realmax) * (b - a)^(n - 1)
x = matrix(0, n, m)
rt = matrix( runif((n-1) * m), (n-1), m)
rs = matrix( runif((n-1) * m), (n-1), m)
s = repmat(s, 1, m)
j = repmat((k + 1), 1, m)
sm = matrix(0, 1, m)
pr = matrix(1, 1, m)
for( i in (n - 1):1) {
e = (rt[(n - i), ] <= t[i, j])
sx = rs[(n - i), ]^(1/i)
sm = sm + (1 - sx) * pr * s / (i+1)
pr = sx * pr
x[(n - i), ] = sm + pr * e
s = s - e
j = j - e
}
x[n, ] = sm + pr * s
rp = matrix( runif(n * m), n, m)
p = apply(rp, 2, order)
x = (b - a) * x[p + repmat(t(seq(0, n * (m - 1), by = n)), n, 1)] + a
x = matrix(x, n, m)
return(t(x))
}
randfixedsum.test <- function()
{
load.packages('car,scatterplot3d')
png(filename = 'plot1.png', width = 500, height = 500, units = 'px', pointsize = 12, bg = 'white')
x = randfixedsum(100, 2, 1, 0.2, 0.8)
y = x[, 2]
x = x[, 1]
par(mar = c(5,4,1,1))
plot(x,y, pch = 20)
dev.off()
png(filename = 'plot2.png', width = 500, height = 500, units = 'px', pointsize = 12, bg = 'white')
layout(matrix(1:4,nrow=2))
x = randfixedsum(1000, 3, 1, 0.2, 0.8)
z = x[, 3]
y = x[, 2]
x = x[, 1]
scatterplot3d(x, y, z, highlight.3d = TRUE, pch = 20, angle=190)
par(mar = c(5,4,1,1))
plot(x, y, pch = 20)
plot(x, z, pch = 20)
plot(y, z, pch = 20)
dev.off()
}
asset.paths <- function(s0, mu, sigma,
nsims = 10000,
periods = c(0, 1)
)
{
s0 = as.vector(s0)
nsteps = len(periods)
dt = c(periods[1], diff(periods))
if( len(s0) == 1 ) {
drift = mu - 0.5 * sigma^2
if( nsteps == 1 ) {
s0 * exp(drift * dt + sigma * sqrt(dt) * rnorm(nsims))
} else {
temp = matrix(exp(drift * dt + sigma * sqrt(dt) * rnorm(nsteps * nsims)), nc=nsims)
for(i in 2:nsteps) temp[i,] = temp[i,] * temp[(i-1),]
s0 * temp
}
} else {
require(MASS)
drift = mu - 0.5 * diag(sigma)
n = len(mu)
if( nsteps == 1 ) {
s0 * exp(drift * dt + sqrt(dt) * t(mvrnorm(nsims, rep(0, n), sigma)))
} else {
temp = array(exp(as.vector(drift %*% t(dt)) + t(sqrt(dt) * mvrnorm(nsteps * nsims, rep(0, n), sigma))), c(n, nsteps, nsims))
for(i in 2:nsteps) temp[,i,] = temp[,i,] * temp[,(i-1),]
s0 * temp
}
}
}
asset.paths.test <- function()
{
S = c(100,105)
X = 98
Time = 0.5
r = 0.05
sigma = c(0.11,0.16)
rho = 0.63
N = 10000
periods = 0:10
prices = asset.paths(S[1], r, sigma[1], N, periods = periods)
png(filename = 'plot1.png', width = 600, height = 600, units = 'px', pointsize = 12, bg = 'white')
matplot(prices[,1:100], type='l', xlab='Years', ylab='Prices',
main='Selected Price Paths')
dev.off()
periods = 0:10
cov.matrix = sigma%*%t(sigma) * matrix(c(1,rho,rho,1),2,2)
prices = asset.paths(S, c(r,r), cov.matrix, N, periods = periods)
png(filename = 'plot2.png', width = 600, height = 600, units = 'px', pointsize = 12, bg = 'white')
layout(1:2)
matplot(prices[1,,1:100], type='l', xlab='Years', ylab='Prices',
main='Selected Price Paths for Asset 1')
matplot(prices[2,,1:100], type='l', xlab='Years', ylab='Prices',
main='Selected Price Paths for Asset 2')
dev.off()
cor(as.vector(prices[1,,] / mlag(prices[1,,])),
as.vector(prices[2,,] / mlag(prices[2,,])),
use='complete.obs', method='pearson')
load.packages('fOptions')
GBSOption(TypeFlag = "c", S = S[1], X = X, Time = Time, r = r, b = r, sigma = sigma[1])
N = 1000000
prices = asset.paths(S[1], r, sigma[1], N, periods = Time)
future.payoff = pmax(0, prices - X)
discounted.payoff = future.payoff * exp(-r * Time)
mean(discounted.payoff)
load.packages('fExoticOptions')
Time = 1/12
TurnbullWakemanAsianApproxOption(TypeFlag = "c", S = S[1], SA = S[1],
X = X, Time = Time, time = Time, tau = 0 , r = r, b = r, sigma = sigma[1])
N = 100000
periods = seq(0,Time,1/360)
n = len(periods)
prices = asset.paths(S[1], r, sigma[1], N, periods = periods)
future.payoff = pmax(0, colSums(prices)/n - X)
discounted.payoff = future.payoff * exp(-r * Time)
mean(discounted.payoff)
Time = 0.5
TwoRiskyAssetsOption(TypeFlag = "cmax", S1 = S[1], S2 = S[2],
X = X, Time = Time, r = r, b1 = r, b2 = r,
sigma1 = sigma[1], sigma2 = sigma[2], rho = rho)
N = 100000
cov.matrix = sigma%*%t(sigma) * matrix(c(1,rho,rho,1),2,2)
prices = asset.paths(S, c(r,r), sigma = cov.matrix, N, periods = Time)
future.payoff = pmax(0, apply(prices,2,max) - X)
discounted.payoff = future.payoff * exp(-r * Time)
mean(discounted.payoff)
Time = 1/12
N = 10000
periods = seq(0,Time,1/360)
n = len(periods)
prices = asset.paths(S, c(r,r), sigma = cov.matrix, N, periods = periods)
future.payoff = pmax(0, colSums(apply(prices,c(2,3),max))/n - X)
discounted.payoff = future.payoff * exp(-r * Time)
mean(discounted.payoff)
}
mat.3d = function(mat, n3dim) { array(mat, c(dim(mat), n3dim)) }
mat.slice.3d = function(mat, slice.at) { aperm(array(t(mat), c(ncol(mat), nrow(mat) / slice.at, slice.at)),c(2,1,3)) }
asset.paths.at.period.ends <- function(
prices,
period.ends,
nsims = 100,
lookback.len = NA
) {
load.packages('MASS')
ret = prices / mlag(prices) - 1
n = ncol(prices)
nperiods = nrow(prices)
index = 1:nperiods
first.index = sapply(1:n, function(i) index[!is.na(prices[,i])][1])
first.value = sapply(1:n, function(i) prices[first.index[i],i])
s0 = rep(1,n)
scenarios = array(NA, c(nperiods, n, nsims))
period.ends.all = sort(unique(c(0, period.ends, len(period.ends))))
for(i in 2:len(period.ends.all)) {
index.hist = index = (period.ends.all[i-1] + 1) : period.ends.all[i]
n.index = len(index)
if(!is.na(lookback.len))
index.hist = max(1, period.ends.all[i] - lookback.len + 1 ) : period.ends.all[i]
hist = ret[index.hist,, drop=F]
include.index = colSums(!is.na(hist)) != 0
if(!any(include.index)) continue
hist = hist[,include.index, drop=F]
mu = apply(hist, 2, mean, na.rm=T)
risk = apply(hist, 2, sd, na.rm=T)
correlation = cor(hist, use = 'complete.obs', method = 'pearson')
cov = correlation * (risk %*% t(risk))
scenarios[index,include.index,] = mat.slice.3d(mvrnorm(n.index * nsims, mu, cov), nsims)
}
scenarios[1,,] = 0
for(i in 1:n)
for(j in 1:nsims) {
if(first.index[i] > 1) scenarios[first.index[i],i,j] = 0
scenarios[first.index[i]:nperiods,i,j] = first.value[i] * cumprod(1 + scenarios[first.index[i]:nperiods,i,j])
}
scenarios = scenarios * mat.3d(1+0*prices, nsims)
scenarios
}
all.permutations <- function(n = 1) {
m = matrix(F,2^n,n)
m[2,1] = T
if (n == 1) return(m)
istart = 2
for(i in 2:n) {
index = (istart+1):(2*istart)
m[index, ] = m[1:istart,]
m[index, i] = T
istart = istart * 2
}
return(m)
}
seasonality.test <- function()
{
load.packages('quantmod')
ticker = 'WMT'
data = getSymbols(ticker, src = 'yahoo', from = '1970-01-01', auto.assign = F)
data = adjustOHLC(data, use.Adjusted=T)
data = data['1980::2012:04:07']
png(filename = 'plot.month.year.seasonality.png', width = 600, height = 500, units = 'px', pointsize = 12, bg = 'white')
month.year.seasonality(data, ticker)
dev.off()
month.ends = endpoints(data, 'months')
month.ends = month.ends[month.ends > 0 & month.ends < nrow(data)]
index = which(format(index(data), '%b')[month.ends] == 'Apr')
png(filename = 'plot.time.seasonality.png', width = 600, height = 500, units = 'px', pointsize = 12, bg = 'white')
layout(1)
time.seasonality(data, 1 + month.ends[index], 20, ticker)
dev.off()
}
pattern.test <- function()
{
load.packages('quantmod')
ticker = 'SPY'
data = getSymbols(ticker, src = 'yahoo', from = '1970-01-01', auto.assign = F)
data = adjustOHLC(data, use.Adjusted=T)
data = data['::2012:04:07']
png(filename = 'plot1.time.series.pattern.matching.png', width = 600, height = 500, units = 'px', pointsize = 12, bg = 'white')
matches = bt.matching.find(Cl(data), main = ticker,	n.query=90, plot=TRUE)
dev.off()
png(filename = 'plot2.time.series.pattern.matching.png', width = 600, height = 500, units = 'px', pointsize = 12, bg = 'white')
out = bt.matching.overlay(matches, plot=TRUE)
dev.off()
png(filename = 'plot.pattern.matching.png', width = 600, height = 500, units = 'px', pointsize = 12, bg = 'white')
plot.patterns(data, 190, ticker)
dev.off()
}
month.year.seasonality <- function
(
data,
ticker,
lookback.len = 20*252
)
{
data = last(data, lookback.len)
nperiods = nrow(data)
month.ends = endpoints(data, 'months')
month.ends = month.ends[month.ends > 0]
prices = Cl(data)[month.ends]
ret = prices / mlag(prices) - 1
ret = ret[-1]
ret.by.month = create.monthly.table(ret)
data_list = lapply(apply(ret.by.month, 2, list), '[[', 1)
group.seasonality(data_list, paste(ticker, 'Monthly', join(format(range(index(data)), '%d-%b-%Y'), ' to\n')))
}
group.seasonality <- function
(
data_list,
smain,
...
)
{
out = compute.stats( data_list,
list(Sharpe=function(x) mean(x,na.rm=T)/sd(x,na.rm=T),
'% Positive'=function(x) sum(x > 0,na.rm=T)/sum(!is.na(x)),
Min=function(x) min(x,na.rm=T),
Max=function(x) max(x,na.rm=T),
Avg=function(x) mean(x,na.rm=T),
Med=function(x) median(x,na.rm=T),
StDev=function(x) sd(x,na.rm=T)
)
)
layout(mat=matrix(1:4, 2, 2, byrow=FALSE))
par(mar=c(4, 3, 2, 2))
col = spl('lightgray,red')
stats.names = spl('Sharpe,% Positive,Min,Avg')
for(i in stats.names) {
barplot(100*out[i,], names.arg = colnames(out),
col=iif(out[i,] > 0, col[1], col[2]),
main=iif(i == stats.names[1], paste(smain,' ', i, sep=''), i),
border = 'darkgray',las=2, ...)
grid(NA,NULL)
abline(h=0, col='black')
}
}
time.seasonality <- function
(
data,
period.starts,
period.len,
ticker
)
{
nperiods = len(period.starts)
dates = index(data)
ndates = len(dates)
prices = Cl(data)
ret = prices / mlag(prices) - 1
ret = c( as.double(ret), rep(NA, period.len))
trading.days = sapply(period.starts, function(i) ret[i : (i + period.len - 1)])
periods.index = 1:nperiods
temp = count(trading.days)
periods.index = periods.index[temp > 0.9 * median(temp)]
last.period = trading.days[, nperiods]
last.period = 100 * ( cumprod(1 + last.period) - 1 )
avg.period = apply(trading.days[, periods.index], 1, mean, na.rm=T)
avg.period = 100 * ( cumprod(1 + avg.period) - 1 )
temp = 100*(apply(1 + trading.days[, periods.index], 2, cumprod) - 1)
quantiles = apply(temp, 1, quantile, probs = c(75, 25)/100, na.rm=T)
cols = spl('blue,red,gray')
par(mar=c(4,4,2,1))
plot(avg.period, type='n', xaxt = 'n', xlim=c(1,period.len),
ylim=range(avg.period, last.period, quantiles, na.rm=T),
main = ticker, xlab = 'Trading Days', ylab = 'Avg % Profit/Loss')
grid()
axis(1, 1:period.len)
lines(quantiles[1,], type='l', lwd=2, col=cols[3])
lines(quantiles[2,], type='l', lwd=2, col=cols[3])
lines(last.period, type='b', lwd=2, col=cols[2], bg=cols[2], pch=24)
lines(avg.period, type='b', lwd=2, col=cols[1], bg=cols[1], pch=22)
first.year = format(dates[period.starts][periods.index[1]], '%Y')
last.year = format(dates[period.starts][last(periods.index)], '%Y')
last.period.start.date = format(dates[period.starts[nperiods]], '%d %b %Y')
last.period.end.date = format(dates[ndates], '%d %b %Y')
if( (period.starts[nperiods] + period.len - 1) < ndates ) {
last.period.end.date = format(dates[period.starts[nperiods] + period.len - 1], '%d %b %Y')
}
plota.legend(c(paste('Avgerage for', first.year, '-', last.year),
paste(last.period.start.date, '-', last.period.end.date),
'Top 25% / Bot 25%'), cols)
}
plot.patterns <- function
(
data,
n,
ticker,
patterns = pattern.db()
)
{
load.packages('sm')
sample = last(data, n)
obj = find.extrema( Cl(sample) )
mhat = obj$mhat
mhat.extrema.loc = obj$mhat.extrema.loc
data.extrema.loc = obj$data.extrema.loc
n.index = len(data.extrema.loc)
plota.control$col.border = 'gray'
plota(sample, type='hl',col='gray')
plota.lines(mhat, col='magenta', lwd=2)
plota.lines(sample, col='blue')
if(n.index > 0) {
plota.lines(sample[data.extrema.loc], type='p', col='blue', lwd=3, pch=19)
out = find.patterns(obj, patterns = patterns, silent=F, plot=T)
}
plota.legend(c(paste(ticker, join(format(range(index(sample)), '%d%b%Y'), ' - ')),
'Close,Kernel,Pattern(s)'),
'gray,blue,magenta,orange')
}
find.extrema <- function(
x
)
{
if(is.xts(x)) {
y = as.vector( Cl(x) )
} else {
y = x
}
n = len(y)
t = 1:n
h = h.select(t, y, method = 'cv')
temp = sm.regression(t, y, h=h, display = 'none')
mhat = approx(temp$eval.points, temp$estimate, t, method='linear')$y
temp = diff(sign(diff(mhat)))
loc = which( temp != 0 ) + 1
loc.dir = -sign(temp[(loc - 1)])
temp = c( y[1], y, y[n] )
temp = cbind(temp[loc], temp[(loc + 1)], temp[(loc + 2)])
max.index = loc + apply(temp, 1, which.max) - 2
min.index = loc + apply(temp, 1, which.min) - 2
data.loc = iif(loc.dir > 0, max.index, min.index)
data.loc = iif(data.loc < 1, 1, iif(data.loc > n, n, data.loc))
if(is.xts(x)) mhat = make.xts(mhat, index(x))
return(list(data = y, mhat = mhat, extrema.dir = loc.dir,
mhat.extrema.loc = loc, data.extrema.loc = data.loc))
}
find.patterns <- function
(
obj,
patterns = pattern.db(),
silent=T,
plot=T
)
{
data = obj$data
mhat = obj$mhat
extrema.dir = obj$extrema.dir
data.extrema.loc = obj$data.extrema.loc
n.index = len(data.extrema.loc)
if(is.xts(mhat)) {
dates = index4xts(obj$mhat)
} else {
dates = 1:len(data)
}
col = col.add.alpha('orange', alpha=150)
out = out.rownames = c()
for(i in 1:n.index) {
for(pattern in patterns) {
if( pattern$start * extrema.dir[i] > 0 ) {
if( i + pattern$len - 1 <= n.index ) {
envir.data = c(data[data.extrema.loc][i:(i + pattern$len - 1)],
data.extrema.loc[i:(i + pattern$len - 1)])
names(envir.data) = c(paste('E', 1:pattern$len, sep=''),
paste('t', 1:pattern$len, sep=''))
envir.data = as.list(envir.data)
envir.data$E = data[data.extrema.loc][-c(1:i)]
envir.data$t = data.extrema.loc[-c(1:i)]
if( eval(pattern$formula, envir = envir.data) ) {
if(!silent) cat('Found', pattern$name, 'at', i, '\n')
if(plot & !is.null(pattern$plot)) {
temp = dates[data.extrema.loc[i:(i + pattern$len - 1)]]
names(temp) = paste('d', 1:pattern$len, sep='')
envir.data = c( envir.data, temp )
envir.data$d = dates[data.extrema.loc[-c(1:i)]]
envir.data$col = col
eval(pattern$plot, envir = envir.data)
}
out.rownames = c(out.rownames, pattern$name)
out = rbind(out, c(data.extrema.loc[i],
iif(is.null(pattern$last.point),
data.extrema.loc[(i + pattern$len - 1)],
eval(pattern$last.point, envir = envir.data)
)))
}
}
}
}
}
if(len(out)>0) {
colnames(out) = spl('start,end')
rownames(out) = out.rownames
}
return(out)
}
find.all.patterns.window <- function()
{
load.packages('quantmod')
ticker = 'SPY'
data = getSymbols(ticker, src = 'yahoo', from = '1970-01-01', auto.assign = F)
data = adjustOHLC(data, use.Adjusted=T)
data = data['2010::']
load.packages('sm')
history = as.vector(coredata(Cl(data)))
window.len = 90
patterns = pattern.db()
found.patterns = c()
for(t in window.len : (len(history)-1)) {
sample = history[(t - window.len + 1):t]
obj = find.extrema( sample )
if(len(obj$data.extrema.loc) > 0) {
out =  find.patterns(obj, patterns = patterns, silent=F, plot=F)
if(len(out)>0) found.patterns = rbind(found.patterns,cbind(t,out,t-window.len+out))
}
if( t %% 10 == 0) cat(t, 'out of', len(history), '\n')
}
colnames(found.patterns) = spl('t,start,end,tstart,tend')
setdiff(unique(names(patterns)), unique(rownames(found.patterns)))
frequency = tapply(rep(1,nrow(found.patterns)), rownames(found.patterns), sum)
barplot(frequency)
index = which(rownames(found.patterns)=='HS')
found.patterns[ index[1:10],]
pattern.index = index[1]
t = found.patterns[pattern.index, 't'];
plot.patterns(data[1:t,], window.len, ticker, patterns)
sample = data[(t - window.len + 1):t]
start.pattern = found.patterns[pattern.index, 'start']
end.pattern = found.patterns[pattern.index, 'end']
abline(v = index(sample)[start.pattern]);
abline(v = index(sample)[end.pattern]);
index(sample)[start.pattern]
index(data[(t - window.len + start.pattern),])
}
bt.patterns.test <- function()
{
load.packages('quantmod')
ticker = 'SPY'
data = getSymbols(ticker, src = 'yahoo', from = '1970-01-01', auto.assign = F)
data = adjustOHLC(data, use.Adjusted=T)
load.packages('sm')
history = as.vector(coredata(Cl(data)))
window.L = 35
window.d = 3
window.len = window.L + window.d
patterns = pattern.db()
found.patterns = c()
for(t in window.len : (len(history)-1)) {
ret = history[(t+1)]/history[t]-1
sample = history[(t - window.len + 1):t]
obj = find.extrema( sample )
if(len(obj$data.extrema.loc) > 0) {
out =  find.patterns(obj, patterns = patterns, silent=F, plot=F)
if(len(out)>0) found.patterns = rbind(found.patterns,cbind(t,out,t-window.len+out, ret))
}
if( t %% 10 == 0) cat(t, 'out of', len(history), '\n')
}
colnames(found.patterns) = spl('t,start,end,tstart,tend,ret')
found.patterns = found.patterns[found.patterns[,'end'] <= window.L,]
pattern.names = unique(rownames(found.patterns))
all.patterns = c()
for(name in pattern.names) {
index = which(rownames(found.patterns) == name)
temp = NA * found.patterns[index,]
i.count = 0
i.start = 1
while(i.start < len(index)) {
i.count = i.count + 1
temp[i.count,] = found.patterns[index[i.start],]
subindex = which(found.patterns[index,'tstart'] > temp[i.count,'tend'])
if(len(subindex) > 0) {
i.start = subindex[1]
} else break
}
all.patterns = rbind(all.patterns, temp[1:i.count,])
}
png(filename = 'plot1.png', width = 600, height = 500, units = 'px', pointsize = 12, bg = 'white')
frequency = tapply(rep(1,nrow(all.patterns)), rownames(all.patterns), sum)
layout(1)
barplot.with.labels(frequency/100, 'Frequency for each Pattern')
dev.off()
png(filename = 'plot2.png', width = 600, height = 500, units = 'px', pointsize = 12, bg = 'white')
all.patterns[,'ret'] = history[(all.patterns[,'t']+20)] / history[all.patterns[,'t']] - 1
data_list = tapply(all.patterns[,'ret'], rownames(all.patterns), list)
group.seasonality(data_list, '20 days after Pattern')
dev.off()
png(filename = 'plot3.png', width = 600, height = 500, units = 'px', pointsize = 12, bg = 'white')
layout(1)
name = 'BBOT'
index = which(rownames(all.patterns) == name)
time.seasonality(data, all.patterns[index,'t'], 20, name)
dev.off()
}
pattern.db <- function()
{
patterns = list()
pattern = list()
pattern$len = 5
pattern$start = 'max'
pattern$formula = expression(
E3 > E1 &
E3 > E5 &
abs(E1 - (E1+E5)/2) < 1.5/100 * (E1+E5)/2 &
abs(E5 - (E1+E5)/2) < 1.5/100 * (E1+E5)/2 &
abs(E2 - (E2+E4)/2) < 1.5/100 * (E2+E4)/2 &
abs(E4 - (E2+E4)/2) < 1.5/100 * (E2+E4)/2
)
pattern$plot = expression({
lines(c(d1,d2,d3,d4,d5), c(E1,E2,E3,E4,E5), lwd=10, col=col)
text(d3, E3, 'HS', adj=c(0.5,-0.5), xpd=TRUE)
})
patterns$HS = pattern
pattern = list()
pattern$len = 5
pattern$start = 'min'
pattern$formula = expression(
E3 < E1 &
E3 < E5 &
abs(E1 - (E1+E5)/2) < 1.5/100 * (E1+E5)/2 &
abs(E5 - (E1+E5)/2) < 1.5/100 * (E1+E5)/2 &
abs(E2 - (E2+E4)/2) < 1.5/100 * (E2+E4)/2 &
abs(E4 - (E2+E4)/2) < 1.5/100 * (E2+E4)/2
)
pattern$plot = expression({
lines(c(d1,d2,d3,d4,d5), c(E1,E2,E3,E4,E5), lwd=10, col=col)
text(d3, E3, 'IHS', adj=c(0.5,1), xpd=TRUE)
})
patterns$IHS = pattern
pattern = list()
pattern$len = 5
pattern$start = 'max'
pattern$formula = expression(
E1 < E3 &
E3 < E5 &
E2 > E4
)
pattern$plot = expression({
beta = ols(cbind(1,c(t1,t3,t5)),c(E1,E3,E5))$coefficients
lines(c(d1,d3,d5), beta[1] + beta[2]*c(t1,t3,t5), lwd=10, col=col)
lines(c(d2,d4), c(E2,E4), lwd=10, col=col)
text(d3, min(E2,E4), 'BTOP', adj=c(0.5,1), xpd=TRUE)
})
patterns$BTOP = pattern
pattern = list()
pattern$len = 5
pattern$start = 'min'
pattern$formula = expression(
E1 > E3 &
E3 > E5 &
E2 < E4
)
pattern$plot = expression({
beta = ols(cbind(1,c(t1,t3,t5)),c(E1,E3,E5))$coefficients
lines(c(d1,d3,d5), beta[1] + beta[2]*c(t1,t3,t5), lwd=10, col=col)
lines(c(d2,d4), c(E2,E4), lwd=10, col=col)
text(d3, max(E2,E4), 'BBOT', adj=c(0.5,0), xpd=TRUE)
})
patterns$BBOT = pattern
pattern = list()
pattern$len = 5
pattern$start = 'max'
pattern$formula = expression(
E1 > E3 &
E3 > E5 &
E2 < E4
)
pattern$plot = expression({
beta = ols(cbind(1,c(t1,t3,t5)),c(E1,E3,E5))$coefficients
lines(c(d1,d3,d5), beta[1] + beta[2]*c(t1,t3,t5), lwd=10, col=col)
lines(c(d2,d4), c(E2,E4), lwd=10, col=col)
text(d3, min(E2,E4), 'TTOP', adj=c(0.5,1), xpd=TRUE)
})
patterns$TTOP = pattern
pattern = list()
pattern$len = 5
pattern$start = 'min'
pattern$formula = expression(
E1 < E3 &
E3 < E5 &
E2 > E4
)
pattern$plot = expression({
beta = ols(cbind(1,c(t1,t3,t5)),c(E1,E3,E5))$coefficients
lines(c(d1,d3,d5), beta[1] + beta[2]*c(t1,t3,t5), lwd=10, col=col)
lines(c(d2,d4), c(E2,E4), lwd=10, col=col)
text(d3, max(E2,E4), 'TBOT', adj=c(0.5,0), xpd=TRUE)
})
patterns$TBOT = pattern
pattern = list()
pattern$len = 5
pattern$start = 'max'
pattern$formula = expression({
avg.top = (E1+E3+E5)/3
avg.bop = (E2+E4)/2
abs(E1 - avg.top) < 0.75/100 * avg.top &
abs(E3 - avg.top) < 0.75/100 * avg.top &
abs(E5 - avg.top) < 0.75/100 * avg.top &
abs(E2 - avg.bop) < 0.75/100 * avg.bop &
abs(E4 - avg.bop) < 0.75/100 * avg.bop &
min(E1,E3,E5) > max(E2,E4)
})
pattern$plot = expression({
avg.top = (E1+E3+E5)/3
avg.bop = (E2+E4)/2
lines(c(d1,d3,d5), rep(avg.top,3), lwd=10, col=col)
lines(c(d2,d4), rep(avg.bop,2), lwd=10, col=col)
text(d3, min(E2,E4), 'RTOP', adj=c(0.5,-0.5), xpd=TRUE)
})
patterns$RTOP = pattern
pattern = list()
pattern$len = 5
pattern$start = 'min'
pattern$formula = expression({
avg.top = (E2+E4)/2
avg.bop = (E1+E3+E5)/3
abs(E2 - avg.top) < 0.75/100 * avg.top &
abs(E4 - avg.top) < 0.75/100 * avg.top &
abs(E1 - avg.bop) < 0.75/100 * avg.bop &
abs(E3 - avg.bop) < 0.75/100 * avg.bop &
abs(E5 - avg.bop) < 0.75/100 * avg.bop &
min(E2,E4) > max(E1,E3,E5)
})
pattern$plot = expression({
avg.top = (E2+E4)/2
avg.bop = (E1+E3+E5)/3
lines(c(d1,d3,d5), rep(avg.bop,3), lwd=10, col=col)
lines(c(d2,d4), rep(avg.top,2), lwd=10, col=col)
text(d3, max(E2,E4), 'RBOT', adj=c(0.5,0), xpd=TRUE)
})
patterns$RBOT = pattern
pattern = list()
pattern$len = 3
pattern$start = 'max'
pattern$formula = expression({
second.top = max(E)
second.top.t = t[which.max(E)]
avg = (E1 + second.top)/2
abs(E1         - avg) < 1.5/100 * avg &
abs(second.top - avg) < 1.5/100 * avg &
second.top.t - t1 > 22
})
pattern$plot = expression({
second.top = max(E)
second.top.d = d[which.max(E)]
avg = (E1 + second.top)/2
points(c(d1, second.top.d), c(E1, second.top), pch=2, lwd=2)
lines(c(d1, second.top.d), rep(avg, 2), lwd=10, col=col)
text(d2, avg, 'DTOP', adj=c(0.5,-0.5), xpd=TRUE)
})
pattern$last.point = expression(t[which.max(E)])
patterns$DTOP = pattern
pattern = list()
pattern$len = 3
pattern$start = 'min'
pattern$formula = expression(
abs(E1 -         (E1+min(E))/2) < 1.5/100 * (E1+min(E))/2 &
abs(max(E[-1]) - (E1+min(E))/2) < 1.5/100 * (E1+min(E))/2 &
t[which.min(E)] - t1 > 22
)
pattern$plot = expression({
second.bot = min(E)
second.bot.d = d[which.min(E)]
avg = (E1 + second.bot)/2
points(c(d1, second.bot.d), c(E1, second.bot), pch=2, lwd=2)
lines(c(d1, second.bot.d), rep(avg, 2), lwd=10, col=col)
text(d2, avg, 'DBOT', adj=c(0.5,1), xpd=TRUE)
})
pattern$last.point = expression(t[which.min(E)])
patterns$DBOT = pattern
for(i in 1:len(patterns)) {
patterns[[i]]$name = names(patterns)[i]
patterns[[i]]$start = iif(patterns[[i]]$start == 'max', 1, -1)
}
return(patterns)
}
bt.cluster.risk.parity.weights.test <- function()
{
load.packages('quantmod')
tickers = spl('GLD,UUP,SPY,QQQ,IWM,EEM,EFA,IYR,USO,TLT')
map = spl('Gold GLD,US Dollar UUP,S&P500 SPY,Nasdaq QQQ,Small Cap IWM,EmergingM EEM,InternationalM EFA,Real Estate IYR,Oil USO,Treasurys TLT')
names(map) = tickers
data <- new.env()
getSymbols(tickers, src = 'yahoo', from = '1900-01-01', env = data, auto.assign = T)
for(i in ls(data)) data[[i]] = adjustOHLC(data[[i]], use.Adjusted=T)
for(i in ls(data)) data[[ map[i] ]] = data[[i]]
rm(list=tickers, envir=data)
bt.prep(data, align='remove.na', dates = '2011:12::2012')
periodicity = 'months'
lookback.len = 250
cluster.group <<- cluster.group.kmeans.90
obj = portfolio.allocation.helper(data$prices,
periodicity = periodicity, lookback.len = lookback.len,
min.risk.fns = list(
EW=equal.weight.portfolio,
RP=risk.parity.portfolio(),
ERC=equal.risk.contribution.portfolio,
G2.EW = distribute.weights(equal.weight.portfolio, cluster.group),
G2.RP=distribute.weights(risk.parity.portfolio(), cluster.group),
G2.MV=distribute.weights(min.var.portfolio, cluster.group),
G2.ERC=distribute.weights(equal.risk.contribution.portfolio, cluster.group)
),
adjust2positive.definite = F,
custom.stats.fn = portfolio.allocation.custom.stats.clusters
)
clusters = coredata(obj$clusters$EW)[13,]
temp = clusters
temp[] = 0
temp[clusters == clusters[names(clusters) == 'Treasurys TLT']] = 1
temp[clusters == clusters[names(clusters) == 'US Dollar UUP']] = 2
temp[clusters == clusters[names(clusters) == 'EmergingM EEM']] = 3
temp[clusters == clusters[names(clusters) == 'Gold GLD']] = 4
clusters = temp
png(filename = 'plot1.png', width = 1200, height = 600, units = 'px', pointsize = 12, bg = 'white')
layout(matrix(1:2,nc=2))
plot.cluster.weights(coredata(obj$weights$ERC)[13,], clusters,
main='ERC Weights')
plot.cluster.weights(coredata(obj$risk.contributions$ERC)[13,], clusters,
main='ERC Risk Contributions')
dev.off()
png(filename = 'plot2.png', width = 1200, height = 600, units = 'px', pointsize = 12, bg = 'white')
layout(matrix(1:2,nc=2))
plot.cluster.weights(coredata(obj$weights$G2.ERC)[13,], clusters,
main='Cluster ERC Weights')
plot.cluster.weights(coredata(obj$risk.contributions$G2.ERC)[13,], clusters,
main='Cluster ERC Risk Contributions')
dev.off()
}
bt.cluster.risk.parity.10.major.assets <- function()
{
tickers = spl('SPY,EFA,EWJ,EEM,IYR,RWX,IEF,TLT,DBC,GLD')
dates='2004:12::'
name = 'ETFs AAA'
load.packages('quantmod')
data <- new.env()
getSymbols(tickers, src = 'yahoo', from = '2000-01-01', env = data, auto.assign = T)
for(i in ls(data)) data[[i]] = adjustOHLC(data[[i]], use.Adjusted=T)
bt.prep(data, align='keep.all', dates=dates, fill.gaps=T)
periodicity = 'weeks'
lookback.len = 250
cluster.group <<- cluster.group.kmeans.90
obj = portfolio.allocation.helper(data$prices,
periodicity = periodicity, lookback.len = lookback.len,
min.risk.fns = list(
EW=equal.weight.portfolio,
RP=risk.parity.portfolio(),
ERC=equal.risk.contribution.portfolio,
Dynamic.EW = distribute.weights(equal.weight.portfolio, cluster.group),
Dynamic.RP=distribute.weights(risk.parity.portfolio(), cluster.group),
Dynamic.ERC=distribute.weights(equal.risk.contribution.portfolio, cluster.group)
),
adjust2positive.definite = F,
custom.stats.fn = portfolio.allocation.custom.stats
)
models = create.strategies(obj, data, dates=(lookback.len):nrow(data$prices))$models
title = paste(name, '(' ,periodicity, ',' , lookback.len, 'days )')
stats = bt.summary.report(models, title, data, obj,
control = list(
plot.weight.transition.maps = F,
plot.risk.contribution.transition.maps = F)
)
}
bt.cluster.risk.parity.dow.30 <- function()
{
load.packages('quantmod')
dates='1995::'
name = 'Dow Jones 30'
data = load.dow.jones(align='keep.all', dates=dates)
sectors = data$sectors
tickers = data$symbolnames
periodicity = 'weeks'
lookback.len = 250
cluster.group <<- cluster.group.kmeans.90
obj = portfolio.allocation.helper(data$prices,
periodicity = periodicity, lookback.len = lookback.len,
min.risk.fns = list(EW=equal.weight.portfolio,
RP=risk.parity.portfolio(),
ERC=equal.risk.contribution.portfolio,
Static.EW = distribute.weights(equal.weight.portfolio, static.group(as.numeric(sectors))),
Static.RP=distribute.weights(risk.parity.portfolio(), static.group(sectors)),
Static.ERC=distribute.weights(equal.risk.contribution.portfolio, static.group(sectors)),
Dynamic.EW = distribute.weights(equal.weight.portfolio, cluster.group),
Dynamic.RP=distribute.weights(risk.parity.portfolio(), cluster.group),
Dynamic.ERC=distribute.weights(equal.risk.contribution.portfolio, cluster.group)
),
adjust2positive.definite = F,
custom.stats.fn = portfolio.allocation.custom.stats
)
models = create.strategies(obj, data, dates=(lookback.len):nrow(data$prices))$models
title = paste(name, '(' ,periodicity, ',' , lookback.len, 'days )')
stats = bt.summary.report(models, title, data, obj,
control = list(
plot.weight.transition.maps = F,
plot.risk.contribution.transition.maps = F)
)
}
load.dow.jones <- function(align='remove.na', dates = NULL) {
tickers = spl('XLY,XLP,XLE,XLF,XLV,XLI,XLB,XLK,XLU')
tickers.desc = spl('ConsumerCyclicals,ConsumerStaples,Energy,Financials,HealthCare,Industrials,Materials,Technology,Utilities')
sector.map = c()
for(i in 1:len(tickers)) {
sector.map = rbind(sector.map,
cbind(sector.spdr.components(tickers[i]), tickers.desc[i])
)
}
colnames(sector.map) = spl('ticker,sector')
load.packages('quantmod')
tickers = dow.jones.components()
sectors = factor(sector.map[ match(tickers, sector.map[,'ticker']), 'sector'])
names(sectors) = tickers
data <- new.env()
getSymbols(tickers, src = 'yahoo', from = '1900-01-01', env = data, auto.assign = T)
for(i in ls(data)) data[[i]] = adjustOHLC(data[[i]], use.Adjusted=T)
bt.prep(data, align=align,dates=dates)
data$sectors = sectors[data$symbolnames]
return(data)
}
portfolio.allocation.custom.stats.clusters <- function(x,ia) {
risk.contributions = portfolio.risk.contribution(x, ia)
clusters = cluster.group(ia)
return(list(
risk.contributions = risk.contributions,
clusters = clusters,
ncluster = max(clusters)
))
}
bt.summary.report <- function(models, title, data, obj=NULL,
control = list(
plot.weight.transition.maps = F,
plot.risk.contribution.transition.maps = !is.null(obj)
)
) {
if(is.null(control$plot.weight.transition.maps)) control$plot.weight.transition.maps = F
if(is.null(control$plot.risk.contribution.transition.maps)) control$plot.risk.contribution.transition.maps = obj != NULL
filename = title
filename.pdf = paste(filename, '.pdf', sep='')
filename.csv = paste(filename, '.csv', sep='')
pdf(file = filename.pdf, width=8.5, height=11)
layout(1:2)
plotbt(models, plotX = T, log = 'y', LeftMargin = 3, main = title)
mtext('Cumulative Performance', side = 2, line = 1)
out = plotbt.strategy.sidebyside(models, return.table=T)
cdi = custom.composite.diversification.indicator(obj, plot.main = F, plot.table = F)
out = rbind(colMeans(cdi, na.rm=T), out)
rownames(out)[1] = 'Composite Diversification Indicator(CDI)'
y = 100 * sapply(models, compute.turnover, data)
out = rbind(y, out)
rownames(out)[1] = 'Portfolio Turnover'
performance.barchart.helper(out, 'Sharpe,Cagr,DVR,MaxDD,Volatility,Portfolio Turnover,Composite Diversification Indicator(CDI)', c(T,T,T,T,F,F,T))
if(control$plot.weight.transition.maps) {
layout(1:4)
for(m in names(models)) {
plotbt.transition.map(models[[m]]$weight, name=m)
legend('topright', legend = m, bty = 'n')
}
}
if(control$plot.risk.contribution.transition.maps) {
dates = index(data$prices)[obj$period.ends]
layout(1:4)
for(m in names(models)) {
plotbt.transition.map(make.xts(obj$risk.contributions[[m]], dates),
name=paste('Risk Contributions',m))
legend('topright', legend = m, bty = 'n')
}
}
dev.off()
load.packages('abind')
append=FALSE
cat(title, '\n', file=filename.csv, append=append)
write.table(out, sep=',',  row.names = , col.names = NA, file=filename.csv, append=TRUE)
cat('\n\n', file=filename.csv, append=TRUE)
if(F) {
out = abind(lapply(models, function(m) m$equity))
colnames(out) = names(models)
write.xts(make.xts(out, index(models[[1]]$equity)), filename.csv, append=TRUE)
}
return(out)
}
pie.labels.fix <- function (x, y, angles, labels, radius = 1, ...) {
par(xpd = TRUE)
xylim <- par("usr")
plotdim <- par("pin")
yradius <- radius * (xylim[4] - xylim[3])/(xylim[2] - xylim[1]) * plotdim[1]/plotdim[2]
xc <- cos(angles) * radius + x
yc <- sin(angles) * yradius + y
text(xc, yc, labels, ...)
par(xpd = FALSE)
}
plot.cluster.weights <- function(weight, clusters, main='') {
load.packages('RColorBrewer,plotrix')
clusters = sort(clusters)
weight = weight[names(clusters)]
weight.cluster = tapply(weight,clusters,sum)
counts = tapply(names(clusters),clusters,len)
ncluster = len(counts)
require(RColorBrewer)
colors = colorRampPalette(brewer.pal(iif(ncluster>9,9,ncluster),'Set1'))(ncluster)
cols = c()
for(i in 1:ncluster) cols = c(cols, col.add.alpha(colors[i], seq(200,100,length.out = counts[i])))
if(F) {
plot(-1:1 ,-1:1, main=main, type='n', yaxt = 'n', xaxt = 'n', xlab = '', ylab = '', axes = F)
bisect.angles = floating.pie(0,0,weight, col=cols, border='white', radius=0.9, cex=0.8)
pie.labels(0,0,bisect.angles,names(weight),radius=1,bg=0,border=F, srt=bisect.angles)
}
par(mar = c(2,2,2,2))
pie(weight, col=cols, border='white', radius=0.9, main=main)
require(plotrix)
bisect.angles = floating.pie(0,0,weight.cluster,radius=0.5,col=colors,border='white')
pie.labels.fix(0,0,bisect.angles,paste(round(100*weight.cluster,0),'%',sep=''),radius=0.2)
}
adaptive.shrinkage.paper.backtests <- function()
{
stats = list()
for(i in 0:8) {
temp = bt.run.data.set(i)
stats[[temp$title]] = temp$stats
}
save(stats, file='stats.Rdata')
if(F) {
load(file='stats.Rdata')
names = c("Portfolio Turnover", "Composite Diversification Indicator(CDI)", "Cagr", "Sharpe", "Volatility", "MaxDD")
custom.order=c(F,T,T,T,F,F)
names = spl('Portfolio Turnover,Sharpe,Volatility,Composite Diversification Indicator(CDI)')
custom.order = c(T,F,T,F)
out = c()
for(n in names) {
temp = sapply(stats, function(x) x[n,])
temp = apply(temp, 1, as.double)
out = rbind(out, rowMeans(apply(temp, 1, rank)))
}
rownames(out) = names
performance.barchart.helper(out, custom.order = !custom.order,
custom.col=c('MV.S_SA_A'='red', 'best.sharpe'='green'))
}
bt.summary.stats(stats)
}
adaptive.shrinkage.paper.backtests.leverage <- function()
{
temp = bt.run.data.set(0)
vol = as.double(temp$stats['Volatility',])
leverage = vol[1] / vol
temp1 = bt.run.data.set(0, leverage = leverage, title = paste(temp$title, 'leverage'))
}
cssa.average.shrinkage <- function(hist) {
n = ncol(hist)
correlation = cor(hist, use='complete.obs', method='pearson')
s0 = apply(hist, 2, sd, na.rm=T)
index = s0==0
if(sum(index) > 0) {
correlation[index,] = 0
correlation[,index] = 0
}
column.avg = (colSums(correlation) - 1) / (n - 1)
new.matrix = outer(column.avg, column.avg, '+') / 2
diag(new.matrix) = 1
new.matrix * (s0 %*% t(s0))
}
shrink.average.cssa <- function(s=NULL) { s=s; function(x, a) cov.shrink(x, cssa.average.shrinkage, s, 1)$sigma }
bt.load.data.set <- function(data.set = 1)
{
if(data.set == 0) {
data <- new.env()
getSymbols.TB(env = data, auto.assign = T, download = F)
bt.prep(data, align='remove.na', dates='1990::')
return(list(data = data, title = 'Futures Forex', tickers = data$symbolnames))
}
data.sets = list()
data.sets[['Selected ETFs']] = list(
tickers = spl('SPY,QQQ,EEM,IWM,EFA,TLT,IYR,GLD'),
align='keep.all', dates='2002:08::')
data.sets[['Selected ETFs 1']] = list(
tickers = spl('GLD,UUP,SPY,QQQ,IWM,EEM,EFA,IYR,USO,TLT'),
align='keep.all', dates='2003:04::')
data.sets[['Selected ETFs 2']] = list(
tickers = spl('VTI,IEV,EEM,EWJ,AGG,GSG,GLD,ICF'),
align='keep.all', dates='2003:10::')
data.sets[['Dow stock (Engle)']] = list(
tickers = spl('AA,AXP,BA,CAT,DD,DIS,GE,IBM,IP,JNJ,JPM,KO,MCD,MMM,MO,MRK,MSFT'),
align='keep.all', dates='1980::')
data.sets[['Nasdaq 100']] = list(
tickers = spl('ATVI,ADBE,ALTR,AMZN,AMGN,APOL,AAPL,AMAT,ADSK,ADP,BBBY,BIIB,BRCM,CHRW,CA,CELG,CEPH,CERN,CHKP,CTAS,CSCO,CTXS,CTSH,CMCSA,COST,DELL,XRAY,DISH,EBAY,EA,EXPD,ESRX,FAST,FISV,FLEX,FLIR,FWLT,GILD,HSIC,HOLX,INFY,INTC,INTU,JBHT,KLAC,LRCX,LIFE,LLTC,LOGI,MAT,MXIM,MCHP,MSFT,MYL,NTAP,NWSA,NVDA,ORLY,ORCL,PCAR,PDCO,PAYX,PCLN,QGEN,QCOM,BBRY,ROST,SNDK,SIAL,SPLS,SBUX,SRCL,SYMC,TEVA,URBN,VRSN,VRTX,VOD,XLNX,YHOO'),
align='keep.all', dates='1995::')
data.sets[['Selected ETFs 10 Major Markets']] = list(
tickers = spl('SPY,EFA,EWJ,EEM,IYR,RWX,IEF,TLT,DBC,GLD'),
align='keep.all', dates='2002:08::')
data.sets[['Sector ETFs']] = list(
tickers = spl('XLY,XLP,XLE,XLF,XLV,XLI,XLB,XLK,XLU'),
align='keep.all', dates='1998:12::')
data.sets[['MSCI Country']] = list(
tickers = spl('SPY,TUR,EIRL,THD,EWL,EWK,EWA,EWD,EWW,EWN,EWJ,EWQ,EWO,EWP,EWH,EWG,MES,EWS,EIDO,EWU,EPOL,EWM,EIS,AFK,EWI,EWT,EWC,FXI,EZA,EWY,VNM,ECH,RSX,EWZ')	,
align='keep.all', dates='2000:08::')
load.packages('quantmod')
title = names(data.sets)[data.set]
info = data.sets[[data.set]]
data <- new.env()
getSymbols(info$tickers, src = 'yahoo', from = '1900-01-01', env = data, auto.assign = T)
data.clean(data, min.ratio=3)
for(i in ls(data)) data[[i]] = adjustOHLC(data[[i]], use.Adjusted=T)
bt.prep(data, align = info$align, dates = info$dates)
prices = coredata(data$prices)
prices[is.na(prices)] = mlag(prices)[is.na(prices)]
prices[is.na(prices)] = mlag(prices)[is.na(prices)]
data$prices[] = prices
return(list(data = data, title = title, tickers = info$tickers))
}
bt.run.data.set <- function(i = 1, dates = '::',
leverage = 1,
title = ''
)
{
data.set = bt.load.data.set(i)
data = data.set$data
if(nchar(title) == 0) title = data.set$title
obj = portfolio.allocation.helper(data$prices,
periodicity = 'months', lookback.len = 60,
min.risk.fns = list(MV=min.var.portfolio),
shrinkage.fns = list(
S=sample.shrinkage,
SA=sample.anchored.shrinkage,
D=shrink.diag(1),
CC=shrink.const.cor(1),
SI=shrink.single.index(1),
A=shrink.average.cssa(1),
AVG=function(h,a) 1/4*(
shrink.diag(1)(h,a) + shrink.const.cor(1)(h,a) +
shrink.single.index(1)(h,a) + shrink.average.cssa(1)(h,a)
),
S_SA_A=function(h,a) 1/3*(shrink.average.cssa(1)(h,a)
+ sample.shrinkage(h,a) + sample.anchored.shrinkage(h,a)),
D_S=shrink.diag(),
CC_S=shrink.const.cor(),
SI_S=shrink.single.index(),
A_S=shrink.average.cssa(),
AVG_S=function(h,a) cov.shrink(h, 1/4 *(
shrink.diag(1)(h,a) + shrink.const.cor(1)(h,a) +
shrink.single.index(1)(h,a) + shrink.average.cssa(1)(h,a)
))$sigma,
S_SA_A_S=function(h,a) cov.shrink(h, 1/3 *(
shrink.average.cssa(1)(h,a) + sample.shrinkage(h,a) + sample.anchored.shrinkage(h,a)
))$sigma
),
adjust2positive.definite = F,
custom.stats.fn = 'portfolio.allocation.custom.stats'
)
obj = bt.shrinkage.best.sharpe('S,SA,A', 252, obj, data)
models = create.strategies(obj, data, dates = dates, leverage = leverage)$models
list(title = title, stats = bt.summary.report(models, title, data, obj,
control = list(plot.weight.transition.maps = F,
plot.risk.contribution.transition.maps = F)
)
)
}
bt.shrinkage.best.sharpe <- function(methods, sharpe.lookback.len, obj, data)
{
models = create.strategies(obj, data)$models
methods = spl(methods)
methods = paste('MV',methods,sep='.')
global.data <- new.env()
for(i in methods) {
temp = models[[i]]$equity
temp[1:min(which(temp != 1))] = NA
global.data[[i]] = make.stock.xts(temp)
}
bt.prep(global.data, align='keep.all')
rets = global.data$prices / mlag(global.data$prices) - 1
sharpe = bt.apply.matrix(rets, runMean, sharpe.lookback.len) / bt.apply.matrix(rets, runSD, sharpe.lookback.len)
sharpe = ifna(sharpe,-Inf)
index.best = unlist(apply(sharpe[obj$period.ends,],1,which.max))
obj$weights$best.sharpe = obj$weights[[1]]
for(m in 1:len(methods))
obj$weights$best.sharpe[ index.best == m, ] = obj$weights[[ methods[m] ]][ index.best == m, ]
n = 1/len(methods)
obj$weights$simple.avg = obj$weights[[1]] * 0
for(m in 1:len(methods))
obj$weights$simple.avg[] = obj$weights$simple.avg[] + n * obj$weights[[ methods[m] ]]
return(obj);
}
bt.summary.report <- function(models, title, data, obj=NULL,
control = list(plot.weight.transition.maps = F,
plot.risk.contribution.transition.maps = !is.null(obj)
)
) {
if(is.null(control$plot.weight.transition.maps)) control$plot.weight.transition.maps = F
if(is.null(control$plot.risk.contribution.transition.maps)) control$plot.risk.contribution.transition.maps = obj != NULL
filename = title
filename.pdf = paste(filename, '.pdf', sep='')
filename.csv = paste(filename, '.csv', sep='')
pdf(file = filename.pdf, width=8.5, height=11)
layout(1:2)
plotbt(models, plotX = T, log = 'y', LeftMargin = 3, main = title)
mtext('Cumulative Performance', side = 2, line = 1)
out = plotbt.strategy.sidebyside(models, return.table=T)
cdi = custom.composite.diversification.indicator(obj, plot.main = F, plot.table = F)
out = rbind(colMeans(cdi, na.rm=T), out)
rownames(out)[1] = 'Composite Diversification Indicator(CDI)'
y = 100 * sapply(models, compute.turnover, data)
out = rbind(y, out)
rownames(out)[1] = 'Portfolio Turnover'
performance.barchart.helper(out, 'Sharpe,Cagr,DVR,MaxDD,Volatility,Portfolio Turnover,Composite Diversification Indicator(CDI)', c(T,T,T,T,F,F,T))
if(control$plot.weight.transition.maps) {
layout(1:len(models))
for(m in names(models)) {
plotbt.transition.map(models[[m]]$weight, name=m)
legend('topright', legend = m, bty = 'n')
}
}
if(control$plot.risk.contribution.transition.maps) {
dates = index(data$prices)[obj$period.ends]
layout(1:len(models))
layout(1:4)
for(m in names(models)) {
plotbt.transition.map(make.xts(obj$risk.contributions[[m]], dates),
name=paste('Risk Contributions',m))
legend('topright', legend = m, bty = 'n')
}
}
dev.off()
load.packages('abind')
append=FALSE
cat(title, '\n', file=filename.csv, append=append)
write.table(out, sep=',',  row.names = , col.names = NA, file=filename.csv, append=TRUE)
cat('\n\n', file=filename.csv, append=TRUE)
if(F) {
out = abind(lapply(models, function(m) m$equity))
colnames(out) = names(models)
write.xts(make.xts(out, index(models[[1]]$equity)), filename.csv, append=TRUE)
}
return(out)
}
bt.summary.stats <- function(stats)
{
names = spl('Portfolio Turnover,Sharpe,Volatility,Composite Diversification Indicator(CDI)')
custom.order = c(T,F,T,F)
out = stats[[1]][names,]
for(i in 1:len(names)) {
temp = apply(sapply(stats, function(x) x[names[i],]), 1, function(x) mean(as.double(x)))
temp = pnorm(temp, mean(temp), sd(temp))
if(custom.order[i]) temp = 1 - temp
out[i,] = temp
}
score = apply(out[1:3,], 2, function(x) mean(as.double(x)))
all.out1 = rbind(out, score)
rownames(all.out1)[5] = 'Composite Standardized Score (higher is better)'
out = stats[[1]][names,]
for(i in 1:len(names)) {
temp = apply(sapply(stats, function(x) x[names[i],]), 2, function(x) rank(iif(custom.order[i],1,-1)*as.double(x)))
out[i,] = rowMeans(temp)
}
score = apply(out[1:3,], 2, function(x) mean(as.double(x)))
all.out2 = rbind(out, score)
rownames(all.out2)[5] = 'Composite Rank Score (lower is better)'
write.table(rbind(all.out1, '', all.out2), sep=',',  row.names = , col.names = NA, file='summary.csv')
out = read.csv('summary.csv', stringsAsFactors=F)
temp = as.matrix(out[,-1,drop=F])
colnames(temp) = colnames(out)[-1]
rownames(temp) = out[,1]
pdf(file = 'summary.pdf', width=8.5, height=11)
performance.barchart.helper(temp, 'Composite Standardized Score (higher is better),Composite Rank Score (lower is better)', c(T,F), nc.plot=1,
custom.col=c('MV.S_SA_A'='red', 'best.sharpe'='green'))
dev.off()
png(filename = 'plot2.png', width = 600, height = 400, units = 'px', pointsize = 12, bg = 'white')
performance.barchart.helper(temp, 'Composite Rank Score (lower is better)', c(F), nc.plot=1,
custom.col=c('MV.S_SA_A'='red', 'best.sharpe'='green'))
dev.off()
}
tableColor <- function(data, header.col='LightGray', row.col='yellow', negative.col='red', ...) {
if (is.null(data))
return("")
add.to.row= NULL
if( nrow(data) > 1) {
temp = as.list(seq(1,nrow(data)-1,by=2))
add.to.row=list(temp, rep("XXX", len(temp)))
}
temp = renderTable(data, add.to.row=add.to.row, ...)
temp = temp()
if(!is.na(negative.col) && !is.null(negative.col))
temp = gsub("(-\\d*\\.\\d*)", paste("<font color=", negative.col, ">\\1</font>"), temp)
if(!is.na(row.col) && !is.null(row.col))
temp = gsub("XXX<TR>", paste("<TR bgcolor=", row.col, ">"), temp)
temp = gsub("XXX", "", temp)
if(!is.na(header.col) && !is.null(header.col))
temp = gsub("<TR>\\s*<TH>\\s*</TH>",paste("<TR bgcolor=", header.col, "><TH></TH>"), temp)
return(temp)
}

createNonReactiveTextInput <- function(id, label, value, 
                                       button.label = '') {
  label = substitute(label)
  button.label = substitute(button.label)
  if(button.label != '')
    list(
      tagList(
        tags$label(label), tags$input(id = id, type = "text", value = value, 
                                      style = "display:none;"), 
        tags$input(id = paste(id, "Temp", sep=''), type = "text", 
                   value = value, style = "display:inline;", 
                   onkeypress = paste("{if (event.keyCode==13) $('#", id, 
                                      "TempChange').click()}", sep = ''))), 
      div(
        tags$button(id = paste(id, "TempChange", sep=''), type = "button", 
                    class = "btn btn-primary", onclick = paste("$('#", id, 
                      "').val($('#", id, "Temp').val()).change();", 
                      sep = ''), button.label)))
  else
    list(
      tagList(
        tags$label(label), tags$input(id = id, type = "text", value = value, 
                                      style = "display:none;"), 
        tags$input(id = paste(id, "Temp", sep = ''), type = "text", 
                   value = value, style = "display:inline;", 
                   onkeypress = paste("{if (event.keyCode==13) $('#", id, 
                                      "').val($('#", id, 
                                      "Temp').val()).change()}", 
                                      sep = ''))))
  }

createNonReactiveTextInputCustom <- function(id, label, tag.label = 'input', 
                                             button.label = '', 
                                             enableEnter = TRUE, opts) {
  label = substitute(label)
  button.label = substitute(button.label)
  
  onkeypress = ''
  if(button.label != '') {
    if(enableEnter)
      onkeypress = paste("{if (event.keyCode==13) $('#", id, 
                         "TempChange').click()}", sep = '')
    
    list(tagList(tags$label(label), 
                 tag(tag.label, c(id = id, style = "display:none;", opts)), 
                 tag(tag.label, c(id = paste(id, "Temp", sep = ''), 
                                  style = "display:inline;", 
                                  onkeypress = onkeypress, opts))),
         div(
           tags$button(id = paste(id, "TempChange", sep = ''), 
                       type = "button", class = "btn btn-primary", onclick = 
                         paste("$('#", id, "').val($('#", id, 
                               "Temp').val()).change();", sep = ''), 
                       button.label)))
    } else {
      if(enableEnter)
        onkeypress = paste("{if (event.keyCode==13) $('#", id, 
                           "').val($('#", id, 
                           "Temp').val()).change()}", sep = '')
      list(
        tagList(
          tags$label(label), tag(tag.label, c(id = id, style = 
                                                "display:none;", opts)), 
          tag(tag.label, c(id = paste(id,"Temp", sep = ''), 
                           style = "display:inline;", 
                           onkeypress = onkeypress, opts))))
    }
  }

strategy.load.historical.data <- function(tickers = spl('DIA,SPY'), 
                                          dates = '1900::', 
                                          align = 'keep.all', 
                                          fill.gaps = T, adjust = T, 
                                          current = F) {
  load.packages('quantmod')
  tickers = spl(tickers)
  data <- new.env()
  
  for(i in 1:len(tickers)) {
    ticker0 = gsub('\\^', '', tickers[i])
    temp = try(getSymbols(tickers[i], src = 'yahoo', from = '1900-01-01', 
                          env = data, auto.assign = T), TRUE)
    if(inherits(temp, 'try-error'))
      cat(i, 'out of', len(tickers), 'Error Reading', tickers[i], '\n', 
          sep='\t')
    else if(is.null(data[[ tickers[i] ]]))
      if( is.null(data[[ ticker0 ]]) )
        cat(i, 'out of', len(tickers), 'Error Reading', tickers[i], '\n', 
            sep = '\t')
    else
      cat(i, 'out of', len(tickers), 'Reading', ticker0, 
          format(range(index(data[[ ticker0 ]])), '%d-%b-%Y'), '\n', 
          sep = '\t')
    else
      cat(i, 'out of', len(tickers), 'Reading', tickers[i], 
          format(range(index(data[[ tickers[i] ]])), '%d-%b-%Y'), '\n', 
          sep = '\t')
  }
  if(adjust) for(i in ls(data)) data[[i]] = adjustOHLC(data[[i]], 
                                                       use.Adjusted = T) 
  if(current) {
    quotes = getQuote(tickers)
    for(i in ls(data))
      if( last(index(data[[i]])) < as.Date(quotes[i, 'Trade Time']) ) {
        data[[i]] = 
          rbind(data[[i]], 
                make.xts(quotes[i, spl('Open,High,Low,Last,Volume,Last')], 
                         as.Date(quotes[i, 'Trade Time'])))
      }
    }
  bt.prep(data, align = align, dates = dates, fill.gaps = fill.gaps)
  return(data)
  }

performance.barchart.helper <- function(out,
names = rownames(out),
custom.order=rep(T,len(spl(names))),
nplots.page = len(spl(names)),
nc.plot = 2,
sort.performance = T,
custom.col = NULL
)
{
layout(mat=matrix(1:(nplots.page + nplots.page %% nc.plot), nc=nc.plot, byrow=FALSE))
par(mar=c(4, 3, 2, 2))
col = spl('ivory2,red')
names = spl(names)
names(custom.order) = names
for(i in names) {
y = as.double(out[i,])
index = iif(sort.performance, order(y, decreasing = custom.order[i]), 1:len(y))
cols = iif(y[index] > 0, col[1], col[2])
names(cols) = colnames(out)[index]
if(!is.null(custom.col))
cols[names(custom.col)] = custom.col
x = barplot(y[index], names.arg = '',
col=cols,
main=i,
border = 'darkgray',las=2)
grid(NA,NULL)
abline(h=0, col='black')
if(y[1] > 0)
text(x, 0 * x, colnames(out)[index], adj=c(-0.1,1), srt=90, xpd = TRUE)
else
text(x, 0 * x, colnames(out)[index], adj=c(1.1,1), srt=90, xpd = TRUE)
if(sort.performance) {
mtext('worst', side = 1,line = 0, outer = F, adj = 1, font = 1, cex = 1)
mtext('best', side = 1,line = 0, outer = F, adj = 0, font = 1, cex = 1)
}
}
}
barplot.with.labels <- function(data, main, plotX = TRUE, label=c('level','name','both')) {
par(mar=c( iif(plotX, 6, 2), 4, 2, 2))
x = barplot(100 * data, main = main, las = 2, names.arg = iif(plotX, names(data), ''))
if(label[1] == 'level') text(x, 100 * data, round(100 * data,1), adj=c(0.5,1), xpd = TRUE)
if(label[1] == 'name') text(x, 0 * data, names(data), adj=c(-0.1,1), srt=90, xpd = TRUE)
if(label[1] == 'both')
text(x, 0 * data, paste(round(100 * data), '% ', names(data), sep=''), adj=c(-0.1,1), srt=90, xpd = TRUE)
}
strategy.performance.snapshoot <- function(models, one.page = F, title = NULL, data = NULL,
control = list(main = T, comparison = T, transition = T, monthly = T),
sort.performance = T
) {
for(n in spl('main,comparison,transition,monthly'))
if(is.null(control[[n]])) control[[n]] = F
out = NULL
if(control$main) {
layout(1:2)
plotbt(models, plotX = T, log = 'y', LeftMargin = 3, main = title)
mtext('Cumulative Performance', side = 2, line = 1)
out = plotbt.strategy.sidebyside(models, return.table=T)
}
if(one.page) return()
if(control$comparison) {
if(is.null(out))
out = plotbt.strategy.sidebyside(models, return.table=T, make.plot = F)
if(!is.null(data)) {
y = 100 * sapply(models, compute.turnover, data)
out = rbind(y, out)
rownames(out)[1] = 'Turnover'
performance.barchart.helper(out, 'Sharpe,Cagr,DVR,MaxDD,Volatility,Turnover', c(T,T,T,T,F,F), sort.performance = sort.performance)
} else
performance.barchart.helper(out, 'Sharpe,Cagr,DVR,MaxDD', c(T,T,T,T), sort.performance = sort.performance)
}
if(control$transition) {
layout(1:min(4,len(models)))
for(m in names(models)) {
plotbt.transition.map(models[[m]]$weight, name=m)
legend('topright', legend = m, bty = 'n')
}
}
if(control$monthly) {
layout(1:min(4,len(models)))
for(n in names(models))
plotbt.monthly.table(models[[n]]$equity, smain=n)
}
}
meom.strategy <- function
(
tickers = spl('DIA,SPY'),
dates = '1900::'
)
{
data = strategy.load.historical.data(tickers, dates, fill.gaps=T)
prices = data$prices
n = ncol(prices)
nperiods = nrow(prices)
month.ends = endpoints(prices, 'months')
month.ends = month.ends[month.ends > 0]
month.ends2 = iif(month.ends + 2 > nperiods, nperiods, month.ends + 2)
month.ends1 = iif(month.ends + 1 > nperiods, nperiods, month.ends + 1)
models = list()
data$weight[] = NA
data$weight[] = ntop(prices, n)
models$equal.weight = bt.run.share(data, clean.signal=F)
data$weight[] = NA
data$weight[month.ends,] = ntop(prices, n)[month.ends,]
data$weight[month.ends2,] = 0
models$meom.equal.weight = bt.run.share(data, clean.signal=F)
buy.rule = prices > bt.apply.matrix(prices, function(x) { WMA(x, 89) } )
buy.rule = ifna(buy.rule, F)
ret2 = ifna(prices / mlag(prices, 2), 0)
position.score = bt.apply.matrix(ret2, SMA, 5) * bt.apply.matrix(ret2, SMA, 40)
position.score[!buy.rule] = NA
data$weight[] = NA;
data$weight[month.ends,] = ntop(position.score[month.ends,], 2)
data$weight[month.ends2,] = 0
models$meom.top2.rank1 = bt.run.share(data, clean.signal=F, trade.summary=T)
position.score = bt.apply.matrix(ret2, SMA, 5) * mlag( bt.apply.matrix(ret2, SMA, 10), 5)
position.score[!buy.rule] = NA
data$weight[] = NA
data$weight[month.ends,] = ntop(position.score[month.ends,], 2)
data$weight[month.ends2,] = 0
models$meom.top2.rank2 = bt.run.share(data, clean.signal=F, trade.summary=T)
data$weight[] = NA
data$weight[month.ends,] = ntop(position.score[month.ends,], 2)
data$weight[month.ends2,] = 0
popen = bt.apply(data, Op)
data$weight[month.ends1,] = iif((prices > popen)[month.ends1,], 0, NA)
models$meom.top2.rank2.hold12 = bt.run.share(data, clean.signal=F, trade.summary=T)
return(rev(models))
}
meom.strategy.test <- function()
{
models = meom.strategy(
tickers = 'DIA,EEM,EFA,EWH,EWJ,EWT,EWZ,FXI,GLD,GSG,IEF,ILF,IWM,IYR,QQQ,SPY,VNQ,XLB,XLE,XLF,XLI,XLP,XLU,XLV,XLY,XLK',
dates='1995::'
)
png(filename = 'plot1.png', width = 600, height = 500, units = 'px', pointsize = 12, bg = 'white')
plotbt.custom.report.part1(models)
dev.off()
png(filename = 'plot2.png', width = 1200, height = 800, units = 'px', pointsize = 12, bg = 'white')
plotbt.custom.report.part2(models)
dev.off()
png(filename = 'plot3.png', width = 600, height = 500, units = 'px', pointsize = 12, bg = 'white')
plotbt.custom.report.part3(models$meom.top2.rank2, trade.summary=T)
dev.off()
}
timing.strategy <- function
(
tickers = spl('DIA,SPY,SHY'),
dates = '1900::',
periodicity = 'months',
ma.len = 200,
cash = 'SHY'
)
{
data = strategy.load.historical.data(tickers, dates)
prices = data$prices
n = ncol(prices)
nperiods = nrow(prices)
period.ends = endpoints(data$prices, periodicity)
period.ends = period.ends[period.ends > 0]
models = list()
position.score = prices
position.score$SHY = NA
weight = ntop(position.score[period.ends,], n)
data$weight[] = NA
data$weight[period.ends,] = weight
models$equal.weight = bt.run.share(data, clean.signal=F)
sma = bt.apply.matrix(prices, SMA, ma.len)
buy.rule = prices > sma
buy.rule = ifna(buy.rule, F)
weight = ntop(position.score[period.ends,], n)
weight[!buy.rule[period.ends,]] = 0
weight$SHY = 1 - rowSums(weight)
data$weight[] = NA
data$weight[period.ends,]  = weight
models$timing = bt.run.share(data, clean.signal=F, trade.summary=T)
return(rev(models))
}
timing.strategy.test <- function()
{
models = timing.strategy(
tickers = 'VTI,EFA,IEF,ICF,DBC,SHY',
dates='2002:08::'
)
plotbt.custom.report.part1(models)
plotbt.custom.report.part2(models)
}
rotation.strategy <- function
(
tickers = spl('DIA,SPY,SHY'),
dates = '1900::',
periodicity = 'months',
top.n = 2,
keep.n = 6
)
{
data = strategy.load.historical.data(tickers, dates)
prices = data$prices
n = ncol(prices)
nperiods = nrow(prices)
period.ends = endpoints(data$prices, periodicity)
period.ends = period.ends[period.ends > 0]
models = list()
data$weight[] = NA
data$weight[period.ends,] = ntop(prices[period.ends,], n)
models$equal.weight = bt.run.share(data, clean.signal=F)
position.score = prices / mlag(prices, 126)
data$weight[] = NA
data$weight[period.ends,] = ntop(position.score[period.ends,], top.n)
models$top = bt.run.share(data, clean.signal=T, trade.summary=T)
data$weight[] = NA
data$weight[period.ends,] = ntop.keep(position.score[period.ends,], top.n, keep.n)
models$top.keep = bt.run.share(data, clean.signal=T, trade.summary=T)
return(rev(models))
}
rotation.strategy.test <- function()
{
models = rotation.strategy(
tickers = 'XLY,XLP,XLE,XLF,XLV,XLI,XLB,XLK,XLU,IWB,IWD,IWF,IWM,IWN,IWO,IWP,IWR,IWS,IWV,IWW,IWZ',
dates='1970::'
)
plotbt.custom.report.part1(models)
plotbt.custom.report.part2(models)
}
create.ia <- function(hist.returns, index=1:ncol(hist.returns), nperiod=nrow(hist.returns))
{
ia = list()
ia$hist.returns = hist.returns
ia$nperiod = nperiod
ia$index = index
ia$n = ncol(ia$hist.returns)
ia$symbols = colnames(ia$hist.returns)
ia$risk = apply(ia$hist.returns, 2, sd, na.rm = T)
ia$correlation = cor(ia$hist.returns, use='complete.obs',method='pearson')
ia$cov = ia$correlation * (ia$risk %*% t(ia$risk))
ia$expected.return = apply(ia$hist.returns, 2, mean, na.rm = T)
return(ia)
}
update.ia <- function(ia, name, cov.shrink)
{
if(name != 'sample') {
ia$cov = cov.shrink
s0 = 1 / sqrt(diag(ia$cov))
ia$correlation = ia$cov * (s0 %*% t(s0))
}
ia
}
create.ia.subset <- function(ia.base, index=1:ia.base$n)
{
ia = list()
ia$hist.returns = ia.base$hist.returns[,index,drop=F]
ia$nperiod = ia.base$nperiod
ia$index = ia.base$index[index]
ia$n = ncol(ia$hist.returns)
ia$symbols = colnames(ia$hist.returns)
ia$risk = ia.base$risk[index]
ia$correlation = ia.base$correlation[index,index,drop=F]
ia$cov = ia.base$cov[index,index,drop=F]
ia$expected.return = ia.base$expected.return[index]
return(ia)
}
create.ia.period <- function
(
prices,
periodicity = 'weeks',
period.ends = endpoints(prices, periodicity)
)
{
prices = prices[period.ends,,drop=F]
ret = coredata(prices / mlag(prices) - 1)
function(hist.returns, index, nperiod)
{
i = nperiod
create.ia(ret[which(
period.ends <= i &
period.ends >= (i - nrow(hist.returns) + 1)
), index, drop=F],
index,
nperiod)
}
}
create.historical.ia <- function
(
hist.returns,
annual.factor
) {
ia = create.ia(hist.returns)
ia$annual.factor = annual.factor
ia$arithmetic.return = apply(hist.returns, 2, mean, na.rm = T)
ia$geometric.return = apply(hist.returns, 2, function(x) prod(1+x)^(1/len(x))-1 )
ia$arithmetic.return = (1 + ia$arithmetic.return)^ia$annual.factor - 1
ia$geometric.return = (1 + ia$geometric.return)^ia$annual.factor - 1
ia$risk = sqrt(ia$annual.factor) * ia$risk
ia$risk = iif(ia$risk == 0, 0.000001, ia$risk)
ia$cov = ia$correlation * (ia$risk %*% t(ia$risk))
ia$expected.return = ia$arithmetic.return
ia
}
ia.build.hist <- function(hist.returns, lookbacks, n.lag)
{
nperiods = nrow(hist.returns)
temp = c()
for (n.lookback in lookbacks)
temp = rbind(temp, hist.returns[(nperiods - n.lookback - n.lag + 1):(nperiods - n.lag), , drop=F])
return(temp)
}
momentum.averaged <- function(prices,
lookbacks = c(20,60,120,250) ,
n.lag = 3
) {
momentum = 0 * prices
for (n.lookback in lookbacks) {
part.mom = mlag(prices, n.lag) / mlag(prices, n.lookback + n.lag) - 1
momentum = momentum + 252 / n.lookback * part.mom
}
momentum / len(lookbacks)
}
create.ia.averaged <- function(lookbacks, n.lag)
{
lookbacks = lookbacks
n.lag = n.lag
function(hist.returns, index, nperiod)
{
temp = ia.build.hist(hist.returns, lookbacks, n.lag)
create.ia(temp, index, nperiod)
}
}
static.weight.portfolio <- function(static.allocation)
{
static.allocation = static.allocation
function
(
ia,
constraints
)
{
return(static.allocation[ia$index])
}
}
equal.weight.portfolio <- function
(
ia,
constraints
)
{
rep(1/ia$n, ia$n)
}
get.risky.asset.index <- function(ia) {
if(is.null(ia$risk)) ia$risk = sqrt(diag(ia$cov))
(ia$risk > 0) & !is.na(ia$risk)
}
set.risky.asset <- function(x, risk.index) {
out = rep(0, len(risk.index))
out[risk.index] = x
return(out)
}
risk.parity.portfolio.basic <- function
(
ia,
constraints
)
{
if(is.null(ia$risk)) ia$risk = sqrt(diag(ia$cov))
risk.index = get.risky.asset.index(ia)
x = 1 / ia$risk[risk.index]
set.risky.asset(x / sum(x), risk.index)
}
risk.parity.portfolio <- function(
risk.fn = function(ia) ia$risk
)
{
algo.map = list(
'cvar' = function(ia) -apply(ia$hist.returns, 2, compute.cvar),
'md' = function(ia) -apply(apply(1+ia$hist.returns, 2, cumprod), 2, compute.max.drawdown),
'cdar' = function(ia) -apply(apply(1+ia$hist.returns, 2, cumprod), 2, compute.cdar)
)
fn = try( match.fun(risk.fn) , silent = TRUE)
if(class(fn)[1] == 'try-error' && is.character(risk.fn) && any(names(algo.map) == tolower(risk.fn)))
fn = algo.map[[ tolower(risk.fn) ]]
if(class(fn)[1] == 'try-error') stop(paste('risk.parity.portfolio', fn))
function
(
ia,
constraints
)
{
if(is.null(ia$risk)) ia$risk = sqrt(diag(ia$cov))
risk.index = get.risky.asset.index(ia)
x = 1 / fn(ia)[risk.index]
x[x < 0] = 0
ifna(set.risky.asset(x / sum(x), risk.index), 0)
}
}
min.var.portfolio <- function
(
ia,
constraints,
cov.matrix = ia$cov,
dvec = rep(0, ia$n)
)
{
risk.index = get.risky.asset.index(ia)
Dmat = cov.matrix[risk.index, risk.index]
sol = try(solve.QP(Dmat=Dmat,
dvec=dvec[risk.index],
Amat=constraints$A[risk.index,,drop=F],
bvec=constraints$b,
meq=constraints$meq), silent = TRUE)
if(inherits(sol, 'try-error'))
sol = try(solve.QP(Dmat=make.positive.definite(Dmat, 0.000000001),
dvec=dvec[risk.index],
Amat=constraints$A[risk.index,,drop=F],
bvec=constraints$b,
meq=constraints$meq), silent = TRUE)
if(inherits(sol, 'try-error')) {
gia <<- ia
stop(sol)
}
set.risky.asset(sol$solution, risk.index)
}
max.div.portfolio <- function
(
ia,
constraints
)
{
risk.index = get.risky.asset.index(ia)
x = min.var.portfolio(ia, constraints, ia$correlation)
x = x[risk.index] / ia$risk[risk.index]
set.risky.asset(x / sum(x), risk.index)
}
equal.risk.contribution.portfolio <- function
(
ia,
constraints
)
{
risk.index = get.risky.asset.index(ia)
cov = ia$cov[risk.index, risk.index]
fn <- function(x){
if (sum(x) == 0) x = x + 1e-2
x  = x / sum(x)
risk.contribution = (x * (cov %*% x))
var(as.double(risk.contribution))
}
x0 = 1/sqrt(diag(cov))
x0 = x0 / sum(x0)
if(!is.null(constraints$x0))
if(all(!is.na(constraints$x0)))
if( sum(constraints$x0) == 1 )
if( fn(x0) > fn(constraints$x0[risk.index]) )
x0 = constraints$x0[risk.index]
load.packages('nloptr')
x = nloptr( x0=x0,eval_f=fn,lb = constraints$lb[risk.index],ub = constraints$ub[risk.index],
opts = list('algorithm'='NLOPT_LN_BOBYQA','xtol_rel'=1.0e-10))
set.risky.asset(x$solution / sum(x$solution), risk.index)
}
ef.portfolio <- function(percent = 0.5)
{
percent = as.double(percent[1])
if(percent > 1) percent = percent / 100
function
(
ia,
constraints
)
{
max.w = max.return.portfolio(ia, constraints)
min.w = min.var.portfolio(ia, constraints)
max.r = sum(max.w * ia$expected.return)
min.r = sum(min.w * ia$expected.return)
target = (max.r - min.r) * percent + min.r
constraints = add.constraints(ia$expected.return,
type='>=', b=target, constraints)
return(min.var.portfolio(ia, constraints))
}
}
min.te.portfolio.test <- function()
{
min.te.portfolio <- function(ia, constraints, index.weight)
min.var.portfolio(ia, constraints, dvec = index.weight %*% ia$cov)
load.packages('quantmod')
tickers = spl('SPY,TLT,XLP')
data <- new.env()
getSymbols.extra(tickers, src = 'yahoo', from = '1980-01-01', env = data, set.symbolnames = T, auto.assign = T)
for(i in data$symbolnames) data[[i]] = adjustOHLC(data[[i]], use.Adjusted=T)
bt.prep(data, align='keep.all', fill.gaps = T)
prices = data$prices
n = ncol(prices)
nperiods = nrow(prices)
ret = prices / mlag(prices) - 1
lookback.len = 120
hist = ret[(nperiods - lookback.len) : nperiods, , drop=F]
ia = create.ia(hist, nperiod=nperiods)
index.weight = coredata(last(prices))
index.weight = index.weight / sum(index.weight)
load.packages('quadprog,corpcor,kernlab')
constraints = create.basic.constraints(n, 0, 1, 1)
x = min.te.portfolio(ia, constraints, index.weight)
x - index.weight
252 * portfolio.return(x, ia)
constraints = create.basic.constraints(n, 0, 0.4, 1)
x = min.te.portfolio(ia, constraints, index.weight)
x - index.weight
252 * portfolio.return(x, ia)
constraints = create.basic.constraints(n, 0, 1, 1)
constraints = add.constraints(ia$expected.return, type = '>=', b=0.95 * max(ia$expected.return), constraints)
x = min.te.portfolio(ia, constraints, index.weight)
x - index.weight
252 * portfolio.return(x, ia)
}
random.hist <- function(
portfolio.fn = min.var.portfolio,
nsamples = 100,
sample.len = 60
)
{
fn = try( match.fun(portfolio.fn) , silent = TRUE)
if(class(fn)[1] == 'try-error') stop(paste('random.hist', fn))
nsamples = nsamples
sample.len = sample.len
function
(
ia,
constraints
)
{
load.packages('MASS')
weight = fn(ia, constraints)
ia.original = ia
for(i in 1:nsamples) {
ia$hist.returns = mvrnorm(sample.len, ia.original$expected.return, Sigma = ia.original$cov)
ia$expected.return = apply(ia$hist.returns, 2, mean)
ia$risk = apply(ia$hist.returns, 2, sd)
ia$correlation = cor(ia$hist.returns, use = 'complete.obs', method = 'pearson')
ia$cov = ia$correlation * (ia$risk %*% t(ia$risk))
weight = weight + fn(ia, constraints)
}
weight / (nsamples + 1)
}
}
random.hist.weight = function(
fn,
data,
period.ends,
nsamples = 100,
sample.len = 120,
silent = F
)
{
prices = data$prices * data$universe
obj = fn(prices, data)
scenarios = asset.paths.at.period.ends(data$prices, period.ends, nsamples, lookback.len=sample.len)
for(j in 1:nsamples) {
prices = make.xts(scenarios[,,j],data$dates)
colnames(prices) = colnames(data$prices)
prices = prices * data$universe
temp = fn(prices, data)
for(i in names(obj$weights))
obj$weights[[i]] = obj$weights[[i]] + temp$weights[[i]]
if (j %% 5 == 0)
if (!silent)
cat(j, 'done', '\n')
}
for(i in names(obj$weights))
obj$weights[[i]] = obj$weights[[i]] / (nsamples + 1)
obj
}
max.sharpe.portfolio.helper <- function
(
ia,
const = spl('long-only,long-short,market-neutral'),
const.sum = 1,
rf = 0
)
{
const = const[1]
n = ia$n
constraints = new.constraints(n)
if( const == 'long-only' )
constraints = add.constraints(diag(n), type='>=', b=0, constraints)
excess.return = ia$expected.return - rf
if( all(excess.return <= 0) )
constraints = add.constraints(excess.return, -1 , type = '=', constraints)
else
constraints = add.constraints(excess.return, 1 , type = '=', constraints)
if( const == 'market-neutral' )
constraints = add.constraints(rep(1,n), 0 , type = '=', constraints)
weight = min.var.portfolio(ia,constraints)
if( const == 'market-neutral' )
return(2*const.sum * weight / sum(abs(weight)))
else
return(const.sum * weight / sum(weight))
}
max.sharpe.portfolio <- function
(
const = spl('long-only,long-short,market-neutral'),
const.sum = 1
)
{
const = const[1]
const.sum = const.sum
function(ia, constraints) { max.sharpe.portfolio.helper(ia, const, const.sum) }
}
max.sharpe.portfolio.test <- function()
{
ia = aa.test.create.ia()
n = ia$n
constraints = create.basic.constraints(n, 0, 1, 1)
ef = portopt(ia, constraints, 50, 'Efficient Frontier')
png(filename = 'plot1.png', width = 500, height = 500, units = 'px', pointsize = 12, bg = 'white')
plot.ef(ia, list(ef), transition.map=F)
max(portfolio.return(ef$weight,ia) /  portfolio.risk(ef$weight,ia))
weight = min.var.portfolio(ia,constraints)
points(100 * portfolio.risk(weight,ia), 100 * portfolio.return(weight,ia), pch=15, col='red')
portfolio.return(weight,ia) /  portfolio.risk(weight,ia)
weight = max.sharpe.portfolio()(ia,constraints)
points(100 * portfolio.risk(weight,ia), 100 * portfolio.return(weight,ia), pch=15, col='orange')
portfolio.return(weight,ia) /  portfolio.risk(weight,ia)
plota.legend('Minimum Variance,Maximum Sharpe','red,orange', x='topright')
dev.off()
weight = max.sharpe.portfolio('long-only')(ia,constraints)
round(weight,2)
round(c(sum(weight[weight<0]), sum(weight[weight>0])),2)
weight = max.sharpe.portfolio('long-short')(ia,constraints)
round(weight,2)
round(c(sum(weight[weight<0]), sum(weight[weight>0])),2)
weight = max.sharpe.portfolio('long-short', -1)(ia,constraints)
round(weight,2)
round(c(sum(weight[weight<0]), sum(weight[weight>0])),2)
weight = max.sharpe.portfolio('market-neutral')(ia,constraints)
round(weight,2)
round(c(sum(weight[weight<0]), sum(weight[weight>0])),2)
}
max.div.portfolio.test <- function()
{
load.packages('quantmod')
tickers = spl('SPY,QQQ,EEM,IWM,EFA,TLT,IYR,GLD')
data = env()
getSymbols.extra(tickers, src = 'yahoo', from = '1980-01-01', env = data, set.symbolnames = T, auto.assign = T)
for(i in data$symbolnames) data[[i]] = adjustOHLC(data[[i]], use.Adjusted=T)
bt.prep(data, align='keep.all', fill.gaps = T)
prices = data$prices
period.ends = date.ends(prices, 'month')
hist.returns = prices[period.ends] / mlag(prices[period.ends]) - 1
ia = create.historical.ia(hist.returns, 12)
constraints = create.basic.constraints(ia$n, 0, 1, 1)
load.packages('quadprog,corpcor,lpSolve,kernlab')
sol1 = max.div.portfolio(ia, constraints)
max.div.portfolio2 <- function
(
ia,
constraints
)
{
risk.index = get.risky.asset.index(ia)
constraints = new.constraints(ia$n, lb = 0)
constraints = add.constraints(diag(ia$n), type = ">=", b = 0, constraints)
constraints = add.constraints(ia$risk, type = '=', b = 1, constraints)
x = min.var.portfolio(ia, constraints)
set.risky.asset(x / sum(x), risk.index)
}
sol2 = max.div.portfolio2(ia, constraints)
round(100*cbind(sol1,sol2),2)
round(100*t(portfolio.risk.contribution(rbind(sol1,sol2), ia)),2)
ia = list(
n = 2,
symbols = spl('A,B'),
risk = c(20, 10) / 100,
correlation = matrix(c(1, 0.5, 0.5, 1), 2, 2),
expected.return = c(1,1)
)
ia$cov = ia$correlation * (ia$risk %*% t(ia$risk))
ia.base = ia
constraints = create.basic.constraints(ia$n, 0, 1, 1)
sol = list(
EW = equal.weight.portfolio(ia, constraints),
ERC = equal.risk.contribution.portfolio(ia, constraints),
MV = min.var.portfolio(ia, constraints),
MDP = max.div.portfolio(ia, constraints),
MDP2 = max.div.portfolio1(ia, constraints)
)
round(100*t(sapply(sol,c)))
round(100*portfolio.risk.contribution(t(sapply(sol,c)), ia))
ia = list(
n = 3,
symbols = spl('A,A1,B'),
risk = c(20, 20, 10) / 100,
correlation = matrix(c(1, 1, 0.5, 1, 1, 0.5, 0.5, 0.5, 1), 3, 3),
expected.return = c(1,1,1)
)
ia$cov = ia$correlation * (ia$risk %*% t(ia$risk))
constraints = create.basic.constraints(ia$n, 0, 1, 1)
sol = list(
EW = equal.weight.portfolio(ia, constraints),
ERC = equal.risk.contribution.portfolio(ia, constraints),
MV = min.var.portfolio(ia, constraints),
MDP = max.div.portfolio(ia, constraints),
MDP2 = max.div.portfolio1(ia, constraints)
)
weights = t(sapply(sol,c))
round(100*weights)
weights.base = cbind(rowSums(weights[,1:2]), weights[,3])
round(100*weights.base)
round(100*portfolio.risk.contribution(weights.base, ia.base))
ia = list(
n = 2,
symbols = spl('LA,B'),
risk = c(5, 10) / 100,
correlation = matrix(c(1, 0.5, 0.5, 1), 2, 2),
expected.return = c(1,1)
)
ia$cov = ia$correlation * (ia$risk %*% t(ia$risk))
constraints = create.basic.constraints(ia$n, 0, 1, 1)
sol = list(
EW = equal.weight.portfolio(ia, constraints),
ERC = equal.risk.contribution.portfolio(ia, constraints),
MV = min.var.portfolio(ia, constraints),
MDP = max.div.portfolio(ia, constraints),
MDP2 = max.div.portfolio1(ia, constraints)
)
weights = t(sapply(sol,c))
round(100*weights)
weights.base = cbind(weights[,1] / 4, weights[,2])
weights.base = weights.base/rowSums(weights.base)
round(100*weights.base)
round(100*portfolio.risk.contribution(weights.base, ia.base))
ia = list(
n = 3,
symbols = spl('A,B,Polico'),
risk = c(20, 10, 11.46) / 100,
correlation = matrix(c(1, 0.5, 0.982, 0.5, 1, 0.655, 0.982, 0.655, 1), 3, 3),
expected.return = c(1,1,1)
)
ia$cov = ia$correlation * (ia$risk %*% t(ia$risk))
constraints = create.basic.constraints(ia$n, 0, 1, 1)
sol = list(
EW = equal.weight.portfolio(ia, constraints),
ERC = equal.risk.contribution.portfolio(ia, constraints),
MV = min.var.portfolio(ia, constraints),
MDP = max.div.portfolio(ia, constraints),
MDP2 = max.div.portfolio1(ia, constraints)
)
weights = t(sapply(sol,c))
round(100*weights)
weights.base = cbind(weights[,1] +  1/2 * weights[,3], weights[,2] +  1/4 * weights[,3])
weights.base = weights.base/rowSums(weights.base)
round(100*weights.base)
round(100*portfolio.risk.contribution(weights.base, ia.base))
}
max.sharpe.nlp.portfolio <- function
(
ia,
constraints
)
{
risk.index = get.risky.asset.index(ia)
cov = ia$cov[risk.index, risk.index]
er = ia$expected.return[risk.index]
fn <- function(x){
(-x %*% er / sqrt( x %*% cov %*% x ))[1]
}
x0 = constraints$ub[risk.index] / sum(constraints$ub[risk.index])
load.packages('Rsolnp')
x = solnp(x0, fn, eqfun = function(w) sum(w), eqB   = 1,
LB = constraints$lb[risk.index], UB = constraints$ub[risk.index],
control = list(trace=0))
set.risky.asset(x$pars, risk.index)
}
find.portfolio.given.risk.test <- function()
{
ia = aa.test.create.ia()
n = ia$n
constraints = create.basic.constraints(n, 0, 1, 1)
load.packages('quadprog,corpcor,lpSolve,kernlab')
weight = max.return.portfolio(ia, constraints)
weight
cat(round(100*c(portfolio.return(weight,ia), portfolio.risk(weight,ia), portfolio.return(weight,ia) /  portfolio.risk(weight,ia)),2), '\n')
weight = min.var.portfolio(ia,constraints)
weight
cat(round(100*c(portfolio.return(weight,ia), portfolio.risk(weight,ia), portfolio.return(weight,ia) /  portfolio.risk(weight,ia)),2), '\n')
weight = target.return.portfolio.helper(ia,constraints, 12/100)
weight
cat(round(100*c(portfolio.return(weight,ia), portfolio.risk(weight,ia), portfolio.return(weight,ia) /  portfolio.risk(weight,ia)),2), '\n')
weight = target.risk.portfolio.helper(ia,constraints, 10/100, silent=F)
weight
cat(round(100*c(portfolio.return(weight,ia), portfolio.risk(weight,ia), portfolio.return(weight,ia) /  portfolio.risk(weight,ia)),2), '\n')
}
min.corr.excel <- function(power.function = 1, final.scale = 'risk')
{
power.function = as.numeric(power.function)
final.scale = final.scale
function
(
ia,
constraints
)
{
risk.index = get.risky.asset.index(ia)
n = sum(risk.index)
x = min.corr.special.case(ia$risk[risk.index])
if(!is.null(x)) return( set.risky.asset(x / sum(x), risk.index) )
cor.matrix = ia$correlation[risk.index, risk.index]
upper.index = upper.tri(cor.matrix)
cor.m = cor.matrix[upper.index]
cor.mu = mean(cor.m)
cor.sd = sd(cor.m)
avg.corr.contribution = (rowSums(cor.matrix) - 1) / (n - 1)
avg.rank = rank(avg.corr.contribution)
avg.rank = avg.rank ^ power.function
rr.adjustment = avg.rank / sum(avg.rank)
norm.dist.m = 0 * cor.matrix
diag(norm.dist.m) = 0
norm.dist.m[upper.index] = 1-pnorm(cor.m, cor.mu, cor.sd)
norm.dist.m = (norm.dist.m + t(norm.dist.m))
rr.norm.dist.m = rep.col(rr.adjustment,n) * norm.dist.m
rr.norm.dist = colSums(rr.norm.dist.m)
rr.weighted = rr.norm.dist / sum(rr.norm.dist)
vol.scale = ia$risk[risk.index]
if(final.scale == 'vol') vol.scale = diag(ia$cov[risk.index, risk.index])
inverse.volatility.weight = (1 / vol.scale) / sum(1 / vol.scale)
x = rr.weighted * inverse.volatility.weight / sum(rr.weighted * inverse.volatility.weight)
set.risky.asset(x / sum(x), risk.index)
}
}
min.corr.excel.portfolio.test <- function()
{
ia = list()
ia$n = 3
ia$risk = c(14, 18, 22) / 100;
ia$correlation = matrix(
c(1, 0.90, 0.85,
0.90, 1, 0.70,
0.85, 0.70, 1), nr=3, byrow=T)
ia$cov = ia$correlation * (ia$risk %*% t(ia$risk))
constraints = create.basic.constraints(ia$n, 0, 1, 1)
min.corr.excel.portfolio(ia,constraints)
}
min.corr.special.case <- function(risk) {
n = len(risk)
if(n == 1) 1
else if(n == 2) 1/risk
else NULL
}
min.corr <- function(power.function = 1)
{
power.function = as.numeric(power.function)
function
(
ia,
constraints
)
{
risk.index = get.risky.asset.index(ia)
x = min.corr.special.case(ia$risk[risk.index])
if(!is.null(x)) return( set.risky.asset(x / sum(x), risk.index) )
cor.matrix = ia$correlation[risk.index, risk.index]
upper.index = upper.tri(cor.matrix)
cor.m = cor.matrix[upper.index]
cor.mu = mean(cor.m)
cor.sd = sd(cor.m)
norm.dist.m = 0 * cor.matrix
diag(norm.dist.m) = NA
norm.dist.m[upper.index] = 1-pnorm(cor.m, cor.mu, cor.sd)
norm.dist.m = (norm.dist.m + t(norm.dist.m))
norm.dist.avg = rowMeans(norm.dist.m, na.rm=T)
norm.dist.rank = rank(-norm.dist.avg)
norm.dist.rank = norm.dist.rank ^ power.function
norm.dist.weight = norm.dist.rank / sum(norm.dist.rank)
diag(norm.dist.m) = 0
weighted.norm.dist.average = norm.dist.weight %*% norm.dist.m
final.weight = weighted.norm.dist.average / sum(weighted.norm.dist.average)
x = final.weight / ia$risk[risk.index]
set.risky.asset(x / sum(x), risk.index)
}
}
min.corr2 <- function(power.function = 1)
{
power.function = as.numeric(power.function)
function
(
ia,
constraints
)
{
risk.index = get.risky.asset.index(ia)
x = min.corr.special.case(ia$risk[risk.index])
if(!is.null(x)) return( set.risky.asset(x / sum(x), risk.index) )
cor.matrix = ia$correlation[risk.index, risk.index]
cor.m = cor.matrix
diag(cor.m) = 0
avg = rowMeans(cor.m)
cor.mu = mean(avg)
cor.sd = sd(avg)
norm.dist.avg = 1-pnorm(avg, cor.mu, cor.sd)
norm.dist.rank = rank(-norm.dist.avg)
norm.dist.rank = norm.dist.rank ^ power.function
norm.dist.weight = norm.dist.rank / sum(norm.dist.rank)
weighted.norm.dist.average = norm.dist.weight %*% (1-cor.m)
final.weight = weighted.norm.dist.average / sum(weighted.norm.dist.average)
x = final.weight / ia$risk[risk.index]
set.risky.asset(x / sum(x), risk.index)
}
}
min.var2111 <- function(power.function = 1)
{
power.function = as.numeric(power.function)
function
(
ia,
constraints
)
{
risk.index = get.risky.asset.index(ia)
x = min.corr.special.case(ia$risk[risk.index])
if(!is.null(x)) return( set.risky.asset(x / sum(x), risk.index) )
data.matrix = ia$cov[risk.index, risk.index]
avg = rowMeans(data.matrix)
data.mu = mean(avg)
data.sd = sd(avg)
norm.dist.avg = 1 - pnorm(avg, data.mu, data.sd)
norm.dist.rank = rank(-norm.dist.avg)
norm.dist.rank = norm.dist.rank ^ power.function
norm.dist.weight = norm.dist.rank / sum(norm.dist.rank)
weighted.norm.dist.average = norm.dist.weight %*% (max(data.matrix) - data.matrix)
final.weight = weighted.norm.dist.average / sum(weighted.norm.dist.average)
x = final.weight / ia$risk[risk.index]
set.risky.asset(x / sum(x), risk.index)
}
}
min.var.excel <- function(power.function = 1)
{
power.function = as.numeric(power.function)
function
(
ia,
constraints
)
{
risk.index = get.risky.asset.index(ia)
x = min.corr.special.case(ia$risk[risk.index])
if(!is.null(x)) return( set.risky.asset(x / sum(x), risk.index) )
data.matrix = ia$cov[risk.index, risk.index]
avg = rowMeans(data.matrix)
data.mu = mean(avg)
data.sd = sd(avg)
norm.dist.avg = 1 - pnorm(avg, data.mu, data.sd)
final.weight = norm.dist.avg / sum(norm.dist.avg)
x = final.weight / diag(data.matrix)
set.risky.asset(x / sum(x), risk.index)
}
}
min.var.excel.portfolio.test <- function()
{
tickers = spl('DBC,EEM,EWJ,GLD')
data = '
-0.004903678  0.005815362  0.006696429  -0.010055275
0.000703977  0.01035895  0.014412417  0.006355806
0.000351741  0.007868383  0.005464481  0.000708299
-0.000351617  -0.002838893  0.006521739  -0.004423735
-0.015124868  -0.015421115  -0.012958963  -0.010782629
-0.004642857  0.009638554  0.014223195  0.003653351
'
ia = create.ia(matrix(scan(text=data), nc = len(tickers)))
data = '
0.000090  0.000044  0.000028  0.000034
0.000044  0.000084  0.000068  0.000039
0.000028  0.000068  0.000101  0.000036
0.000034  0.000039  0.000036  0.000039
'
ia$cov = matrix(scan(text=data), nc = ia$n)
ia$risk = sqrt(diag(ia$cov))
constraints = create.basic.constraints(ia$n, 0, 1, 1)
min.var.excel.portfolio(ia,constraints)
}
min.var2.portfolio <- function(ia,constraints) { min.corr.excel(final.scale = 'vol')(ia,constraints) }
min.var2 <- function(power.function = 1)
{
power.function = as.numeric(power.function)
function(ia,constraints) min.corr.excel(power.function, final.scale = 'vol')(ia,constraints)
}
min.var.excel.portfolio <- function(ia,constraints) { min.var.excel()(ia,constraints) }
min.corr.excel.portfolio <- function(ia,constraints) { min.corr.excel()(ia,constraints) }
min.corr.portfolio <- function(ia,constraints) { min.corr()(ia,constraints) }
min.corr2.portfolio <- function(ia,constraints) { min.corr2()(ia,constraints) }
min.cvar <- function(alpha = 0.95) {
alpha = alpha
function(ia,constraints) {
ia$parameters.alpha = as.numeric(alpha)
min.cvar.portfolio(ia,constraints)
}
}
min.cdar <- function(alpha = 0.95) {
alpha = alpha
function(ia,constraints) {
ia$parameters.alpha = as.numeric(alpha)
min.cdar.portfolio(ia,constraints)
}
}
min.risk.downside <- function(mar = 0) {
mar = mar
function(ia,constraints) {
ia$parameters.mar = as.numeric(mar)
min.risk.downside.portfolio(ia,constraints)
}
}
min.mad.downside <- function(mar = 0) {
mar = mar
function(ia,constraints) {
ia$parameters.mar = as.numeric(mar)
min.mad.downside.portfolio(ia,constraints)
}
}
sample.shrinkage <- function( hist, hist.all ) {
cov(hist, use='complete.obs', method='pearson')
}
sample.anchored.shrinkage <- function( hist, hist.all ) {
cov(hist.all, use='complete.obs', method='pearson')
}
sample.mix.shrinkage <- function( hist, hist.all ) {
0.5 * sample.shrinkage(hist, hist.all) +
0.5 * sample.anchored.shrinkage(hist, hist.all)
}
exp.sample.shrinkage <- function( hist, hist.all ) {
hist = na.omit(hist)
lam = 0.9
i = 0:(nrow(hist)-1)
wt = lam^i
wt = wt/sum(wt)
cov.wt(hist, wt=rev(wt))$cov
}
diagonal.shrinkage <- function( hist, hist.all ) {
n = ncol(hist)
s0 = apply(hist, 2, sd, na.rm=T)
diag(n) * (s0 %*% t(s0))
}
average.shrinkage <- function( hist, hist.all ) {
n = ncol(hist)
correlation = cor(hist, use='complete.obs', method='pearson')
avg.correlation = (sum(correlation) - n) / (n*n - n)
create.cov.matrix(avg.correlation, hist)
}
min.shrinkage <- function( hist, hist.all ) {
correlation = cor(hist, use='complete.obs', method='pearson')
create.cov.matrix(min(correlation), hist)
}
max.shrinkage <- function( hist, hist.all ) {
n = ncol(hist)
correlation = cor(hist, use='complete.obs', method='pearson')
create.cov.matrix(max(correlation-diag(n)), hist)
}
create.cov.matrix <- function( value, hist ) {
n = ncol(hist)
s0 = apply(hist, 2, sd, na.rm=T)
temp = diag(n)
((matrix(1,n,n) - temp) * value + temp) * (s0 %*% t(s0))
}
ledoit.wolf.shrinkage <- function( hist, hist.all ) {
require(BurStFin)
var.shrink.eqcor(hist, 1, compatible = T)
}
factor.model.shrinkage <- function( hist, hist.all ) {
require(BurStFin)
factor.model.stat(hist, 1)
}
cov.shrink <- function(h, prior = NULL, shrinkage = NULL, roff.method = 1) {
require(tawny)
T = nrow(h)
S = cov.sample(h)
if( is.function(prior) ) prior = prior(h)
if( is.null(prior) ) prior = tawny::cov.prior.cc(S)
if( is.null(shrinkage) ) {
p = tawny::shrinkage.p(h, S)
if( roff.method == 0 )
r = sum(p$diags, na.rm=TRUE)
else
r = tawny::shrinkage.r(h, S, p)
c = tawny::shrinkage.c(prior, S)
k = (p$sum - r) / c
shrinkage = max(0, min(k/T, 1))
}
return(list(sigma = shrinkage * prior + (1 - shrinkage) * S, shrinkage = shrinkage))
}
cov.sample <- function(h) {
T = nrow(h)
x = h - rep.row(colMeans(h), T)
(t(x) %*% x) / T
}
cov.const.cor <- function(h) {
sample = cov.sample(h)
n = ncol(h)
var = diag(sample)
cov.mat = sqrt(var) %*% t(sqrt(var))
avg.rho = (sum(sample / iif(cov.mat==0,1,cov.mat))-n)/(n*(n-1))
prior = avg.rho * cov.mat
diag(prior) = var
prior
}
cov.diag <- function(h) {
S = cov.sample(h)
diag(diag(S))
}
cov.market <- function(h) {
T = nrow(h)
x = h - rep.row(colMeans(h), T)
xmkt = rowMeans(x)
n = ncol(h)
sample=cov(cbind(x, xmkt))*(T - 1)/T
covmkt=sample[1:n,n+1]
varmkt=sample[n+1,n+1]
prior=covmkt %*% t(covmkt) / varmkt
diag(prior) = diag(sample[1:n,1:n])
prior
}
cov.2param <- function(h) {
sample = cov.sample(h)
n = ncol(h)
meanvar=mean(diag(sample))
meancov=(sum(sample) -  sum(diag(sample)))/(n*(n-1))
meanvar*diag(n) + meancov*(1-diag(n))
}
shrink.diag <- function(s=NULL) { s=s; function(x, a) { cov.shrink(x, cov.diag, s, 0)$sigma }}
shrink.const.cor <- function(s=NULL) { s=s; function(x, a) { cov.shrink(x, cov.const.cor, s, 1)$sigma }}
shrink.single.index <- function(s=NULL) { s=s; function(x, a) { cov.shrink(x, cov.market, s, 1)$sigma }}
shrink.two.parameter <- function(s=NULL) { s=s; function(x, a) { cov.shrink(x, cov.2param, s, 1)$sigma }}
empty.group <- function
(
ia
)
{
return( rep(1,ia$n) )
}
static.group <- function(group)
{
group = as.numeric(group)
function
(
ia
)
{
group[ia$index]
}
}
cluster.group.hclust <- function
(
ia
)
{
if(ia$n <= 2) return(c(1,1)[1:ia$n])
dissimilarity = 1 - ia$correlation
distance = as.dist(dissimilarity)
fit =  hclust(distance, method='ward')
minh = min(fit$height)
maxh = max(fit$height)
group = cutree(fit, h=minh + (maxh - minh) /3)
return( group )
}
cluster.group.kmeans.90 <- function
(
ia
)
{
if(ia$n <= 2) return(c(1,1)[1:ia$n])
dissimilarity = 1 - cor(ia$hist.returns, use='complete.obs',method='spearman')
distance = as.dist(dissimilarity)
n = ncol(ia$correlation)
n = ceiling(n*2/3)
xy = cmdscale(distance)
for (i in 2:n) {
fit = kmeans(xy, centers=i, iter.max=100, nstart=100)
p.exp = 1- fit$tot.withinss / fit$totss
if(p.exp > 0.9) break
}
group = fit$cluster
return( group )
}
cluster.group.kmeans.elbow <- function
(
ia
)
{
if(ia$n <= 2) return(c(1,1)[1:ia$n])
dissimilarity = 1 - cor(ia$hist.returns, use='complete.obs',method='spearman')
distance = as.dist(dissimilarity)
n = ncol(ia$correlation)
n = ceiling(n*2/3)
xy = cmdscale(distance)
p.exp = rep(NA, n)
for (i in 2:n) {
fit = kmeans(xy, centers=i, iter.max=100, nstart=100)
p.exp[i] = 1- fit$tot.withinss / fit$totss
}
icluster = find.maximum.distance.point(p.exp[-1]) + 1
fit = kmeans(xy, centers=icluster, iter.max=100, nstart=100)
group = fit$cluster
return( group )
}
cluster.group.FTCA <- function
(
threshold = 0.5
)
{
function
(
ia
)
{
n = ia$n
map.index = 1:n
min.cluster.group = 1
if (threshold >= 1) return(map.index)
group = rep(0, n)
names(group) = names(ia$risk)
index = rep(TRUE, n)
names(index) = names(ia$risk)
while (n > 0) {
if (n == 1) {
group[index] = min.cluster.group
break
} else {
cor.matrix = ia$correlation[index, index]
if (n == 2) {
if (cor.matrix[1,2] > threshold)
group[index] = min.cluster.group
else
group[index] = c(min.cluster.group, min.cluster.group + 1)
break
} else {
avg.corr.contribution = (rowSums(cor.matrix) - 1) / (n - 1)
avg.rank = rank(avg.corr.contribution)
tip = which.min(avg.rank)
top = which.max(avg.rank)
if (cor.matrix[tip,top] > threshold) {
group[index] = min.cluster.group
break
} else {
index.top = map.index[index][cor.matrix[,top] > threshold]
index.tip = map.index[index][cor.matrix[,tip] > threshold]
group[index.tip] = min.cluster.group
group[index.top] = min.cluster.group + 1
index[index.tip] = F
index[index.top] = F
min.cluster.group = min.cluster.group + 2
n = sum(index)
}
}
}
}
return(group)
}
}
cluster.group.FTCA.test <- function() {
load.packages('quantmod')
tickers = spl('GLD,UUP,SPY,QQQ,IWM,EEM,EFA,IYR,USO,TLT')
tickers = spl('GLD,TLT,SPY,IWM,QQQ,EFA,EEM,IYR')
tickers = spl('XLY,XLP,XLE,XLF,XLV,XLI,XLB,XLK,XLU')
data <- new.env()
getSymbols(tickers, src = 'yahoo', from = '1900-01-01', env = data, auto.assign = T)
for(i in ls(data)) data[[i]] = adjustOHLC(data[[i]], use.Adjusted=T)
bt.prep(data, align='keep.all')
portfolio.allocation.custom.stats.clusters <- function(x,ia) {
return(list(
clusters.FTCA = cluster.group.FTCA(0.5)(ia)
))
}
periodicity = 'months'
lookback.len = 252
obj = portfolio.allocation.helper(data$prices,
periodicity = periodicity, lookback.len = lookback.len,
min.risk.fns = list(EW=equal.weight.portfolio),
custom.stats.fn = portfolio.allocation.custom.stats.clusters
)
clusters = obj$clusters.FTCA$EW
clusters['2012:05::']
temp1 = clusters['2011::']
plot.data = coredata(temp1)
rownames(plot.data) = format(index.xts(temp1), '%Y%m')
plot.table(plot.data, highlight = plot.data + 1)
obj = portfolio.allocation.helper(data$prices,
periodicity = periodicity, lookback.len = lookback.len,
min.risk.fns = list(
C.EW.kmeans = distribute.weights(equal.weight.portfolio, cluster.group.kmeans.90),
C.EW.FTCA = distribute.weights(equal.weight.portfolio, cluster.group.FTCA(0.5)),
C.RP.kmeans = distribute.weights(risk.parity.portfolio(), cluster.group.kmeans.90),
C.RP.FTCA = distribute.weights(risk.parity.portfolio(), cluster.group.FTCA(0.5)),
C.MD.kmeans = distribute.weights(max.div.portfolio, cluster.group.kmeans.90),
C.MD.FTCA = distribute.weights(max.div.portfolio, cluster.group.FTCA(0.5)),
C.MV.kmeans = distribute.weights(min.var.portfolio, cluster.group.kmeans.90),
C.MV.FTCA = distribute.weights(min.var.portfolio, cluster.group.FTCA(0.5)),
C.MVE.kmeans = distribute.weights(min.var.excel.portfolio, cluster.group.kmeans.90),
C.MVE.FTCA = distribute.weights(min.var.excel.portfolio, cluster.group.FTCA(0.5)),
C.MCE.kmeans = distribute.weights(min.corr.excel.portfolio, cluster.group.kmeans.90),
C.MCE.FTCA = distribute.weights(min.corr.excel.portfolio, cluster.group.FTCA(0.5)),
C.MS.kmeans = distribute.weights(max.sharpe.portfolio(), cluster.group.kmeans.90),
C.MS.FTCA = distribute.weights(max.sharpe.portfolio(), cluster.group.FTCA(0.5)),
C.ERC.kmeans = distribute.weights(equal.risk.contribution.portfolio, cluster.group.kmeans.90),
C.ERC.FTCA = distribute.weights(equal.risk.contribution.portfolio, cluster.group.FTCA(0.5))
)
)
models = create.strategies(obj, data)$models
png(filename = 'plot1.png', width = 600, height = 500, units = 'px', pointsize = 12, bg = 'white')
strategy.performance.snapshoot(models, T)
dev.off()
png(filename = 'plot2.png', width = 600, height = 500, units = 'px', pointsize = 12, bg = 'white')
barplot.with.labels(sapply(models, compute.turnover, data), 'Average Annual Portfolio Turnover')
dev.off()
}
distribute.weights <- function
(
fn,
group.fn = NA,
fn.within = NA
)
{
fn = match.fun(fn)
if(!is.function(group.fn)) if(!is.na(group.fn)) group.fn = match.fun(group.fn)
if(!is.function(fn.within)) if(!is.na(fn.within)) fn.within = match.fun(fn.within)
if(!is.function(fn.within)) fn.within = fn
function
(
ia,
constraints
)
{
if(!is.function(group.fn)) return(fn(ia, constraints))
group = as.numeric(group.fn(ia))
groups = unique(group[!is.na(group)])
ngroups = len(groups)
if(ngroups == 1) return(fn.within(ia, constraints))
weight0 = rep(0, ia$n)
if(ngroups == 0) return(weight0)
group[is.na(group)] = Inf
hist.g = NA * ia$hist.returns[,1:ngroups]
for(g in 1:ngroups) {
index = group == groups[g]
if( sum(index) == 1 ) {
weight0[index] = 1
hist.g[,g] = ia$hist.returns[, index, drop=F]
} else {
ia.temp = create.ia.subset(ia, index)
constraints.temp = create.basic.constraints(ia.temp$n, 0, 1, 1)
w0 = match.fun(fn.within)(ia.temp, constraints.temp)
weight0[index] = w0
hist.g[,g] = ia.temp$hist.returns %*% w0
}
}
ia.g = create.ia(hist.g)
constraints.g = create.basic.constraints(ngroups, 0, 1, 1)
group.weights = match.fun(fn)(ia.g, constraints.g)
for(g in 1:ngroups)
weight0[group == groups[g]] = weight0[group == groups[g]] * group.weights[g]
weight0
}
}
get.algo <- function(algo.name, has.param = F) {
algo.map = list(
'cluster' = distribute.weights,
'max.sharpe' = max.sharpe.portfolio,
'risk.parity' = risk.parity.portfolio
)
if(any(names(algo.map) == algo.name))
if(has.param)
algo.map[[ algo.name ]]
else
algo.map[[ algo.name ]]()
else {
if(has.param) {
fn = try( match.fun(algo.name) , silent = TRUE)
if(class(fn)[1] == 'try-error')	fn = try( match.fun(paste0(algo.name, '.portfolio')) , silent = TRUE)
} else {
fn = try( match.fun(paste0(algo.name, '.portfolio')) , silent = TRUE)
if(class(fn)[1] == 'try-error')	fn = try( match.fun(algo.name) , silent = TRUE)
}
if(class(fn)[1] == 'try-error') stop(paste('get.algo', fn))
fn
}
}
get.group <- function(group.name) {
group.map = list(
'none' = empty.group,
'hclust' = cluster.group.hclust,
'kmeans90' = cluster.group.kmeans.90
)
if(any(names(group.map) == group.name))
group.map[[ group.name ]]
else
stop(paste('Unknown group', group.name))
}
map.min.risk.fns <- function(strategys) {
strategys = spl(strategys,';')
min.risk.fns = list()
strategys = strategys[ nchar(strategys) > 0]
for(i in 1:len(strategys)) {
temp = spl(strategys[i])
f.name = temp[1]
if(nchar(f.name) == 0 || f.name == 'Empty') f.name = paste(temp[-1], collapse='-')
temp = tolower(temp)[-1]
if(len(temp) == 1)
min.risk.fns[[ f.name ]] = get.algo(temp[1])
else {
if(temp[1] == 'cluster') {
params = trim(spl(temp[2], ':'))
min.risk.fns[[ f.name ]] = get.algo(temp[1], T)( get.algo(params[2]), get.group(params[1]) )
} else
min.risk.fns[[ f.name ]] = get.algo(temp[1],T)(temp[-1])
}
}
min.risk.fns
}
rso.portfolio <- function
(
weight.fn,
k,
s,
const.lb = 0,
const.ub = 1,
const.sum = 1
)
{
weight.fn = match.fun(weight.fn)
k = k
s = s
const.lb = const.lb
const.ub = const.ub
const.sum = const.sum
constraints0 = create.basic.constraints(k, const.lb, const.ub, const.sum)
function
(
ia,
constraints
)
{
constraints1 = constraints0
k1 = k
if(k > ia$n) {
k1 = ia$n
constraints1 = create.basic.constraints(k1, const.lb, const.ub, const.sum)
}
space = seq(1:ia$n)
index.samples =t(replicate(s, sample(space, size=k1)))
weight = matrix(NA, nrow = s, ncol = ia$n)
for(i in 1:s){
ia.temp = create.ia.subset(ia, index.samples[i,])
weight[i,index.samples[i,]] = weight.fn(ia.temp, constraints1)
}
final.weight = ifna(colMeans(weight, na.rm=T), 0)
final.weight / sum(final.weight)
}
}
portfolio.allocation.helper.parallel <- function
(
cores = 1,
prices,
periodicity = 'weeks',
period.ends = endpoints(prices, periodicity),
lookback.len = 60,
n.skip = 1,
universe = prices[period.ends,,drop=F]>0,
prefix = '',
min.risk.fns = 'min.var.portfolio',
custom.stats.fn = NULL,
shrinkage.fns = 'sample.shrinkage',
create.ia.fn = create.ia,
update.ia.fn = update.ia,
adjust2positive.definite = T,
silent = F,
log = log.fn(),
log.frequency = 10,
const.lb = 0,
const.ub = 1,
const.sum = 1
)
{
cores = round(cores)
if(cores <= 1)
return(portfolio.allocation.helper(prices, periodicity, period.ends, lookback.len, n.skip,
universe, prefix,
min.risk.fns, custom.stats.fn, shrinkage.fns,
create.ia.fn, update.ia.fn,
adjust2positive.definite, silent, log,
const.lb, const.ub, const.sum))
load.packages('foreach,doParallel')
cl<-makeCluster(cores)
registerDoParallel(cl, cores = cores)
period.ends = period.ends[period.ends > 0]
start.i = which(period.ends >= (lookback.len + n.skip))[1]
chunks = c(1, floor(seq(start.i, len(period.ends)+1, length.out = cores + 1)[-1]))
temp = 1:len(period.ends)
if( nrow(universe) != len(period.ends) ) {
if( nrow(universe) == nrow(prices) )
universe = universe[period.ends,,drop=F]
else
stop("universe incorrect number of rows")
}
universe[is.na(universe)] = F
out <- foreach(i=1:cores, .packages=spl('quantmod,SIT')) %dopar% {
new.universe = universe
new.universe[temp[-c(chunks[i] : (chunks[i+1]-1))],]=F
portfolio.allocation.helper(prices, periodicity, period.ends, lookback.len, n.skip,
universe = new.universe, prefix,
min.risk.fns, custom.stats.fn, shrinkage.fns,
create.ia.fn, update.ia.fn,
adjust2positive.definite, silent, log, log.frequency,
const.lb, const.ub, const.sum)
}
stopCluster(cl)
base.out = out[[1]]
for(i in 2:cores) {
include.index = temp[c(chunks[i] : (chunks[i+1]-1))]
for(v in names(out[[i]])) {
if(is.list(out[[i]][[v]]))
for(n in names(out[[i]][[v]]))
base.out[[v]][[n]][include.index,] = out[[i]][[v]][[n]][include.index,]
if(is.matrix(out[[i]][[v]]))
base.out[[v]][include.index,] = out[[i]][[v]][include.index,]
}
}
base.out
}
portfolio.allocation.helper <- function
(
prices,
periodicity = 'weeks',
period.ends = endpoints(prices, periodicity),
lookback.len = 60,
n.skip = 1,
universe = prices[period.ends,,drop=F]>0,
prefix = '',
min.risk.fns = 'min.var.portfolio',
custom.stats.fn = NULL,
shrinkage.fns = 'sample.shrinkage',
create.ia.fn = create.ia,
update.ia.fn = update.ia,
adjust2positive.definite = T,
silent = F,
log = log.fn(),
log.frequency = 10,
const.lb = 0,
const.ub = 1,
const.sum = 1
)
{
load.packages('quadprog,corpcor')
load.packages('quadprog,corpcor,lpSolve,kernlab')
period.ends = period.ends[period.ends > 0]
if( nrow(universe) != len(period.ends) ) {
if( nrow(universe) == nrow(prices) )
universe = universe[period.ends,,drop=F]
else
stop("universe incorrect number of rows")
}
universe[is.na(universe)] = F
if(len(const.lb) == 1) const.lb = rep(const.lb, ncol(prices))
if(len(const.ub) == 1) const.ub = rep(const.ub, ncol(prices))
if(is.character(min.risk.fns)) {
min.risk.fns = spl(min.risk.fns)
names(min.risk.fns) = min.risk.fns
min.risk.fns = as.list(min.risk.fns)
}
for(i in 1:len(min.risk.fns)) {
f = spl(names(min.risk.fns)[i], '_')
f.name = paste(prefix, gsub('\\.portfolio', '', f[1]),sep='')
if(is.character(min.risk.fns[[i]])) {
if(len(f) == 1) {
min.risk.fns[[ i ]] = match.fun(f[1])
} else {
f.name = paste(f.name, f[-1], sep='_')
min.risk.fns[[ i ]] = match.fun(f[1])(f[-1])
}
}
names(min.risk.fns)[i] = f.name
}
if(is.character(shrinkage.fns)) {
shrinkage.fns = spl(shrinkage.fns)
names(shrinkage.fns) = shrinkage.fns
shrinkage.fns = as.list(shrinkage.fns)
}
for(i in 1:len(shrinkage.fns)) {
f = names(shrinkage.fns)[i]
f.name = gsub('\\.shrinkage', '', f[1])
if(is.character(shrinkage.fns[[ i ]]))
shrinkage.fns[[ i ]] = match.fun(f)
names(shrinkage.fns)[i] = f.name
}
dates = index(prices)[period.ends]
weight = NA * prices[period.ends,,drop=F]
prices = coredata(prices)
ret = prices / mlag(prices) - 1
start.i = which(period.ends >= (lookback.len + n.skip))[1]
weight[] = 0
weights = list()
for(f in names(min.risk.fns))
for(c in names(shrinkage.fns))
weights[[ paste(f,c,sep='.') ]] = weight
custom = list()
if( !is.null(custom.stats.fn) ) {
custom.stats.fn = match.fun(custom.stats.fn)
dummy = matrix(NA, nr=nrow(weight), nc=len(weights))
colnames(dummy) = names(weights)
dummy = make.xts(dummy, dates)
temp = ret
temp[] = rnorm(prod(dim(ret)))
temp = custom.stats.fn(1:ncol(ret), create.ia(temp))
for(ci in names(temp)) {
temp1 = NA * dummy
if(len(temp[[ ci ]]) > 1) {
temp1 = list()
for(w in names(weights))
temp1[[w]] = NA * weights[[w]]
}
custom[[ ci ]] = temp1
}
}
index.map = 1:ncol(ret)
if(!is.na(start.i)) {
for( j in start.i:len(period.ends) ) {
i = period.ends[j]
hist = ret[ (i- lookback.len +1):i,, drop=F]
include.index = count(hist)== lookback.len
index = universe[j,] & include.index
n = sum(index)
if(n > 0) {
hist = hist[ , index, drop=F]
hist.all = ret[ 1:i, index, drop=F]
if(n > 1) {
constraints = create.basic.constraints(n, const.lb[index], const.ub[index], const.sum)
ia.base = create.ia.fn(hist, index.map[index], i)
for(c in names(shrinkage.fns)) {
cov.shrink = shrinkage.fns[[c]](hist, hist.all)
ia = update.ia.fn(ia.base, c, cov.shrink)
if(adjust2positive.definite) {
temp = try(make.positive.definite(ia$cov, 0.000000001), TRUE)
if(!inherits(temp, 'try-error')) ia$cov = temp
temp = try(make.positive.definite(ia$correlation, 0.000000001), TRUE)
if(!inherits(temp, 'try-error')) ia$correlation = temp
}
for(f in names(min.risk.fns)) {
fname = paste(f,c,sep='.')
if (j > 1) constraints$x0 = as.vector( weights[[ fname ]][(j-1),index] )
weights[[ fname ]][j,index] = min.risk.fns[[f]](ia, constraints)
}
}
} else {
ia = create.ia.fn(hist, index.map[index], i)
for(c in names(shrinkage.fns)) {
for(f in names(min.risk.fns)) {
fname = paste(f,c,sep='.')
weights[[ fname ]][j,index] = 1
}
}
}
if( !is.null(custom.stats.fn) ) {
for(w in names(weights)) {
x = as.vector(weights[[ w ]][j, index])
temp = custom.stats.fn(x, ia)
for(ci in names(temp)) {
if(is.list(custom[[ ci ]]))
custom[[ ci ]][[ w ]][j, index] = temp[[ ci ]]
else
custom[[ ci ]][j, w] = temp[[ ci ]]
}
}
}
}
if( j %% log.frequency == 0) if(!silent) log(j, percent = (j - start.i) / (len(period.ends) - start.i))
}
}
if( len(shrinkage.fns) == 1 ) {
names(weights) = gsub( paste('\\.', names(shrinkage.fns), '$', sep=''), '', names(weights) )
for(ci in names(custom))
names(custom[[ ci ]]) = gsub( paste('\\.', names(shrinkage.fns), '$', sep=''), '', names(custom[[ ci ]]) )
}
return(c(list(weights = weights, period.ends = period.ends,
periodicity = periodicity, lookback.len = lookback.len), custom))
}
portfolio.allocation.custom.stats <- function(x,ia) {
risk.contributions = portfolio.risk.contribution(x, ia)
return(list(
risk.contributions = risk.contributions,
degree.diversification = 1 - sqrt(x %*% ia$cov %*% x) / (ia$risk %*% x),
risk.gini = 1 - portfolio.concentration.gini.coefficient(risk.contributions)
))
}
portfolio.allocation.helper.basic <- function
(
prices,
periodicity = 'weeks',
period.ends = endpoints(prices, periodicity),
lookback.len = 60,
prefix = '',
universe = prices[period.ends,]>0,
min.risk.fns = 'min.var.portfolio',
silent = F
)
{
load.packages('quadprog,corpcor')
period.ends = period.ends[period.ends > 0]
universe[is.na(universe)] = F
if(is.character(min.risk.fns)) min.risk.fns = spl(min.risk.fns)
dates = index(prices)[period.ends]
prices = coredata(prices)
ret = prices / mlag(prices) - 1
start.i = which(period.ends >= (lookback.len + 1))[1]
weight = NA * prices[period.ends,]
weight[] = 0
weights = list()
for(f in min.risk.fns)
weights[[f]] = weight
for( j in start.i:len(period.ends) ) {
i = period.ends[j]
hist = ret[ (i- lookback.len +1):i, ]
include.index = count(hist)== lookback.len
index = universe[j,] & include.index
n = sum(index)
if(n > 0) {
if(n > 1) {
hist = hist[ , index]
constraints = create.basic.constraints(n, 0, 1, 1)
ia = list()
ia$index = index
ia$n = n
ia$hist.returns = hist
ia$expected.return = apply(hist, 2, mean)
ia$risk = apply(hist, 2, sd)
ia$correlation = cor(hist, use='complete.obs', method='pearson')
ia$cov = ia$correlation * (ia$risk %*% t(ia$risk))
for(f in min.risk.fns)
weights[[ f ]][j,index] = match.fun(f)(ia, constraints)
} else {
for(f in min.risk.fns)
weights[[ f ]][j,index] = 1
}
}
if( j %% 10 == 0) if(!silent) cat(j, '\n')
}
return(list(weights = weights, period.ends = period.ends,
periodicity = periodicity, lookback.len = lookback.len))
}
create.strategies <- function
(
obj,
data,
leverage = 1,
min.weight = NA,
round.weight = NA,
execution.price = NA,
close.all.positions.index = NULL,
silent = F,
log = log.fn(),
prefix = '',
suffix = '',
clean.signal = F,
...
)
{
if(len(leverage) > 1 || leverage[1] != 1) {
if(len(leverage) == 1) leverage = rep(leverage, len(obj$weights))
for(i in 1:len(obj$weights)) obj$weights[[i]] = leverage[i] * obj$weights[[i]]
}
if(!is.na(min.weight) && min.weight != 0) for(i in names(obj$weights))
obj$weights[[i]][] = bt.apply.min.weight(coredata(obj$weights[[i]]), min.weight)
if(!is.na(round.weight) && round.weight != 0) for(i in names(obj$weights)) {
obj$weights[[i]][] = bt.apply.round.weight(coredata(obj$weights[[i]]), round.weight)
obj$weights[[i]][] = bt.apply.round.weight(coredata(obj$weights[[i]]), round.weight)
obj$weights[[i]][] = bt.apply.round.weight(coredata(obj$weights[[i]]), round.weight)
}
models = list()
n = len(names(obj$weights))
for(j in 1:n) {
i = names(obj$weights)[j]
i = paste(prefix, i, suffix, sep='')
if(!silent) log(i, percent = j / n)
data$weight[] = NA
data$execution.price[] = execution.price
data$weight[obj$period.ends,] = obj$weights[[j]]
if( !is.null(close.all.positions.index) ) data$weight[close.all.positions.index,] = 0
models[[i]] = bt.run.share(data, clean.signal = clean.signal, silent = silent, ...)
models[[i]]$period.weight = obj$weights[[j]]
}
obj$models = models
return(obj)
}
asset.allocation.strategy.test <- function()
{
tickers = 'SPY,QQQ,EEM,IWM,EFA,TLT,IYR,GLD'
dates='2000::'
data = strategy.load.historical.data(tickers, dates)
obj = portfolio.allocation.helper(data$prices,
periodicity = 'months', lookback.len = 60,
min.risk.fns = list(
EW=equal.weight.portfolio,
RP=risk.parity.portfolio(),
MV=min.var.portfolio
),
custom.stats.fn = 'portfolio.allocation.custom.stats'
)
models = create.strategies(obj, data)$models
plotbt.custom.report.part1(models)
plotbt.custom.report.part2(models)
strategy.performance.snapshoot(models)
}
target.vol.strategy <- function(model, weight,
target = 10/100,
lookback.len = 21,
max.portfolio.leverage = 100/100,
annual.periods = 252
)
{
ret = diff(log(model$equity))
hist.vol.model = sqrt(annual.periods) * runSD(ret, n = lookback.len)
hist.vol.model = as.vector(hist.vol.model)
weight.target = weight * (target / hist.vol.model)
rs = rowSums(abs(weight.target), na.rm=T)
weight.target = weight.target / iif(rs > max.portfolio.leverage, rs/max.portfolio.leverage, 1)
return(weight.target)
}
calendar.signal <- function(key.date, ...) {
offsets = list( ... )
if( is.list(offsets[[1]]) ) offsets = offsets[[1]]
else {
default.names = as.character(substitute(c(...))[-1])
default.names = paste0('L(', default.names, ')')
if(is.null(names(offsets))) names(offsets) = default.names
else names(offsets) = iif(nchar(names(offsets))==0, default.names, names(offsets))
}
signals = list()
for(n in names(offsets)) {
offset = offsets[[n]]
signal = mlag(key.date, offset[1])
for(i in offset) signal = signal | mlag(key.date, i)
signals[[n]] = signal
}
signals
}
calendar.strategy <- function(data, ..., universe = data$prices > 0, do.lag.universe = 1, commission = 0) {
signals = list( ... )
if( is.list(signals[[1]]) ) signals = signals[[1]]
else {
default.names = as.character(substitute(c(...))[-1])
if(is.null(names(signals))) names(signals) = default.names
else names(signals) = iif(nchar(names(signals))==0, default.names, names(signals))
}
universe = mlag(universe, do.lag.universe)
models = list()
nassets = ncol(data$prices)
for(n in names(signals)) {
data$weight[] = NA
temp = ifna(universe & signals[[n]], F)
if(nassets == 1)
data$weight[] = temp
else
data$weight[] = ifna(temp / rowSums(temp),0)
models[[n]] = bt.run.share(data, do.lag = 0, trade.summary=T, clean.signal=T, commission = commission, silent=T)
}
models
}
last.trades <- function(..., n=20, make.plot=T, return.table=F, smain = NULL) {
models = variable.number.arguments( ... )
model = models[[1]]
name=ifnull(names(models),NULL)[1]
if(!is.null(model$trade.summary)) {
ntrades = min(n, nrow(model$trade.summary$trades))
trades = last(model$trade.summary$trades, ntrades)
if(!is.null(smain) || !is.null(name)) colnames(trades)[1] = iif(is.null(smain),name,smain)
if(make.plot) {
layout(1)
plot.table(trades)
}
if(return.table) trades
}
}
last.signals <- function(..., n=20, make.plot=T, return.table=F, smain = NULL) {
models = variable.number.arguments( ... )
model = models[[1]]
name=ifnull(names(models),NULL)[1]
if(!is.null(model$period.weight)) {
data = round(100*model$period.weight,0)
ntrades = min(n, nrow(data))
trades = last(data, ntrades)
if(!is.null(smain) || !is.null(name))
smain = iif(is.null(smain),name,smain)
else
smain = 'Date'
if(make.plot) {
layout(1)
plot.table(as.matrix(trades))
}
if(return.table) trades
}
}
ifna.prev <- function(y)
{
y1 = !is.na(y)
y1[1]=T
return( y[cummax( (1:length(y)) * y1 )]	)
}
ifna.prevx <- function(y) {
y1 = !is.na(y)
y1[1]=T
return( cummax( (1:length(y)) * y1 ) )
}
ifna.prevx.rev <- function(y) {
y1 = !is.na(y)
y1[length(y)] = T
y1[!y1] = Inf
rev(cummin(rev((1:length(y)) * y1)))
}
ifna.prevx.test <- function() {
y = c(NA,1,1,NA,2,2,NA,NA)
y[ifna.prevx(y)]
y[ifna.prevx.rev(y)]
}
cross <- function( array1, array2, eq=F ) {
iif(eq, array1 >= array2, array1 > array2) & iif(len(array1) > 1, mlag(array1), array1) < iif(len(array2) > 1, mlag(array2), array2)
}
cross.up <- function( array1, array2 ) { cross( array1, array2 ) }
cross.dn <- function( array1, array2 ) { cross( array2, array1 ) }
cross.up.eq <- function( array1, array2 ) { cross( array1, array2, T ) }
cross.dn.eq <- function( array1, array2 ) { cross( array2, array1, T ) }
percent.rank <- function
(
data,
n=252
)
{
pctRank <- function(x,i) sum(x[i,1] >= x[(i- (n-1) ):i,])
out = data
data = coredata(data)
if( is.null(dim(data)) ) dim(data) = c(len(data),1)
rng = n:len(data)
out[] = c( rep(NA,(n-1)), sapply(rng, function(i) pctRank(data, i) / n) )
return(out)
}
percent.rankM <- function
(
...,
n = 252
)
{
data = variable.number.arguments( ... )
out = data[[1]]
for(j in 1:len(data)) data[[j]] = coredata(data[[j]])
rank.data = data[[ len(data) ]]
pctRank <- function(x,i) sum(rank.data[i] >= x[(i- (n-1) ):i])
rng = n:len(rank.data)
out[] = 0
for(j in 1:len(data))
out[] = out[] + c( rep(NA,(n-1)), sapply(rng, function(i) pctRank(data[[j]], i) / n) )
return(out/len(data))
}
DV <- function
(
HLC,
n=2,
bounded=FALSE
)
{
hlMean = rowMeans( HLC[,-3] )
res = runMean( HLC[,3] / hlMean, n ) - 1
if(bounded) res = percent.rank(res, 252)
return(res)
}
DVI <- function
(
x,
n=250
)
{
ColumnC = ( x / runMean(x,3) ) - 1
ColumnD = ( runMean( ColumnC , 5 ) + ( runMean( ColumnC , 100 ) / 10 ) ) / 2
ColumnE = runMean( ColumnD , 5 )
ColumnF = iif( x > mlag(x) , 1 , -1 )
ColumnG = ( runSum( ColumnF , 10 ) + ( runSum( ColumnF , 100 ) / 10 ) ) / 2
ColumnH = runMean( ColumnG , 2 )
DVI.magnitude = percent.rank( ColumnE , n )
DVI.stretch = percent.rank( ColumnH, n )
DVI = ( 0.8 * DVI.magnitude ) + ( 0.2 * DVI.stretch )
return(list(DVI=DVI, DVI.magnitude=DVI.magnitude, DVI.stretch=DVI.stretch))
}
TSI <- function
(
HLC,
n=10
)
{
HLC = apply(HLC, 2, ifna.prev)
ratio = ( HLC[,3] - mlag(HLC[,3], n) ) / ATR( HLC , n )[, "atr"]
out = SMA( SMA( ratio , n ), 100 )
return(out)
}
ulcer.index <- function
(
x,
n=14
)
{
sqrt(runSum(( x / runMax(x,n) -1 )^2, n) / n)
}
ev.ratio <- function
(
data,
n = 252
)
{
ret = coredata(data / mlag(data) - 1)
rng = n:len(data)
out = data
out[] = c( rep(NA,(n-1)), sapply(rng,
function(i) {
r = ret[(i- (n-1) ):i]
-sum(r > 0) / n * sum(r[r > 0]) / sum(r[r < 0])
}))
return(out)
}
ntop <- function
(
data,
topn = 1,
dirMaxMin = TRUE
)
{
temp = coredata(data)
if(is.logical(temp)) temp[] = iif(!temp,NA,temp)
if(topn == ncol(data)) {
index = is.na(temp)
temp[index] = 0
temp[!index] = 1
out = data
out[] = ifna(temp / rowSums(temp),0)
return( out )
}
index.n = rowSums(!is.na(temp))
for( i in 1:nrow(data) ) {
if( index.n[i] > 0 ) {
o = sort.list(temp[i,], na.last = TRUE, decreasing = dirMaxMin)
temp[i,] = 0
n = min(topn, index.n[i])
temp[i,o[1:n]] = 1/n
} else temp[i,] = 0
}
out = data
out[] = temp
out
}
ntop.helper <- function
(
x,
n=1,
dirMaxMin = TRUE
)
{
x = as.vector(x)
index.n = sum(!is.na(x))
if( index.n > 0 ) {
o = sort.list(x, na.last=TRUE, decreasing = dirMaxMin)
x[] = 0
n = min(n, index.n)
x[o[1:n]] = 1/n
} else x[] = 0
x
}
ntop.speed.test <- function()
{
load.packages('quantmod')
tickers = spl('XLY,XLP,XLE,XLF,XLV,XLI,XLB,XLK,XLU,IWB,IWD,IWF,IWM,IWN,IWO,IWP,IWR,IWS,IWV,IWW,IWZ')
data <- new.env()
getSymbols(tickers, src = 'yahoo', from = '1970-01-01', env = data, auto.assign = T)
for(i in ls(data)) data[[i]] = adjustOHLC(data[[i]], use.Adjusted=T)
bt.prep(data, align='keep.all', dates='1970::2011')
prices = data$prices
n = len(tickers)
a = coredata(prices)
b = a
c = a
tic(12)
for( i in 1:nrow(b) ) {
b[i,] = ntop.helper(b[i,], n, T)
}
toc(12)
tic(12)
d = prices
for( i in 1:nrow(c) ) {
d[i,] = ntop.helper(c[i,], n, T)
}
toc(12)
range(b-d)
}
ntop.keep <- function
(
data,
topn = 1,
keepn = 1,
dirMaxMin = TRUE
)
{
temp = coredata(data)
if(is.logical(temp)) temp[] = iif(!temp,NA,temp)
index.n = rowSums(!is.na(temp))
for( i in 1:nrow(temp) ) {
if( index.n[i] > 0 ) {
x = temp[i,]
o = sort.list(x, na.last = TRUE, decreasing = dirMaxMin)
x[] = 0
n = min(topn, index.n[i])
x[o[1:n]] = 1
if( i >= 2 ) {
y = temp[(i-1),]
n1 = min(keepn, index.n[i])
y[-o[1:n1]] = 0
index1 = y != 0
index1.n = sum(index1)
if( index1.n > 0 ) {
x[] = 0
x[index1] = 1
for( j in 1:n ) {
if( sum(x) == topn ) break
x[o[j]] = 1
}
}
}
temp[i,] = x/sum(x)
} else temp[i,] = 0
}
out = data
out[] = temp
out
}
br.rank <- function(x)
{
t(apply(coredata(-x), 1, rank, na.last='keep'))
}
bt.rank <- function(x, dirMaxMin = TRUE, do.sort = F )
{
res = NA * x
res[] = {
if(!do.sort) {
if(dirMaxMin)
t(apply(coredata(-x), 1, rank, na.last='keep'))
else
t(apply(coredata(x), 1, rank, na.last='keep'))
} else {
index = 1:ncol(x)
x = coredata(x)
out = t(apply(x, 1, function(y) {
temp = sort.list(y, na.last = TRUE, decreasing = dirMaxMin)
temp[temp] = index
temp
}
))
out[is.na(x)] = NA
out
}
}
res
}
super.smoother.filter <- function(x) {
a1 = exp( -1.414 * pi / 10.0 )
b1 = 2.0 * a1 * cos( (1.414*180.0/10.0) * pi / 180.0 )
c2 = b1
c3 = -a1 * a1
c1 = 1.0 - c2 - c3
x = c1 * (x + mlag(x)) / 2
x[1] = x[2]
out = x * NA
out[] = filter(x, c(c2, c3), method='recursive', init=c(0,0))
out
}
roofing.filter <- function(x) {
alpha1 = (cos((0.707*360 / 48) * pi / 180.0 ) + sin((0.707*360 / 48) * pi / 180.0 ) - 1) / cos((0.707*360 / 48) * pi / 180.0 )
x = (1 - alpha1 / 2)*(1 - alpha1 / 2)*( x - 2*mlag(x) + mlag(x,2))
x[1] = x[2] = x[3]
HP = x * NA
HP[] = filter(x, c(2*(1 - alpha1), - (1 - alpha1)*(1 - alpha1)), method='recursive', init=c(0,0))
super.smoother.filter(HP)
}
roofing.stochastic.indicator  <- function(x, lookback = 20) {
Filt = roofing.filter(x)
HighestC = runMax(Filt, lookback)
HighestC[1:lookback] = as.double(HighestC[lookback])
LowestC = runMin(Filt, lookback)
LowestC[1:lookback] = as.double(LowestC[lookback])
Stoc = (Filt - LowestC) / (HighestC - LowestC)
super.smoother.filter(Stoc)
}
runQuantile = function(x, k, probs) {
load.packages('caTools')
temp = rep.col(x * NA, len(probs))
temp[k:len(x),] = runquantile(as.vector(coredata(x)), k, probs, endrule='trim')
temp
}
spl <- function
(
s,
delim = ','
)
{
unlist(strsplit(s,delim))
}
join <- function
(
v,
delim = ''
)
{
paste(v,collapse=delim)
}
trim <- function
(
s
)
{
s = sub(pattern = '^\\s+', replacement = '', x = s)
sub(pattern = '\\s+$', replacement = '', x = s)
}
len <- function
(
x
)
{
length(x)
}
lst <- function(
...
)
{
values = list( ... )
if(len(values) == 0) return(values)
values.names = names(values)
names = as.character(substitute(c(...))[-1])
if( is.null(values.names) )
names(values) = names
else
names(values) = iif(nchar(values.names) > 0, values.names, names)
values
}
vars2list <- function(...) {
warning('vars2list is depreciated as of Feb 29, 2016 please use lst function instead')
lst(...)
}
variable.number.arguments <- function(...) {
out = lst(...)
if( is.list(out[[1]]) && is.list(out[[1]][[1]]) )
out[[1]]
else
out
}
env <- function
(
...,
hash = TRUE,
parent = parent.frame(),
size = 29L
)
{
temp = new.env(hash = hash, parent = parent, size = size)
values = lst(...)
if(len(values) == 0) return(temp)
if(len(values) == 1 && is.environment(values[[1]]))
list2vars(values[[1]], temp)
else
list2vars(values, temp)
temp
}
env.del = function(names, env) {
warning('env.del is depreciated as of Apr 25, 2016 please use env.rm function instead')
env.rm(names, env)
}
env.rm = function(names, env) {
missing = setdiff(names, ls(env))
if( len(missing) > 0)
warning('Following names are missing in environment:', missing, '\n, names available in environment:', ls(env))
rm(list=intersect(names, ls(env)), envir=env)
}
list2vars <- function(data, env = parent.frame()) {
for(n in ls(data, all.names=T))
env[[n]] = data[[n]]
}
debug.save = function() {
gall <<- parent.frame()
}
debug.load = function() {
list2vars(gall,  parent.frame())
}
check.args = function(default.args, args=NULL) {
if(is.null(args)) return(default.args)
for(n in setdiff(ls(default.args, all.names=T), ls(args, all.names=T)))
args[[n]] = default.args[[n]]
args
}
test.equality = function(..., eps = 1e-10, na.rm=T, type=c('first', 'all')) {
values = list(...)
n = len(values)
if(n == 0) return
values.names = names(values)
names = as.character(substitute(c(...))[-1])
names = iif(nchar(values.names) > 0, values.names, names)
out = c()
if(type[1] == 'all')
for(i in 1:(n-1))
for(j in (i+1):n)
out = rbind(out, c(names[i], names[j], all(abs(values[[i]] - values[[j]]) < 1e-10, na.rm=na.rm)))
else
for(i in 2:n)
out = rbind(out, c(names[1], names[i], all(abs(values[[1]] - values[[i]]) < 1e-10, na.rm=na.rm)))
colnames(out) = spl('item1, item2,equal')
out
}
iif <- function
(
cond,
truepart,
falsepart
)
{
if(len(cond) == 1) { if(cond) truepart else falsepart }
else {
if(length(falsepart) == 1) {
temp = falsepart
falsepart = cond
falsepart[] = temp
}
if(length(truepart) == 1)
falsepart[cond] = truepart
else {
cond = ifna(cond,F)
if(requireNamespace('xts', quietly = T) && xts::is.xts(truepart))
falsepart[cond] = coredata(truepart)[cond]
else
falsepart[cond] = truepart[cond]
}
falsepart
}
}
ifna <- function
(
x,
y
)
{
return(iif(is.na(x) | is.nan(x) | is.infinite(x), y, x))
}
fast.na.omit <- function
(
x
)
{
x[!is.na(x)]
}
ifnull <- function
(
x,
y
) {
return(iif(is.null(x), y, x))
}
count <- function(
x,
side = 2
)
{
if( is.null(dim(x)) ) {
sum( !is.na(x) )
} else {
apply(!is.na(x), side, sum)
}
}
run.count <- function
(
x,
window.len
)
{
n    = length(x)
xcount = cumsum( !is.na(x) )
ycount = xcount[-c(1 : (k-1))] - c(0, xcount[-c((n-k+1) : n)])
return( c( xcount[1:(k-1)], ycount))
}
map2monthly <- function(equity)
{
if(compute.annual.factor(equity) >= 12) return(equity)
dates = index(equity)
equity = coredata(equity)
temp = as.Date(c('', 10000*date.year(dates) + 100*date.month(dates) + 1), '%Y%m%d')[-1]
new.dates = seq(temp[1], last(temp), by = 'month')
map = match( 100*date.year(dates) + date.month(dates), 100*date.year(new.dates) + date.month(new.dates) )
temp = rep(NA, len(new.dates))
temp[map] = equity
return( make.xts( ifna.prev(temp), new.dates) )
}
create.monthly.table <- function(monthly.data)
{
nperiods = nrow(monthly.data)
years = date.year(index(monthly.data[c(1,nperiods)]))
years = years[1] : years[2]
temp = matrix( double(), len(years), 12)
rownames(temp) = years
colnames(temp) = spl('Jan,Feb,Mar,Apr,May,Jun,Jul,Aug,Sep,Oct,Nov,Dec')
index = date.month(index(monthly.data[c(1,nperiods)]))
temp[] = matrix( c( rep(NA, index[1]-1), monthly.data, rep(NA, 12-index[2]) ), ncol=12, byrow = T)
return(temp)
}
load.packages <- function
(
packages,
repos = "http://cran.r-project.org",
dependencies = c("Depends", "Imports"),
...
)
{
packages = spl(packages)
for( ipackage in packages ) {
if(!require(ipackage, quietly=TRUE, character.only = TRUE)) {
install.packages(ipackage, repos=repos, dependencies=dependencies, ...)
if(!require(ipackage, quietly=TRUE, character.only = TRUE)) {
stop("package", sQuote(ipackage), 'is needed.  Stopping')
}
}
}
}
tic <- function
(
identifier
)
{
assign(paste('saved.time', identifier, sep=''), proc.time()[3], envir = .GlobalEnv)
}
toc <- function
(
identifier
)
{
if( exists(paste('saved.time', identifier, sep=''), envir = .GlobalEnv) ) {
prevTime = get(paste('saved.time', identifier, sep=''), envir = .GlobalEnv)
diffTimeSecs = proc.time()[3] - prevTime
print(paste('Elapsed time is', round(diffTimeSecs, 2), 'seconds\n'))
} else {
print('Toc error\n')
}
}
test.tic.toc <- function()
{
tic(10)
for( i in 1 : 100 ) {
temp = runif(100)
}
toc(10)
}
mlag <- function
(
m,
nlag = 1
)
{
if( is.null(dim(m)) ) {
n = len(m)
if(nlag > 0) {
m[(nlag+1):n] = m[1:(n-nlag)]
m[1:nlag] = NA
} else if(nlag < 0) {
m[1:(n+nlag)] = m[(1-nlag):n]
m[(n+nlag+1):n] = NA
}
} else {
n = nrow(m)
if(nlag > 0) {
m[(nlag+1):n,] = m[1:(n-nlag),]
m[1:nlag,] = NA
} else if(nlag < 0) {
m[1:(n+nlag),] = m[(1-nlag):n,]
m[(n+nlag+1):n,] = NA
}
}
return(m);
}
mlast = function(m, nlast=1) {
n = iif(is.null(dim(m)), len(m), nrow(m))
if(nlast >= n)
m
else
iif(is.null(dim(m)), m[(n-nlast+1):n], m[(n-nlast+1):n,,drop=F])
}
repmat <- function
(
v,
n,
m
)
{
if( is.null(dim(v)) ) v = matrix(v,1)
kronecker( matrix(1, n, m), v )
}
vec <- function
(
x
)
{
if( !is.null(dim(x)) && len(dim(x)) != 1)
dim(x) = len(x)
x
}
mat <- function
(
x,
col = T
)
{
if( is.null(dim(x)) || len(dim(x)) != 2) {
n = names(x)
if( col ) {
dim(x) = c(len(x), 1)
if( !is.null(n) ) rownames(x) = n
} else {
dim(x) = c(1, len(x))
if( !is.null(n) ) colnames(x) = n
}
}
x
}
rep.row <- function
(
m,
nr
)
{
if(nr == 1) matrix(m, 1)
else matrix(m, nr=nr, nc=len(m), byrow=T)
}
rep.col <- function
(
m,
nc
)
{
if(nc == 1) m
else {
if(is.xts(m))
make.xts(matrix(coredata(m), nr=len(m), nc=nc, byrow=F), index(m))
else
matrix(m, nr=len(m), nc=nc, byrow=F)
}
}
col.name2row = function(x, move.name=T) {
if(move.name && !is.null(colnames(x))) {
x = rbind(colnames(x), x)
colnames(x) = NULL
x
} else {
colnames(x) = x[1,]
x[-1,]
}
}
row.name2col = function(x, move.name=T) {
if(move.name && !is.null(rownames(x))) {
x = cbind(rownames(x), x)
rownames(x) = NULL
x
} else {
rownames(x) = x[,1]
x[,-1]
}
}
lookup.index <- function
(
data,
i,
details = F
)
{
n = nrow(data)
irow = ((i - 1) %% n) + 1
icol = ((i - 1) %/% n) +1
if(details)
list(irow=irow,icol=icol,obs=data[irow,icol],obsr=data[max(0,irow-5):min(nrow(data),irow+5),icol])
else
list(irow=irow,icol=icol)
}
beta.degree <- function(beta)
{
atan(beta)*360/(2*pi)
}
to.nice = function(out,nround=2,sprefix='',eprefix='') {
if( !is.null(dim(out)) ) {
temp = matrix('', nrow(out),ncol(out))
rownames(temp) = iif(is.xts(out), paste(index(out)),rownames(out))
colnames(temp) = colnames(out)
temp.n = apply(out,2,as.double)
index = is.na(temp.n)
temp[] = paste(sprefix,format(round( temp.n ,nround),big.mark=",", scientific=FALSE),eprefix ,sep='')
temp[index] = coredata(out)[index]
temp
} else {
temp.n = as.double(out)
index = is.na(temp.n)
temp = paste(sprefix,format(round( temp.n ,nround),big.mark=",", scientific=FALSE),eprefix ,sep='')
temp[index] = out[index]
temp
}
}
to.percent = function(x, nround=2) to.nice(100*x,nround,'','%')
to.cash = function(x, nround=2) to.nice(x,nround,'$')
Sys.setenv(TZ = 'GMT')
XTSFunctions <- function() {}
make.xts <- function
(
x,
order.by
)
{
tzone = Sys.getenv('TZ')
orderBy = class(order.by)
index = as.numeric(as.POSIXct(order.by, tz = tzone))
if( is.null(dim(x)) ) {
if( len(order.by) == 1 )
x = t(as.matrix(x))
else
dim(x) = c(len(x), 1)
}
x = as.matrix(x)
x = structure(.Data = x,
index = structure(index, tzone = tzone, tclass = orderBy),
class = c('xts', 'zoo'), .indexCLASS = orderBy, tclass=orderBy, .indexTZ = tzone, tzone=tzone)
if (!is.null(attributes(x)$dimnames[[1]]))
dimnames(x) <- dimnames(x)
x
}
as.xts.list <- function(data) { for(n in names(data)) colnames(data[[n]])=n; do.call(cbind, data)}
xts2ts = function(x) {
annual.factor = compute.annual.factor(x)
map = c(date.day, date.week, date.month, date.quarter)
names(map) = trim(spl('252, 52, 12, 4'))
date.fn = map[[paste(annual.factor)]]
first.date = index(first(x))
last.date = index(last(x))
start = date.year(first.date)
end = date.year(last.date)
if( !is.null(date.fn) ) {
start = c(start, date.fn(first.date))
end = c(end, date.fn(last.date))
}
ts(coredata(x[,1]), start = start, end = end, deltat = 1 / annual.factor)
}
flip.xts <- function(x)
{
dates = index(x)
dates.index = nrow(x):1
out = make.xts(coredata(x)[dates.index,], dates[dates.index])
indexClass(out) = indexClass(x)
return( out )
}
write.xts <- function
(
x,
filename,
append = FALSE,
...
)
{
cat('Date', file = filename, append = append)
write.table(x, sep=',',  row.names = format(index(x), ...),
col.names = NA, file = filename, append = T, quote = F)
}
read.xts <- function
(
x,
date.fn = paste,
index.class = 'Date',
decreasing = FALSE,
sep = ',',
date.column = 1,
skip = 0L,
...
)
{
load.packages('data.table')
if (is.matrix(x) || (is.data.frame(x) && !is.data.table(x)) ) {
data = x
dates = as.matrix(data[,date.column,drop=F])
data  = data[,-date.column,drop=F]
} else {
filename = x
if(!is.data.table(x)) {
out = fread(filename, stringsAsFactors=F, sep=sep, autostart=2, skip=skip)
setnames(out,gsub(' ', '_', trim(colnames(out))))
} else out = x
rest.columns.expr = parse(text = paste('list(', paste(names(which(sapply(out,class)[-(1:date.column)] != 'character')),collapse=','),')'))
dates = as.matrix(out[,date.column,with=FALSE])
index = which(sapply(out,class) != 'character')
index = index[ index > date.column ]
data =  as.matrix(out[,index,with=FALSE]+0)
}
dates = as.POSIXct(match.fun(date.fn)(dates), tz = Sys.getenv('TZ'), ...)
dates.index = iif(is.null(decreasing), 1:nrow(data), order(dates, decreasing = decreasing) )
out = make.xts(data[dates.index,,drop=F], dates[dates.index])
indexClass(out) = index.class
return( out )
}
read.xts.old <- function
(
filename,
date.fn = paste,
index.class = 'Date',
decreasing = FALSE,
...
)
{
out = read.csv(filename, stringsAsFactors=F)
dates = as.POSIXct(match.fun(date.fn)(out[,1]), tz = Sys.getenv('TZ'), ...)
dates.index = order(dates, decreasing = decreasing)
out = make.xts(out[dates.index,-1,drop=F], dates[dates.index])
indexClass(out) = index.class
return( out )
}
read.xts.yahoo.old <- function
(
filename,
date.fn = paste,
index.class = 'Date',
decreasing = FALSE,
...
)
{
temp = scan(filename, what=list('',double(0), double(0),double(0),double(0),double(0),double(0)), skip=1, sep=',', quiet =T)
dates = as.POSIXct(match.fun(date.fn)(temp[[1]]), tz = Sys.getenv('TZ'), ...)
dates.index = order(dates, decreasing = decreasing)
out = matrix(double(1),len(dates), 6)
colnames(out) = spl('Open,High,Low,Close,Volume,Adjusted')
out[,1] = temp[[2]]
out[,2] = temp[[3]]
out[,3] = temp[[4]]
out[,4] = temp[[5]]
out[,5] = temp[[6]]
out[,6] = temp[[7]]
out = make.xts(out[dates.index,],  dates[dates.index])
indexClass(out) = index.class
return( out )
}
read.xts.test <- function() {
load.packages('rbenchmark')
filename = 'c:/stocks/SPY.csv'
test1 <- function() {
out = read.csv(filename, stringsAsFactors=F)
}
test2 <- function() {
out1 = fread(filename, stringsAsFactors=F)
}
test3 <- function() {
out2 = scan(filename, what=list('',double(0), double(0),double(0),double(0),double(0),double(0)), skip=1, sep=',', quiet =T)
}
library(rbenchmark)
benchmark(
test1(),
test2(),
test3(),
columns = c("test", "replications", "elapsed", "relative"),
order = "relative",
replications = 20
)
test1 <- function() {
out = read.xts(filename, format = '%Y-%m-%d')
}
test2 <- function() {
out1 = read.xts.old(filename, format = '%Y-%m-%d')
}
test3 <- function() {
out2 = read.xts.yahoo.old(filename, format = '%Y-%m-%d')
}
library(rbenchmark)
benchmark(
test1(),
test2(),
test3(),
columns = c("test", "replications", "elapsed", "relative"),
order = "relative",
replications = 20
)
}
index.xts <- function
(
x
)
{
temp = attr(x, 'index')
class(temp) = c('POSIXct', 'POSIXt')
type = attr(x, '.indexCLASS')[1]
if( type == 'Date' || type == 'yearmon' || type == 'yearqtr')
temp = as.Date(temp)
return(temp)
}
index4xts <- function
(
x
)
{
temp = attr(x, 'index')
class(temp)='POSIXct'
return(temp)
}
index2date.time <- function(temp) {
class(temp)='POSIXct'
if( attr(x, '.indexCLASS')[1] == 'Date') {
as.Date(temp)
} else {
as.POSIXct(temp, tz = Sys.getenv('TZ'))
}
}
dates2index <- function( x, dates = 1:nrow(x) ) {
dates.index = dates
if(!is.numeric(dates)) {
temp = x[,1]
temp[] = 1:nrow(temp)
dates.index = as.numeric(temp[dates])
}
return(dates.index)
}
find.names <- function(names, data, return.index = T)
{
names = spl(names)
all.names = colnames(data)
out = list()
for(n in names) {
loc = grep(n, all.names, ignore.case = TRUE)
if(len(loc) == 0 && ncol(data) == 1 &&
(grepl(n,'close',ignore.case = TRUE) || grepl(n,'adjusted',ignore.case = TRUE))
) loc = 1
if(len(loc) > 0) out[[n]] = iif(return.index, loc, all.names[loc])
}
iif(len(names) == 1 && len(out) == 1, out[[1]][1], out)
}
make.stock.xts <- function(out, column=1) {
names = spl('Open,High,Low,Close,Volume,Adjusted')
names.index = find.names(names, out)
temp= list()
for(n in names)
if( !is.null(names.index[[n]]) )
temp[[n]] = out[,names.index[[n]][1]]
if(is.null(temp$Close) && is.null(temp$Adjusted))
temp$Close = temp$Adjusted = out[,column]
if(is.null(temp$Adjusted)) temp$Adjusted = temp$Close
if(is.null(temp$Close)) temp$Close = temp$Adjusted
if(is.null(temp$Open)) temp$Open = temp$Close
if(is.null(temp$High)) temp$High = temp$Close
if(is.null(temp$Low)) temp$Low = temp$Close
if(is.null(temp$Volume)) temp$Volume = 0
out = cbind(temp$Open,temp$High,temp$Low,temp$Close,temp$Volume,temp$Adjusted)
colnames(out) = names
out
}
scale.one <- function
(
x,
overlay = F,
main.index = which(!is.na(x[1,]))[1]
)
{
index = 1:nrow(x)
if( overlay )
x / rep.row(apply(x, 2,
function(v) {
i = index[!is.na(v)][1]
v[i] / as.double(x[i,main.index])
}
), nrow(x))
else
x / rep.row(apply(x, 2, function(v) v[index[!is.na(v)][1]]), nrow(x))
}
get.extension <- function(x)
{
trim( tail(spl(x,'\\.'),1) )
}
get.full.filename <- function(x)
{
trim( tail(spl(gsub('\\\\','/',x),'/'),1) )
}
get.filename <- function(x)
{
temp = spl(get.full.filename(x),'\\.')
join(temp[-len(temp)])
}
getSymbols.sit <- function
(
Symbols,
env = .GlobalEnv,
auto.assign = TRUE,
stock.folder = 'c:/temp/Seasonality/stocks',
stock.date.format = '%Y-%m-%d',
...
)
{
require(quantmod)
for(i in 1:len(Symbols)) {
s = Symbols[i]
temp = list()
temp[[ s ]] = list(src='csv', format=stock.date.format, dir=stock.folder)
setSymbolLookup(temp)
temp = quantmod::getSymbols(s, env = env, auto.assign = auto.assign)
if (!auto.assign) {
cat(s, format(range(index(temp)), '%d-%b-%Y'), '\n', sep='\t')
return(temp)
}
if(!is.null(env[[ s ]]))
cat(i, 'out of', len(Symbols), 'Reading', s, format(range(index(env[[ s ]])), '%d-%b-%Y'), '\n', sep='\t')
else
cat(i, 'out of', len(Symbols), 'Missing', s, '\n', sep='\t')
}
}
parse.expr = function(expr) {
if (is.character(expr))
expr = spl(expr)
expr = gsub('<newline>','\n', expr)
expr = trim(spl(gsub('\n', ',', join(expr, ','))))
expr = expr[nchar(expr) > 0 & substring(expr, 1, 1) != "#"]
sapply(expr, function(x) spl(x,'#')[1])
}
map.symbols = function(Symbols) {
Symbols = toupper(Symbols)
map = list()
for(s in Symbols) {
temp = spl(s, "=")
if ( len(temp) > 1 ) {
name = temp[1]
values = trim(spl(temp[2], '\\+'))
value1 = values[1]
value1.name = grepl('\\[', value1)
value1 = gsub('\\]','',gsub('\\[','',value1))
value1 = trim(spl(value1,';'))
values = values[-1]
for(n in trim(spl(name,';')))
map[[ n  ]] = c(value1[1], values)
if( len(value1) > 1 || value1.name)
for(n in value1)
map[[ n  ]] = c(n, values)
} else {
temp = spl(temp, '\\+')
name = temp[1]
values = trim(temp[-1])
for(n in trim(spl(name,';')))
map[[ n  ]] = c(n, values)
}
}
map
}
getSymbols.extra <- function
(
Symbols = NULL,
env = parent.frame(),
getSymbols.fn = getSymbols,
raw.data = new.env(),
set.symbolnames = F,
auto.assign = T,
try.extend = T,
...
)
{
Symbols = parse.expr(Symbols)
if(len(Symbols) < 1) return(Symbols)
map = map.symbols(Symbols)
Symbols = unique(unlist(map))
Symbols = setdiff(Symbols, ls(raw.data, all.names=T))
data = new.env()
if(len(Symbols) > 0) match.fun(getSymbols.fn)(Symbols, env=data, auto.assign = T, ...)
for(n in ls(raw.data, all.names=T)) data[[n]] = raw.data[[n]]
if (set.symbolnames) env$symbolnames = names(map)
for(s in names(map)) {
env[[ s ]] = data[[ gsub('\\^', '', map[[ s ]][1]) ]]
if(try.extend)
if( len(map[[ s ]]) > 1)
for(i in 2:len(map[[ s ]]))
if(is.null(data[[ gsub('\\^', '', map[[ s ]][i]) ]]))
cat('Not Downloaded, main =', s, 'missing' , gsub('\\^', '', map[[ s ]][i]), '\n', sep='\t')
else
env[[ s ]] = extend.data(env[[ s ]], data[[ gsub('\\^', '', map[[ s ]][i]) ]], scale=T)
if (!auto.assign)
return(env[[ s ]])
}
}
getSymbols.extra.test <- function()
{
tickers = spl('REIT=RWX, RWX+VNQ, REIT.LONG=RWX+VNQ+VGSIX')
data <- new.env()
getSymbols.extra(tickers, src = 'yahoo', from = '1980-01-01', env = data, auto.assign = T)
bt.start.dates(data)
data$symbolnames = spl('REIT.LONG,RWX,REIT')
for(i in data$symbolnames) data[[i]] = adjustOHLC(data[[i]], use.Adjusted=T)
bt.prep(data, align='keep.all', fill.gaps = T)
plota.matplot(data$prices)
raw.data <- new.env()
raw.data$GOLD = bundes.bank.data.gold()
tickers = spl('GLD, GLD.LONG=GLD+GOLD')
data <- new.env()
getSymbols.extra(tickers, src = 'yahoo', from = '1980-01-01', env = data, raw.data = raw.data, auto.assign = T)
bt.start.dates(data)
data$symbolnames = spl('GLD.LONG,GLD')
for(i in data$symbolnames) data[[i]] = adjustOHLC(data[[i]], use.Adjusted=T)
bt.prep(data, align='keep.all', fill.gaps = T)
plota.matplot(data$prices)
}
getSymbols.intraday <- function
(
Symbols = NULL,
env = parent.frame(),
getSymbols.fn = getSymbols,
auto.assign = T,
...
)
{
if(len(Symbols) > 0) {
match.fun(getSymbols.fn)(Symbols, env = env, auto.assign = auto.assign, ...)
Symbols = parse.expr(Symbols)
if(len(Symbols) < 1) return(Symbols)
map = map.symbols(Symbols)
map = sapply(map, first)
data.today = getQuote.yahoo.today(unique(map))
data.today[data.today == 'N/A'] = NA
lookup = 1:nrow(data.today)
names(lookup) = toupper(trim(data.today$Symbol))
data.today = data.today[lookup[map],]
data.today$Symbol = names(map)
bt.append.today(env, data.today)
}
}
getSymbols.intraday.test <- function() {
tickers = '
LQD + VWESX
DBC + CRB
VTI +VTSMX
ICF + VGSIX
CASH = SHY
'
data.intraday  = env()
getSymbols.extra(tickers, src = 'yahoo', from = '2012-01-01', env = data.intraday,
raw.data = data.proxy.raw, set.symbolnames = T, auto.assign = T,
getSymbols.fn = getSymbols.intraday)
last(data.intraday$CASH,5)
last(data.intraday$VTI,5)
data = env()
getSymbols.extra(tickers, src = 'yahoo', from = '2012-01-01', env = data,
raw.data = data.proxy.raw, set.symbolnames = T, auto.assign = T)
last(data$CASH,5)
last(data$VTI,5)
}
convert2expr = function(expr) {
if(class(substitute(expr)) == '{') {
if(F) {
return(as.expression(substitute(expr)))
} else {
temp = deparse(substitute(expr))
return(parse(text = temp[-c(1,length(temp))]))
}
}
if(is.character(expr)) return(parse(text = expr))
expr
}
convert2expr.test = function() {
convert2expr({x=2+y})
convert2expr({x=2+y; a=b})
convert2expr({
x=2+y
a=b
})
convert2expr(expression(x=2+y))
convert2expr(expression(x=2+y,a=b))
convert2expr('x=2+y')
convert2expr('x=2+y; a=b')
convert2expr('
x=2+y
a=b
')
a = convert2expr({x=2+y})
expr.symbols(a)
}
remove.operators = function(tokens) {
tokens = unique(trim( tokens[-grep('[=\\+\\-\\*/><\\(\\)\\{\\}]',tokens)] ))
tokens[nchar(tokens) > 0 & tokens != 'expression' & tokens != 'convert2expr']
}
expr.symbols = function(expr) {
if(mode(substitute(expr)) == 'call')
return(remove.operators(
c(names(as.pairlist(substitute(expr))),
all.names(substitute(expr)))
))
if(mode(substitute(expr)) == 'name')
if(is.expression(expr)) {
return(remove.operators(
c(names(as.pairlist(expr)),
all.names(expr))
))
} else {
return(remove.operators(
all.names(expr)
))
}
if(is.character(expr))
return(remove.operators(
all.names(parse(text=expr))
))
expr
}
expr.symbols.test = function() {
expr.symbols({
x = 2 + y+z
a=2+b
})
expr.symbols({x=y+2})
expr.symbols('x=y+2')
expr.symbols(expression(x=2+y))
expr.symbols(expression(x=2+y, a=b))
a = expression(x=2+y+zzasd)
expr.symbols(a)
}
log.fn <- function(p.start=0, p.end=1) {
p.start = p.start
p.end = p.end
function(..., percent=NULL) {
cat(..., iif(is.null(percent),'',paste(', percent = ', round(100 * (p.start + percent * (p.end - p.start)), 1), '%', sep='')), '\n')
}
}
log.fn.msg <- function(msg, log = log.fn()) {
log = log
msg = msg
function(..., percent=NULL) { log(paste(msg, ...), percent=percent) }
}
asc <- function(x) { strtoi(charToRaw(x),16L) }
chr <- function(n) { rawToChar(as.raw(n)) }
make.random.string <- function(nbits = 256) { chr( runif(nbits/8, 1, 255) ) }
random.string <- function(lenght = 12) { join(sample(c(0:9, letters, LETTERS),lenght, replace=TRUE)) }
ls.f <- function(env=sys.frame(-1))unlist(lapply(ls(env=env),function(x)if(is.function(get(x)))x))
ls.v <- function(env=sys.frame(-1))unlist(lapply(ls(env=env),function(x)if(!is.function(get(x)))x))
parse.number <- function(x) {
as.numeric(gsub('[^0-9\\+-\\.]', '', x) )
}
map2vector = function(expr, labels, default = 0) {
if( is.xts(labels) ) labels = names(labels)
if( is.character(labels) ) labels = spl(labels)
n = iif( is.numeric(labels), labels, len(labels) )
if( len(expr) == 0 ) return(rep(default, n))
if(mode(substitute(expr)) == 'call') {
value = rep(default, n)
names(value) = labels
e = evalq(environment(), as.list(value))
eval(substitute(expr), e)
temp = unlist(as.list(e))
value[names(temp)] = temp
return(value)
}
if( is.list(expr) ) {
out = rep(default, n)
out[ match(toupper(names(expr)), toupper(labels)) ] = unlist(expr)
out
} else
ifna(iif( len(expr) == 1, rep(expr, n), expr), default)
}
rev.map = function(map) {
value = names(map)
names(value) = map
value
}
write.file = function(..., file) cat(..., file=file)
read.file = function(file) readChar(file, file.info(file)$size)
string.buffer = function() structure(list(file = rawConnection(raw(0L), open='w')), class = 'StringBuffer')
add = function(x,...,sep,end.sep) UseMethod('add',x)
add.StringBuffer = function(x,...,sep=',',end.sep='\n') {
cat(..., file = x$file, sep = sep)
if(nchar(end.sep) > 0) cat(end.sep, file = x$file)
}
string = function(x) UseMethod('string',x)
string.StringBuffer = function(x) rawToChar(rawConnectionValue(x$file))
close = function(x) UseMethod('close',x)
close.StringBuffer = function(x) {close(x$file); x$file = NULL}
string.buffer.test = function() {
file = rawConnection(raw(0L), open="w")
write('asbcd', file)
write('234543', file)
res =rawToChar(rawConnectionValue(file))
close(file)
file = NULL;
sb = string.buffer()
add(sb, 'asbcd')
add(sb, '234543')
string(sb)
close(sb)
sb=NULL
test.base = function() {
s =''
for(i in 1:10000)
s = paste(s,'abcdef',sep='')
nchar(s)
}
test.string.buffer = function() {
sb = string.buffer()
for(i in 1:10000)
add(sb, 'abcdef', '')
s = string(sb)
sb=NULL
nchar(s)
}
library(rbenchmark)
benchmark(
test.base(),
test.string.buffer(),
columns = c("test", "replications", "elapsed", "relative"),
order = "relative",
replications = 1
)
}
.onLoad <- function(libname = find.package("SIT"), pkgname = "SIT") {
Sys.setenv(TZ = 'GMT')
library(grDevices)
plota.control <<- new.env()
plota.control$col.border = 'black'
plota.control$col.up = 'green'
plota.control$col.dn = 'red'
plota.control$col.x.highlight = 'orange'
plota.control$col.y.highlight = 'orange'
plota.control$xaxis.ticks = c()
plota.theme.green.orange();
}
len = function(x) length(x)
mlast = function(m, n=1) m[(nrow(m)-n+1), ,drop=F]
spl = function(s, delim = ',') unlist(strsplit(s,delim))
rep.row = function(m, nr) matrix(m, nr=nr, nc=len(m), byrow=T)
trim = function(s) {
s = sub(pattern = '^\\s+', replacement = '', x = s)
sub(pattern = '\\s+$', replacement = '', x = s)
}
index.xts = function(x) {
temp = attr(x, 'index')
class(temp) = c('POSIXct', 'POSIXt')
type = attr(x, '.indexCLASS')[1]
if( type == 'Date' || type == 'yearmon' || type == 'yearqtr')
temp = as.Date(temp)
return(temp)
}
custom.date.bus = function(expr, dates, calendar = NULL) {
apply.business.days(dates, function(x) custom.date(expr, x), calendar)
}
custom.date = function(expr, dates) {
if( xts::is.xts(dates) ) dates = index(dates)
dates = as.Date(dates)
expr = gsub('the ', '', tolower(expr))
expr = gsub(' in every ', ' every ', expr)
expr = gsub(' of every ', ' every ', expr)
tokens = trim(spl(spl(spl(expr, ' in '), ' every '), ' of '))
stack = list(
splits = date.all(dates),
dates.index = 1:len(dates)
)
selected = rep.row(c(1,len(dates)), len(dates))
selected.n = 1
for(token in rev(tokens[nchar(tokens) > 0])) {
selected0 = selected[1:selected.n, , drop=F]
selected.n = 0
for(i in 1:nrow(selected0)) {
temp = custom.date.token(token, dates, stack, selected0[i,])
selected[(selected.n+1):(selected.n+nrow(temp)),] = temp
selected.n = selected.n + nrow(temp)
}
}
selected[1:selected.n,1]
}
custom.date.token = function(expr, dates, stack, selected) {
tokens = trim(spl(tolower(expr), ' '))
tokens = tokens[nchar(tokens) > 0]
n.tokens = len(tokens)
freq = custom.date.map(tokens[n.tokens])
periods = date.periodicity.map(freq$freq)
if( is.null(periods) )
warning('unknown freq', freq$freq)
if( periods == 'days' ) {
temp = cbind(selected[1]:selected[2], selected[1]:selected[2])
rownames(temp) =  stack$splits$dayofweek[selected[1]:selected[2]]
} else
temp = custom.date.extract(selected[1], selected[2], periods, stack)
if( !is.null(freq$pick) )
temp = temp[rownames(temp) == freq$pick,,drop=F]
if( n.tokens == 1 ) return(temp)
if( n.tokens == 2 ) {
if( tokens[1] == 'last' )
return(mlast(temp))
if( tokens[1] == 'first' )
return(temp[1,,drop=F])
}
offset = stringr::str_match(tokens[1], '^[0-9]+')[1]
if( is.na(offset) ) warning('unknown offset', tokens[1])
offset = as.numeric(offset)
if( offset > nrow(temp) ) {
if( n.tokens == 2 )
mlast(temp)
else
temp[1,,drop=F]
} else {
if( n.tokens == 2 )
temp[offset,,drop=F]
else
temp[nrow(temp)-offset,,drop=F]
}
}
month.map.abbr = 1:12
names(month.map.abbr) = spl('jan,feb,mar,apr,may,jun,jul,aug,sep,oct,nov,dec')
month.map = 1:12
names(month.map) = spl('january,february,march,april,may,june,july,august,september,october,november,december')
day.map.abbr = 0:6
names(day.map.abbr) = spl('sun,mon,tue,wed,thu,fri,sat')
day.map = 0:6
names(day.map) = spl('sunday,monday,tuesday,wednesday,thursday,friday,saturday')
custom.date.map = function(token) {
if( !is.na(month.map.abbr[token]) )
return(list(freq='months', pick = month.map.abbr[token]))
if( !is.na(month.map[token]) )
return(list(freq='months', pick = month.map[token]))
if( !is.na(day.map.abbr[token]) )
return(list(freq='day', pick = day.map.abbr[token]))
if( !is.na(day.map[token]) )
return(list(freq='day', pick = day.map[token]))
match = stringr::str_match(token, '^q[1-4]{1}')[1]
if( !is.na(match) ) return(list(freq='quarter', pick = substring(match,2,2) ))
match = stringr::str_match(token, '^m[0-9]+')[1]
if( !is.na(match) ) return(list(freq='month', pick = substring(match,2,3) ))
match = stringr::str_match(token, '^s[1-2]{1}')[1]
if( !is.na(match) ) return(list(freq='semiannual', pick = substring(match,2,2) ))
list(freq=token)
}
custom.date.extract = function(i0, i1, freq, stack) {
label = stack$splits[[freq]][i0 : i1]
label.n = len(label)
temp = unique(c( 0, stack$dates.index[1:label.n][diff( label ) != 0], label.n ))
temp.n = len(temp)
temp = cbind(1 + temp[1:(temp.n - 1)], temp[2:temp.n])
rownames(temp) =  label[temp[,1]]
(i0 - 1) + temp
}
custom.date.test = function () {
dates = seq(Sys.Date()-1000, Sys.Date(), 1)
dates[custom.date('last day in Apr', dates)]
dates[custom.date('first day in Apr', dates)]
dates[custom.date('last day in first week in Apr', dates)]
dates[custom.date('last Mon in Apr', dates)]
dates[custom.date('last Fri in Apr', dates)]
dates[custom.date('first day in Apr', dates)]
dates[custom.date('1st day in Apr', dates)]
dates[custom.date('10th day in Apr', dates)]
dates[custom.date('50th day in Apr', dates)]
dates[custom.date('10th to last day in Apr', dates)]
dates[custom.date('3rd Mon in Q', dates)]
dates[custom.date('3rd Mon in 1st Q', dates)]
dates[custom.date('3rd Mon in Q1', dates)]
dates[custom.date('3rd Mon in last M in Q1', dates)]
format(dates[custom.date('3rd Fri in Q', dates)], '%Y %b %d %w')
dates = seq(as.Date('1-jan-2010','%d-%b-%Y'), as.Date('29-apr-2015','%d-%b-%Y'), 1)
dates[custom.date('last day in Apr', dates)]
dates = seq(as.Date('1-jan-2010','%d-%b-%Y'), as.Date('30-apr-2015','%d-%b-%Y'), 1)
dates[custom.date('last day in Apr', dates)]
dates = seq(as.Date('1-jan-2010','%d-%b-%Y'), as.Date('29-apr-2015','%d-%b-%Y'), 1)
dates[custom.date.bus('last day in Apr', dates)]
dates = seq(as.Date('1-jan-2010','%d-%b-%Y'), as.Date('20-oct-2015','%d-%b-%Y'), 1)
dates[custom.date('last day in Apr', dates)]
expect_identical(
dates[custom.date('last day in Apr', dates)],
as.Date(
c("2010-04-30", "2011-04-30", "2012-04-30", "2013-04-30", "2014-04-30", "2015-04-30")
,'%Y-%m-%d')
)
}
custom.date.debug = function () {
dates[custom.date('last day in first week in Apr', dates)]
dates = seq(Sys.Date()-1000, Sys.Date(), 1)
custom.date('last day in Apr', dates)
custom.date('3rd Mon in 1st Q', dates)
custom.date('Mon in 3rd W in 1st Q', dates)
dates = seq(Sys.Date()-1000, Sys.Date(), 1)
stack = env(splits, dates.index = 1:len(dates))
i0 = 1
i1 = len(dates)
freq = 'month'
temp = custom.date.extract(i0, i1, freq, stack)
temp = temp[rownames(temp) == month.map.abbr['apr'],]
dates[temp[,1]]
dates[temp[,2]]
i0 = temp[1,1]
i1 = temp[1,2]
freq = 'week'
temp1 = custom.date.extract(i0, i1, freq, stack)
temp1 = temp1[1,,drop=F]
dates[temp1[,1]]
dates[temp1[,2]]
}
date.dayofweek <- function(dates) {
as.POSIXlt(dates)$wday
}
date.day <- function(dates) {
as.POSIXlt(dates)$mday
}
date.week <- function(dates) {
dates = as.POSIXlt(dates)
offset = (7 + dates$wday - dates$yday %% 7) %%7
(dates$yday +  offset)%/% 7
}
date.month <- function(dates) {
as.POSIXlt(dates)$mon + 1
}
quarter.map = c(1,1,1,2,2,2,3,3,3,4,4,4)
date.quarter <- function(dates) {
quarter.map[date.month(dates)]
}
semiannual.map = c(1,1,1,1,1,1,2,2,2,2,2,2)
date.semiannual = function (dates) {
semiannual.map[date.month(dates)]
}
date.year = function (dates) {
as.POSIXlt(dates)$year + 1900
}
date.all = function(dates)
{
dates = as.POSIXlt(dates)
offset = (7 + dates$wday - dates$yday %% 7) %%7
list(
dayofweek = dates$wday,
mday = dates$mday,
yday = dates$yday,
weeks = (dates$yday +  offset)%/% 7,
months = dates$mon + 1,
quarters = quarter.map[dates$mon + 1],
semiannual = semiannual.map[dates$mon + 1],
years = dates$year + 1900
)
}
date.period.test = function() {
date.dayofweek0 <- function(dates) {
return(as.double(format(dates, '%w')))
}
date.day0 <- function(dates) {
return(as.double(format(dates, '%d')))
}
date.week0 <- function(dates) {
return(as.double(format(dates, '%U')))
}
date.month0 <- function(dates) {
return(as.double(format(dates, '%m')))
}
date.quarter0 <- function(dates) {
(((date.month(dates))-1) %/% 3)+1
}
date.year0 = function (dates) {
return(as.double(format(dates, '%Y')))
}
dates1 = seq(Sys.Date()-100000, Sys.Date(), 1)
all.equal(diff(date.week0(dates1))!=0 , diff(date.week(dates1))!=0 )
dates = seq(Sys.Date()-10000, Sys.Date(), 1)
library(rbenchmark)
benchmark(
test1 = diff(date.week0(dates))!=0,
test2 = diff(date.week(dates))!=0,
columns = c("test", "replications", "elapsed", "relative"),
order = "relative",
replications = 200
)
all.equal(diff(date.dayofweek0(dates1))!=0 , diff(date.dayofweek(dates1))!=0 )
benchmark(
test1 = diff(date.dayofweek0(dates))!=0,
test2 = diff(date.dayofweek(dates))!=0,
columns = c("test", "replications", "elapsed", "relative"),
order = "relative",
replications = 200
)
all.equal(diff(date.day0(dates1))!=0 , diff(date.day(dates1))!=0 )
benchmark(
test1 = diff(date.day0(dates))!=0,
test2 = diff(date.day(dates))!=0,
columns = c("test", "replications", "elapsed", "relative"),
order = "relative",
replications = 200
)
all.equal(diff(date.month0(dates1))!=0 , diff(date.month(dates1))!=0 )
benchmark(
test1 = diff(date.month0(dates))!=0,
test2 = diff(date.month(dates))!=0,
columns = c("test", "replications", "elapsed", "relative"),
order = "relative",
replications = 200
)
all.equal(diff(date.quarter0(dates1))!=0 , diff(date.quarter(dates1))!=0 )
benchmark(
test1 = diff(date.quarter0(dates))!=0,
test2 = diff(date.quarter(dates))!=0,
columns = c("test", "replications", "elapsed", "relative"),
order = "relative",
replications = 200
)
all.equal(diff(date.year0(dates1))!=0 , diff(date.year(dates1))!=0 )
benchmark(
test1 = diff(date.year0(dates))!=0,
test2 = diff(date.year(dates))!=0,
columns = c("test", "replications", "elapsed", "relative"),
order = "relative",
replications = 200
)
}
date.week.ends <- function(dates, last.date=T)
{
ends = which(diff( 100*date.year(dates) + date.week(dates) ) != 0)
ends.add.last.date(ends, len(dates), last.date)
}
date.month.ends <- function(dates, last.date=T)
{
ends = which(diff( 100*date.year(dates) + date.month(dates) ) != 0)
ends.add.last.date(ends, len(dates), last.date)
}
date.quarter.ends <- function(dates, last.date=T)
{
ends = which(diff( 10*date.year(dates) + date.quarter(dates) ) != 0)
ends.add.last.date(ends, len(dates), last.date)
}
date.semiannual.ends = function(dates, last.date=T)
{
ends = which(diff( 10*date.year(dates) + date.semiannual(dates) ) != 0)
ends.add.last.date(ends, len(dates), last.date)
}
date.year.ends <- function(dates, last.date=T)
{
ends = which(diff( date.year(dates) ) != 0)
ends.add.last.date(ends, len(dates), last.date)
}
ends.add.last.date = function(ends, last.date, action=T)
{
if(action)
unique(c(ends, last.date))
else
ends
}
date.periodicity.map = function(periodicity) {
switch(periodicity,
days = 'days',
day = 'days',
daily = 'days',
d = 'days',
weeks = 'weeks',
week = 'weeks',
weekly = 'weeks',
w = 'weeks',
months = 'months',
month = 'months',
monthly = 'months',
m = 'months',
quarters = 'quarters',
quarter = 'quarters',
quarterly = 'quarters',
q = 'quarters',
semiannual = 'semiannual',
semiannually = 'semiannual',
s = 'semiannual',
years = 'years',
year = 'years',
yearly = 'years',
annual = 'years',
annually = 'years',
y = 'years',
NULL)
}
date.ends.fn = function(periodicity) {
switch(date.periodicity.map(periodicity),
weeks = date.week.ends,
months = date.month.ends,
quarters = date.quarter.ends,
semiannual = date.semiannual.ends,
years = date.year.ends,
NULL)
}
apply.business.days = function(dates, dates.fn = NULL, calendar = NULL, base = T) {
if( xts::is.xts(dates) ) {
dates = index(dates)
apply.business.days.internal(dates, dates.fn, calendar, base)
} else {
ok.index = business.days(dates = dates, calendar = calendar, return.index=T)
index = apply.business.days.internal(dates[ok.index], dates.fn, calendar, base)
if(base)
(1:len(dates))[ok.index][index]
else
index
}
}
apply.business.days.internal = function(dates, dates.fn = NULL, calendar = NULL, base = T) {
if( xts::is.xts(dates) ) dates = index(dates)
dates = as.Date(dates)
n = len(dates)
holidays = NULL
if( !is.null(calendar) ) {
if( requireNamespace('RQuantLib', quietly = T) )
holidays = RQuantLib::getHolidayList(calendar, dates[1] - 60, dates[1] - 1)
else
warning('RQuantLib could not be loaded')
}
before = business.days(dates[1] - 60, dates[1] - 1, holidays)
n.before = len(before)
if( !is.null(holidays) )
holidays = RQuantLib::getHolidayList(calendar, dates[n] + 1, dates[n] + 60)
after = business.days(dates[n] + 1, dates[n] + 60, holidays)
dates = c(before, dates, after)
if( !is.null(dates.fn) )
index = dates.fn(dates)
else
index = 1:len(dates)
if( base ) {
index = index[index > n.before & index <= (n.before + n)]
index = index - n.before
return(index)
}
original.dates.index = (n.before + 1) : (n.before + n)
temp.cum = cumsum(rep(1, len(dates)))
temp = temp.cum * NA
temp[index] = temp.cum[index]
days.since = temp.cum - ifna.prev(temp)
days.till = temp[ifna.prevx.rev(temp)] - temp.cum
list(days.since = days.since[original.dates.index], days.till = days.till[original.dates.index])
}
business.days.location.end = function(dates, dates.fn = date.month.ends, calendar = NULL) {
apply.business.days(dates, dates.fn, calendar, F)
}
date.ends.index <- function(out, timing) {
if(timing <= 0)
which(out$days.till == (-timing))
else
which(out$days.since == (timing))
}
date.ends = function(dates, periodicity, by=1, skip=0, last.date=T, calendar = NULL) {
periodicity = trim(tolower(periodicity))
bi.flag = substr(periodicity,1,2) == 'bi'
if(bi.flag) periodicity = substr(periodicity,3,1000)
by = if(bi.flag) 2 else by
periodicity = trim(gsub('-','',periodicity))
ends = apply.business.days(dates, calendar = calendar,
dates.fn = function(x) {
fn = date.ends.fn(periodicity)
if( is.null(fn) )
xts::endpoints(xts::xts(1:len(x), x), periodicity)
else
fn(x, last.date=F)
})
if( xts::is.xts(dates) ) dates = index(dates)
ends = ends.add.last.date(ends, len(dates), last.date)
if( skip > 0) ends = ends[-c(1:skip)]
if( by > 1) ends = ends[seq(1, len(ends), by=by)]
ends
}
date.end <- function(date = Sys.Date(), periodicity = 'months', date.format = '%Y-%m-%d') {
date = as.Date(paste(date), date.format)
temp = seq(date, date + 40, 1)
temp[date.ends.fn(periodicity)(temp)[1]]
}
business.days = function(
from = Sys.Date(),
to = as.Date(from) + 31,
holidays = NULL,
dates = NULL,
calendar = NULL,
return.index = F
) {
if( is.null(dates) )
dates = seq(as.Date(from), as.Date(to), 1)
else if( xts::is.xts(dates) )
dates = index(dates)
dates = as.Date(dates)
rm.index = date.dayofweek(dates) == 6 | date.dayofweek(dates) == 0
if( !is.null(calendar) ) {
if( requireNamespace('RQuantLib', quietly = T) )
holidays = RQuantLib::getHolidayList(calendar, dates[1], dates[len(dates)])
else
warning('RQuantLib could not be loaded')
}
if( !is.null(holidays) ) {
holidays = as.Date(holidays)
rm.index = rm.index | !is.na(match(dates, holidays))
}
if( return.index )
!rm.index
else
dates[!rm.index]
}
business.days.till.end <- function(from, holidays = NULL, fn.ends = date.month.ends) {
from = as.Date(from)
dates = business.days(from - 10, from, holidays)
from = dates[len(dates)]
dates = business.days(from, from + 40, holidays)
index = match.fun(fn.ends)(dates, F)
index[1] - 1
}
business.days.since.end <- function(from, holidays = NULL, fn.ends = date.month.ends) {
from = as.Date(from)
dates = business.days(from - 10, from, holidays)
from = dates[len(dates)]
dates = business.days(from - 40, from + 10, holidays)
index = match.fun(fn.ends)(dates, F)
last.index = index[len(index)]
if( dates[last.index] == from) return(0)
from.index = sum(dates <= from)
if( dates[last.index] < from) return(from.index - last.index)
last.index = index[(len(index) - 1)]
return(from.index - last.index)
}
next.business.day <- function(from, holidays = NULL, offset = 0) {
from = as.Date(from)
dates = business.days(from + offset, from + 10, holidays)
dates[1]
}
last.business.day <- function(from, holidays = NULL, offset = 0) {
from = as.Date(from)
dates = business.days(from - 10, from - offset, holidays)
dates[1]
}
third.friday.month <- function(years, months)
{
helper <- function(year, month) {
day = date.dayofweek( as.Date(c('', 10000*year + 100*month + 1), '%Y%m%d')[-1] )
day = c(20,19,18,17,16,15,21)[1 + day]
as.Date(c('', 10000*year + 100*month + day), '%Y%m%d')[-1]
}
if(len(years) > 1 && len(months) > 1) {
out = c()
for(month in months)
out = c(out, helper(years,month))
as.Date(out)
} else
helper(years,months)
}
map.spx.expiration <- function(data, backfill = T, offset = 0) {
dates = as.Date(index(data))
years = date.year(range(dates))
friday = third.friday.month(years[1]:(years[2]+1), 1:12)
friday.future = friday[friday > dates[len(dates)]]
friday = friday[friday <= dates[len(dates)]]
key.date.index = match(friday, dates)
na.index = which(is.na(key.date.index))
if(backfill && len(na.index)>0)
key.date.index[na.index] = match(friday[na.index]-1, dates)
if(offset != 0) {
friday = c(dates[key.date.index], friday.future)
offset.date = friday - offset
key.date.index = match(offset.date, dates)
}
key.date.index = na.omit(key.date.index)
key.date = NA * data[,1]
key.date[key.date.index,] = T
key.date
}