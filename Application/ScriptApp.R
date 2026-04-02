# **************************************************************************** #
# Application
# **************************************************************************** #
library(igraph); library(GIGrvg); library(Rcpp); library(RcppArmadillo); library(mvtnorm); library(RcppDist)
library(truncnorm); library(abind); library(statmod)

sourceCpp('Code/STPPM.cpp'); source("Code/utils_STPPM.R")

load("dtaApp.RData")

narea = nrow(dta$X); ntime = ncol(dta$X)
seasons = as.numeric(dta$season); seasons_id = unique(seasons); nseason = length(seasons_id)

shape_theta = rate_theta = 1
shape_zeta = rate_zeta = 1

nbeta = dim(dta$X)[3]; mu_beta = rep(0,nbeta); sigma_beta = diag(10, nbeta)
ndelta = dim(dta$V)[3]; mu_delta = rep(0,ndelta); sigma_delta = diag(10, ndelta)

shape_upsilon = 10; rate_upsilon = 1
shape_kappa = 100; rate_kappa = 1

niter <- 15e+03; nburn <- 2*niter/3; nthin <- (niter - nburn)/1000

# Dependent
for(j in 1:12){
  time_start = Sys.time()
  result = STPPM_PIG(dta$Y, dta$O, dta$X, dta$V, seasons, dta$adj_mat, qq = j,
                     shape_theta, rate_theta, shape_upsilon, rate_upsilon,
                     shape_kappa, rate_kappa, shape_zeta, rate_zeta,
                     mu_beta, sigma_beta, mu_delta, sigma_delta, niter, nburn, nthin)
  time_end = Sys.time()
  time = time_end - time_start
  result$time = time
    
  save(result, file = paste0("Res_q",j,".RData"))
}

# Independent
time_start = Sys.time()
result = STPPM_PIG_iid(dta$Y, dta$O, dta$X, dta$V, seasons, dta$adj_mat,
                       shape_theta, rate_theta, shape_upsilon, rate_upsilon,
                       shape_kappa, rate_kappa, mu_beta, sigma_beta, mu_delta, sigma_delta,
                       niter, nburn, nthin)
time_end = Sys.time()
time = time_end - time_start
result$time = time

save(result, file = "Res_iid.RData")