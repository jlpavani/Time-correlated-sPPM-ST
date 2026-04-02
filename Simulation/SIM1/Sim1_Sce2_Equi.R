# **************************************************************************** #
# Scenario 1 - Partition constant over time
# **************************************************************************** #
library(igraph); library(GIGrvg); library(Rcpp); library(RcppArmadillo)
library(mvtnorm); library(RcppDist)

sourceCpp('Code/STPPM.cpp'); source("Code/utils_STPPM.R")

load("DtaSIM1.RData")

narea = nrow(dta$X); ntime = ncol(dta$X)
seasons = dta$season; seasons_id = unique(seasons); nseason = length(seasons_id)

# **************************************************************************** #
partition = matrix(NA, narea, nseason)
partition[,c(1,5,9,13,17)] = c(3,3,3,3,3,1,1,1,1,2,2,1,1,1,1,3,1,1,1,3,3,1,3,3,3,3,4,3,4,3,4,4,4,3,4,
                               4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,2,2,2,2,2,2,2,2,2,2,4,4,4,4,4,4)
partition[,c(2,6,10,14,18)] = c(2,2,2,2,2,1,1,1,1,1,1,1,2,2,1,2,1,1,1,2,2,2,2,2,2,2,3,2,3,2,3,3,3,2,3,
                                3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,1,1,1,1,1,1,1,1,1,1,3,3,3,3,3,2)
partition[,c(3,7,11,15,19)] = c(2,2,2,1,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,1,2,1,2,2,2,1,2,
                                2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2)
partition[,c(4,8,12,16,20)] = c(2,2,2,2,2,1,1,1,1,1,1,1,1,2,1,2,1,1,1,2,2,2,1,2,2,2,1,1,1,2,1,1,1,1,1,
                                1,1,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,1,2,2,2,2)
ncl = as.numeric(apply(partition, 2, max)); rho = (ncl - 1)/(narea-1)

theta_star = c(1,3,5,7, 1,3,6, 5,7, 1,3, 1,3,5,7, 1,3,6, 5,7, 1,3,
               1,3,5,7, 1,3,6, 5,7, 1,3, 1,3,5,7, 1,3,6, 5,7, 1,3,
               1,3,5,7, 1,3,6, 5,7, 1,3)
theta = matrix(NA, narea, ntime); l = 0
for(ss in 1:nseason){
  for(j in 1:ncl[ss]){
    l = l + 1
    theta[which(partition[,ss] == j),which(seasons == seasons_id[ss])] = theta_star[l]
  }
}

lambda = dta$O * exp(dta$X[,,1] * dta$beta[1] + dta$X[,,2] * dta$beta[2]) * theta

# **************************************************************************** #
shape_theta = rate_theta = 1
shape_zeta = rate_zeta = 1

mu_beta = c(0,0); sigma_beta = diag(10, 2)
mu_delta = c(0,0,0); sigma_delta = diag(10, 3)

shape_upsilon = 10; rate_upsilon = 1
shape_kappa = 100; rate_kappa = 1

niter <- 10e+03; nburn <- 0.7*niter; nthin <- (niter - nburn)/1000

O = dta$O[,1:156]; X = dta$X[,1:156,]; V = dta$V[,1:12,]; lambda = lambda[,1:156]; ntime = 156
season = as.numeric(dta$season[1:156]); season[which(season == 202101)] <- 202004

for(j in 1:100){
  
  set.seed(j); Y = matrix(rpois(narea*ntime, lambda = lambda), narea, ntime)
  
  time_start = Sys.time()
  result = STPPM_PIG(Y, O, X, V, season, dta$adj_mat, qq = 1,
                     shape_theta, rate_theta, shape_upsilon, rate_upsilon,
                     shape_kappa, rate_kappa, shape_zeta, rate_zeta,
                     mu_beta, sigma_beta, mu_delta, sigma_delta, niter, nburn, nthin)
  time_end = Sys.time()
  time = time_end - time_start
  result$time = time
  
  result$true_values <- Y
  save(result, file = paste0("Res_st1_sc2",j,"PIG.RData"))
  
  time_start = Sys.time()
  result = STPPM_Poi(Y, O, X, season, dta$adj_mat, qq = 1,
                     shape_theta, rate_theta, shape_upsilon, rate_upsilon,
                     shape_kappa, rate_kappa, shape_zeta, rate_zeta,
                     mu_beta, sigma_beta, niter, nburn, nthin)
  time_end = Sys.time()
  time = time_end - time_start
  result$time = time
  
  result$true_values <- Y
  save(result, file = paste0("Res_st1_sc2",j,"Poi.RData"))
}