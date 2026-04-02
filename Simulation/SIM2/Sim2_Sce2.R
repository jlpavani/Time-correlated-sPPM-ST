# **************************************************************************** #
# Scenario 2 - The number of clusters remains consistent, with 4, 3, 1, and 2 
# clusters for summer, autumn, winter, and spring, respectively. Great variation on theta
# **************************************************************************** #
library(igraph); library(GIGrvg); library(Rcpp); library(RcppArmadillo)
library(mvtnorm); library(RcppDist); library(statmod)

sourceCpp('Code/STPPM.cpp'); source("Code/utils_STPPM.R")

load("DtaSIM2.RData")

narea = nrow(dta$X); ntime = ncol(dta$X)
seasons = dta$season; seasons_id = unique(seasons); nseason = length(seasons_id)

# **************************************************************************** #
partition = matrix(NA, narea, nseason)
partition[,c(1,5,9,13,17)] = c(3,3,3,3,3,1,1,1,1,2,2,1,1,1,1,3,1,1,1,3,3,1,3,3,3,3,4,3,4,3,4,4,4,3,4,
                               4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,2,2,2,2,2,2,2,2,2,2,4,4,4,4,4,4)
partition[,c(2,6,10,14,18)] = c(2,2,2,2,2,1,1,1,1,1,1,1,2,2,1,2,1,1,1,2,2,2,2,2,2,2,3,2,3,2,3,3,3,2,3,
                                3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,1,1,1,1,1,1,1,1,1,1,3,3,3,3,3,2)
partition[,c(3,7,11,15,19)] = rep(1, narea)
partition[,c(4,8,12,16,20)] = c(2,2,2,1,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,1,2,1,2,2,2,1,2,
                                2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2)
ncl = as.numeric(apply(partition, 2, max)); rho = (ncl - 1)/(narea-1)

theta_star = matrix(NA, max(ncl), nseason)
theta_star[1:ncl[1], c(1,5,9,13,17)] = c(1,25,45,10)
theta_star[1:ncl[2], c(2,6,10,14,18)] = c(1,25,10)
theta_star[1:ncl[3], c(3,7,11,15,19)] = 1
theta_star[1:ncl[4], c(4,8,12,16,20)] = c(1,10)

theta = matrix(NA, narea, ntime)
for(ss in 1:nseason){
  for(j in 1:ncl[ss]){
    theta[which(partition[,ss] == j),which(seasons == seasons_id[ss])] = theta_star[j, ss]
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

niter <- 1e+04; nburn <- 0.7*niter; nthin <- (niter - nburn)/1000

season = as.numeric(dta$season); seasons_id = unique(season); nseason = length(seasons_id)
z = matrix(1, narea, nseason); zt = matrix(1, narea, ntime); psi = dta$psi
for(j in 1:nsim){
  
  for(i in 1:narea){
    for(ss in 1:nseason){
      set.seed(j*i*ss); z[i, ss] = statmod::rinvgauss(1, mean=1, shape=NULL, dispersion=1/psi[i,ss])
      zt[i, which(season == seasons_id[ss])] = z[i,ss]
    }
  } 
  set.seed(j); Y = matrix(rpois(narea*ntime, lambda = lambda*zt), narea, ntime)
  
  for(qq in 5:1){
    time_start = Sys.time()
    result = STPPM_PIG(Y, dta$O, dta$X, dta$V, season, dta$adj_mat, qq = qq,
                       shape_theta, rate_theta, shape_upsilon, rate_upsilon,
                       shape_kappa, rate_kappa, shape_zeta, rate_zeta,
                       mu_beta, sigma_beta, mu_delta, sigma_delta, niter, nburn, nthin)
    time_end = Sys.time()
    time = time_end - time_start
    result$time = time
    
    result$true_values <- Y; result$true_z <- z
    save(result, file = paste0("Res_st2_sce2_q",qq,"_",j,"_PIG.RData"))
  }
  
  time_start = Sys.time()
  result = STPPM_PIG_iid(Y, dta$O, dta$X, dta$V, season, dta$adj_mat, shape_theta, 
                         rate_theta, shape_upsilon, rate_upsilon, shape_kappa, rate_kappa, 
                         mu_beta, sigma_beta, mu_delta, sigma_delta, niter, nburn, nthin)
  time_end = Sys.time()
  time = time_end - time_start
  result$time = time
  
  result$true_values <- Y; result$true_z <- z
  save(result, file = paste0("Res_st2_sce2_iid_",j,"_PIG.RData"))
}