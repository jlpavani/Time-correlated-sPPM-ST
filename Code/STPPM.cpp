// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppDist)]]

#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h> 
#include <Rcpp.h>
#include <stdio.h>
#include <iostream>
#include <armadillo>
#include <mvnorm.h>
#include <truncnorm.h>
#include <vector>

using namespace Rcpp;
using namespace std;
using namespace arma;

// Check if k belongs to the vector vec
bool belong(arma::uword k, arma::uvec vec){
  
    double su = 0.0;
	int si = vec.size();
    for(int i = 0; i < si; i++){
		if(k == vec(i)){
			su += 1;
		}
	}
    if(su != 0){
        return true;
    }else{
        return false;
    }
}

// Logposterior distribution of beta
double logpost_beta(arma::mat Y, arma::mat O, arma::cube X, arma::mat theta_t, arma::mat z_t, arma::colvec betaa, arma::colvec mu_beta, arma::mat sigma_beta, int narea, int ntime){
	
	arma::mat XB(narea, ntime, arma::fill::zeros);
	for(int i = 0; i < narea; i++){
		for(int t = 0; t < ntime; t++){
			arma::rowvec X1 = X(arma::span(i), arma::span(t), arma::span::all);
			XB(i,t) = arma::as_scalar(X1 * betaa);	
		}
	}
	
	double post_value = arma::accu(Y % XB) - arma::accu(O % z_t % theta_t % exp(XB)) - 0.5 * (arma::as_scalar((betaa - mu_beta).t() * arma::inv(sigma_beta) * (betaa - mu_beta)));

	return post_value;
}
	
// Logposterior distribution of delta
double logpost_delta(arma::cube V, arma::colvec deltaa, arma::mat z, arma::colvec mu_delta, arma::mat sigma_delta, int narea, int nseason){
	
	arma::mat VD(narea, nseason, arma::fill::zeros);
	for(int i = 0; i < narea; i++){
		for(int s = 0; s < nseason; s++){
			arma::rowvec V1 = V(arma::span(i), arma::span(s), arma::span::all);
			VD(i,s) = arma::as_scalar(V1 * deltaa);
		}
	}
  
	double post_value = 0.5 * (arma::accu(VD) - arma::accu(exp(VD) % pow(z - 1, 2)/z) - arma::as_scalar((deltaa - mu_delta).t() * arma::inv(sigma_delta) * (deltaa - mu_delta)));
  
	return post_value;
}

// Logposterior distribution of cc
double logpost_cc(double upsilon, double kappaa, arma::colvec cc, int cc0, arma::colvec uu, double zeta, arma::colvec rho, double w, int qq, int ss){
  
	cc(ss) = cc0;
	int c_size = cc.size();
	
	double sum_aux = 0.0;
	for(int h = 0; h <= qq; h++){
		
		int i_start = ss + h - qq; if(i_start < 0){ i_start = 0; }
		int i_end = ss + h; if(i_end >= (c_size - 1)){ i_end = c_size - 1; }

		sum_aux += lgamma(upsilon + kappaa + sum(cc.rows(i_start, i_end))) - lgamma(kappaa + sum(cc.rows(i_start, i_end) - uu.rows(i_start, i_end)));
	}

	int r_end = ss + qq; if(r_end > (c_size - 1)){ r_end = c_size - 1; }
	double post_value = sum_aux - internal::lfactorial(cc(ss) - uu(ss)) + cc(ss) * (log(zeta) + log(1 - w) + sum(log(1 - rho.rows(ss, r_end))));

	return post_value; 
}

// Logposterior distribution of uu
double logpost_uu(double upsilon, double kappaa, arma::colvec cc, arma::colvec uu, int uu0, arma::colvec rho, double w, int qq, int ss){
  
	uu(ss) = uu0;
	int u_size = uu.size();
	
	double sum_aux = 0.0;
	for(int h = 0; h <= qq; h++){

		int i_start = ss + h - qq; if(i_start < 0){ i_start = 0; }
		int i_end = ss + h; if(i_end >= (u_size - 1)){ i_end = u_size - 1; }
	
		sum_aux += lgamma(upsilon + sum(uu.rows(i_start, i_end))) + lgamma(kappaa + sum(cc.rows(i_start, i_end) - uu.rows(i_start, i_end)));
	}
  
	int r_end = ss + qq; if(r_end > (u_size - 1)){ r_end = u_size - 1; }
	double post_value = - sum_aux - internal::lfactorial(uu(ss)) - internal::lfactorial(cc(ss) - uu(ss)) + uu(ss) * ( log(w/(1 - w)) + sum( log(rho.rows(ss, r_end)/(1 - rho.rows(ss, r_end))) ) );
  
	return post_value; 
}

// Logposterior distribution of upsilon
double logpost_upsilon(double upsilon, double kappaa, arma::colvec cc, arma::colvec uu, arma::colvec rho, double w, double shape_upsilon, double rate_upsilon, int qq, int nseason){
  
	double sum_aux = 0.0;
	for(int s = 0; s < nseason; s++){
		
		int i_start = s - qq; if(i_start < 0){ i_start = 0; }
		sum_aux += lgamma(upsilon + kappaa + sum(cc.rows(i_start, s))) - lgamma(upsilon + sum(uu.rows(i_start, s)));
	}
	
	double post_value = sum_aux + upsilon * (sum(log(rho)) + log(w) - rate_upsilon) + lgamma(upsilon + kappaa) - lgamma(upsilon) + (shape_upsilon - 1) * log(upsilon);
  
	return post_value; 
}

// Logposterior distribution of kappa
double logpost_kappa(double upsilon, double kappaa, arma::colvec cc, arma::colvec uu, arma::colvec rho, double w, double shape_kappa, double rate_kappa, int qq, int nseason){
  
	double sum_aux = 0.0;
	for(int s = 0; s < nseason; s++){
    
		int i_start = s - qq; if(i_start < 0){ i_start = 0; }
		sum_aux += + lgamma(upsilon + kappaa + sum(cc.rows(i_start, s))) - lgamma(kappaa + sum(cc.rows(i_start, s) - uu.rows(i_start, s)));
	}
  
	double post_value = sum_aux + kappaa * ( sum(log(1 - rho)) + log(1 - w) - rate_kappa ) + lgamma(upsilon + kappaa) - lgamma(kappaa) + (shape_kappa - 1) * log(kappaa);
  
	return post_value; 
}

// GIG for z
double C_rgig(double p1, double p2, double p3){
	
    Function f("rgig");   
    SEXP res = f(1, Named("lambda")=p1, _["chi"]=p2, _["psi"]= p3);
	
	return *REAL(res);
}

// sample
int C_sample(arma::colvec x){
	
    Function f("sample");   
    SEXP res = f(Named("x")=x, _["size"]=1);
	
	return *REAL(res);
}

// Partition update
List C_partition_up(arma::mat adj_mat, int narea, int ncl_t, arma::colvec cl_lab, arma::colvec y_t, arma::colvec OXB_t, arma::colvec z_t, 
					double shape_theta, double rate_theta, double upsilon, double kappaa, arma::colvec uu, arma::colvec cc, int tt, int qq){
	
    Function f("partition_up");   

    List L = f(Named("adj_mat")=adj_mat, _["narea"]=narea, _["ncl_t"]=ncl_t, _["cl_lab"]=cl_lab, _["y_t"]=y_t, _["OXB_t"]=OXB_t, _["z_t"]=z_t,
			 _["shape_theta"]=shape_theta, _["rate_theta"]=rate_theta, _["upsilon"]=upsilon, _["kappaa"]=kappaa, _["uu"]=uu, _["cc"]=cc, _["tt"]=tt, _["qq"]=qq);
    
    return L;
}

// Partition update
List C_partition_up2(arma::mat adj_mat, int narea, int ncl_t, arma::colvec cl_lab, double rho_s){
	
    Function f("partition_up2");   

    List L = f(Named("adj_mat")=adj_mat, _["narea"]=narea, _["ncl_t"]=ncl_t, _["cl_lab"]=cl_lab, _["rho_s"]=rho_s);
    
    return L;
}

// Main function - PIG
// [[Rcpp::export]]
List STPPM_PIG(arma::mat Y, arma::mat O, arma::cube X, arma::cube V, arma::colvec season, arma::mat adj_mat, int qq, double shape_theta, double rate_theta, 
			   double shape_upsilon, double rate_upsilon, double shape_kappa, double rate_kappa, double shape_zeta, double rate_zeta,
			   arma::colvec mu_beta, arma::mat sigma_beta, arma::colvec mu_delta, arma::mat sigma_delta, 
			   int niter, int nburn, int nthin){

	RNGScope scope;
	
	// *************** 
	// Sample size
	// *************** 
	
	arma::colvec season_id = unique(season);
	int narea = Y.n_rows, ntime = Y.n_cols, nseason = season_id.n_rows, nbeta = X.n_slices, ndelta = V.n_slices;
	
	// *************** 
	// Objects to save samples
	// *************** 
	
	arma::uvec seq_thin = arma::regspace<arma::uvec>(nburn, nthin, niter); // Generate a seq to save results according to nthin - (start, delta, end)
	
	int nsample = ((niter - nburn)/nthin), jj = 0; 
	double acceptance_beta = 0.0, acceptance_delta = 0.0, acceptance_uu = 0.0, acceptance_cc = 0.0, acceptance_upsilon = 0.0, acceptance_kappa = 0.0;
	arma::colvec vec_w_post(nsample), vec_zeta_post(nsample), vec_upsilon_post(nsample), vec_kappa_post(nsample);
	arma::mat mat_ncl_post(nsample, nseason), mat_beta_post(nsample, nbeta), mat_delta_post(nsample, ndelta), 
			  mat_rho_post(nseason, nsample), mat_uu_post(nseason, nsample), mat_cc_post(nseason, nsample);
	arma::cube cube_partition_post(narea, nseason, nsample), cube_z_post(narea, nseason, nsample), cube_theta_post(narea, nseason, nsample);
	
	// *************** 
	// Initialize partition (with only one cluster)
	// *************** 
  
	arma::rowvec ncl(nseason + qq, arma::fill::ones);
	arma::mat cl_id(narea, nseason + qq, arma::fill::ones);
  
	// *************** 
	// Initialize parameters according to their prior distribution
	// *************** 
	
	arma::colvec rho(nseason + 2*qq, arma::fill::value(0.05)), cc(nseason + 2*qq, arma::fill::zeros), uu(nseason + 2*qq, arma::fill::zeros);
	double w = 0.05, zeta = 1.0, upsilon = 5.0, kappaa = 100.0;
  
	arma::mat theta(narea, nseason, arma::fill::ones), z(narea, nseason, arma::fill::ones), theta_t(narea, ntime, arma::fill::ones), z_t(narea, ntime, arma::fill::ones);
	arma::colvec betaa(nbeta, arma::fill::zeros), deltaa(ndelta, arma::fill::zeros);
	arma::mat i_sigma_beta = arma::inv(sigma_beta), i_sigma_delta = arma::inv(sigma_delta);
 
	// *************** Auxiliary structures
	arma::mat XB(narea, ntime, arma::fill::zeros), VD(narea, nseason, arma::fill::zeros), OXB_s(narea, nseason, arma::fill::zeros), Y_s(narea, nseason, arma::fill::zeros);
	for(int i = 0; i < narea; i++){
		for(int s = 0; s < nseason; s++){
			for(int t = 0; t < ntime; t++){
			
				arma::rowvec X1 = X(arma::span(i), arma::span(t), arma::span::all);
				XB(i,t) = arma::as_scalar(X1 * betaa);	
			
				if(season(t) == season_id(s)){
					Y_s(i,s) += Y(i,t);
					OXB_s(i,s) += O(i,t) * exp(XB(i,t));
				}
			}
			arma::rowvec V1 = V(arma::span(i), arma::span(s), arma::span::all);
			VD(i,s) = arma::as_scalar(V1 * deltaa);
		}
	}
	arma::mat psi(narea, nseason);
	psi = exp(VD);
	
	arma::mat sigma_beta_cand(nbeta, nbeta, arma::fill::eye), sigma_delta_cand(ndelta, ndelta, arma::fill::eye);
	sigma_beta_cand = sigma_beta_cand * 1e-5, sigma_delta_cand = sigma_delta_cand * 5e-3;
	double sd_ups = 2.0, sd_kap = 20.0;
	
	// ***************
	// Gibbs Sampling
	// ***************
	for(int iter = 1; iter < niter; iter++){
	
		double upsilon_candidate = r_truncnorm(upsilon, sd_ups, 0.0, datum::inf); 
		double accept_prob_upsilon = logpost_upsilon(upsilon_candidate, kappaa, cc.rows(0, nseason - 1), uu.rows(0, nseason - 1), rho.rows(0, nseason - 1), w, shape_upsilon, rate_upsilon, qq, nseason) 
									+ d_truncnorm(upsilon, upsilon_candidate, sd_ups, 0.0, datum::inf, 1)
									- logpost_upsilon(upsilon, kappaa, cc.rows(0, nseason - 1), uu.rows(0, nseason - 1), rho.rows(0, nseason - 1), w, shape_upsilon, rate_upsilon, qq, nseason)
									- d_truncnorm(upsilon_candidate, upsilon, sd_ups, 0.0, datum::inf, 1);				
									
		if((log(R::runif(0, 1)) <= std::min(0.0, accept_prob_upsilon))){
			
			upsilon = upsilon_candidate;
			acceptance_upsilon += 1;
		}

		double kappa_candidate = r_truncnorm(kappaa, sd_kap, 0.0, datum::inf);
		double accept_prob_kappa = logpost_kappa(upsilon, kappa_candidate, cc.rows(0, nseason - 1), uu.rows(0, nseason - 1), rho.rows(0, nseason - 1), w, shape_kappa, rate_kappa, qq, nseason)
								   + d_truncnorm(kappaa, kappa_candidate, sd_kap, 0, datum::inf, 1)
								   - logpost_kappa(upsilon, kappaa, cc.rows(0, nseason - 1), uu.rows(0, nseason - 1), rho.rows(0, nseason - 1), w, shape_kappa, rate_kappa, qq, nseason)
								   - d_truncnorm(kappa_candidate, kappaa, sd_kap, 0, datum::inf, 1);			 
    
		if((log(R::runif(0, 1)) <= std::min(0.0, accept_prob_kappa))){
			
			kappaa = kappa_candidate;
			acceptance_kappa += 1;
		}

		double shape1_w_post = upsilon + sum(uu.rows(0, nseason - 1));
		double shape2_w_post = kappaa + sum(cc.rows(0, nseason - 1) - uu.rows(0, nseason - 1));
		w = Rf_rbeta(shape1_w_post, shape2_w_post);

		double shape_zeta_post = shape_zeta + sum(cc.rows(0, nseason - 1));
		double scale_zeta_post = 1/(rate_zeta + nseason);
		zeta = Rf_rgamma(shape_zeta_post, scale_zeta_post);

		for(int ss = 0; ss < nseason; ss++){
	
			arma::uvec ts = arma::find(season == season_id(ss));

			// Update partition
			List new_part = C_partition_up(adj_mat, narea, ncl(ss), cl_id.col(ss), Y_s.col(ss), OXB_s.col(ss), z.col(ss), shape_theta, rate_theta, upsilon, kappaa, uu.rows(0, nseason - 1), cc.rows(0, nseason - 1), ss, qq);
			arma::colvec cl_aux = new_part["cl_id"];
			int ncl_aux = new_part["ncl"];
			cl_id.col(ss) = cl_aux; 
			ncl(ss) = ncl_aux;

			// Update cc
			arma::colvec opts_cc_aux(12, arma::fill::zeros); 
			int l = 0;
			for(int j = std::max(uu(ss) + 0.0, cc(ss) - 5.0); j <= (cc(ss) + 5.0); j++){
				opts_cc_aux(l) = j;
				l += 1;
			}
			arma::colvec opts_cc(l, arma::fill::zeros);
			opts_cc = opts_cc_aux.rows(0, (l - 1)); 

			int cc_candidate = C_sample(opts_cc); 
			double accept_prob_cc = logpost_cc(upsilon, kappaa, cc, cc_candidate, uu, zeta, rho, w, qq, ss) - logpost_cc(upsilon, kappaa, cc, cc(ss), uu, zeta, rho, w, qq, ss);

	  		if((log(R::runif(0, 1)) <= std::min(0.0, accept_prob_cc))){
			
				cc(ss) = cc_candidate;
				acceptance_cc += 1;
			}

			// Update uu
			arma::colvec opts_uu_aux(12, arma::fill::zeros);
			l = 0;
			for(int j = std::max(0.0, uu(ss) - 5.0); j <= std::min(cc(ss), uu(ss) + 5.0); j++){
				opts_uu_aux(l) = j;
				l += 1;
			} 
			arma::colvec opts_uu(l, arma::fill::zeros);
			opts_uu = opts_uu_aux.rows(0, (l - 1));
			
			int uu_candidate = C_sample(opts_uu);
			double accept_prob_uu = logpost_uu(upsilon, kappaa, cc, uu, uu_candidate, rho, w, qq, ss) - logpost_uu(upsilon, kappaa, cc, uu, uu(ss), rho, w, qq, ss);
   
	  		if((log(R::runif(0, 1)) <= std::min(0.0, accept_prob_uu))){
			
				uu(ss) = uu_candidate;
				acceptance_uu += 1;
			}

			// Update rho
			int i_start = ss - qq; if(i_start < 0){ i_start = 0; }
			double shape1_rho_post = ncl_aux + upsilon - 1 + sum(uu.rows(i_start, ss));
			double shape2_rho_post = narea - ncl_aux + kappaa + sum(cc.rows(i_start, ss) - uu.rows(i_start, ss));
			rho(ss) = Rf_rbeta(shape1_rho_post, shape2_rho_post);

			// Update cluster-specific parameters
			arma::colvec Ys_aux = Y_s.col(ss), z_aux = z.col(ss), OXBs_aux = OXB_s.col(ss), theta_aux(narea, arma::fill::zeros);
			for(int j = 1; j < (ncl_aux + 1); j++){
        
				arma::uvec Sj = arma::find(cl_id.col(ss) == j);
				int size_cl = Sj.size();
				
				double Y_sj = arma::accu( Ys_aux.elem( Sj ) );
				double zOXB_sj = arma::accu( z_aux.elem( Sj ) % OXBs_aux.elem( Sj ));
	
				double shape_theta_post = Y_sj + shape_theta;
				double scale_theta_post = 1/(zOXB_sj + rate_theta);
        
				double theta_star = Rf_rgamma(shape_theta_post, scale_theta_post);
				
				arma::colvec theta_aux1( size_cl, arma::fill::value(theta_star) );
				theta_aux.elem(Sj) = theta_aux1;
			}
			theta.col(ss) = theta_aux;
			theta_t.each_col(ts) = theta_aux;

			// Update z
			arma::colvec z_aux1(narea, arma::fill::zeros);
			for(int i = 0; i < narea; i++){
        
				double par1 = Y_s(i,ss) - 0.5;
				double par2 = 2 * theta(i,ss) * OXB_s(i,ss) + psi(i,ss);

				z(i,ss) = C_rgig(par1, psi(i,ss), par2);
				z_aux1(i) = z(i,ss);
			}
			z_t.each_col(ts) = z_aux1; 
		}
		for(int ss = nseason; ss < (nseason + qq); ss++){
	 
			// Update partition
			List new_part = C_partition_up2(adj_mat, narea, ncl(ss), cl_id.col(ss), rho(ss));
			arma::colvec cl_aux = new_part["cl_id"];
			int ncl_aux = new_part["ncl"];
			cl_id(arma::span::all, arma::span(ss)) = cl_aux; 
			ncl(ss) = ncl_aux;

			// Update cc
			arma::colvec opts_cc_aux(12, arma::fill::zeros); 
			int l = 0;
			for(int j = std::max(uu(ss) + 0.0, cc(ss) - 5.0); j <= (cc(ss) + 5.0); j++){
				opts_cc_aux(l) = j;
				l += 1;
			}
			arma::colvec opts_cc(l, arma::fill::zeros);
			opts_cc = opts_cc_aux.rows(0, (l - 1));
			
			int cc_candidate = C_sample(opts_cc);
			double accept_prob_cc = logpost_cc(upsilon, kappaa, cc, cc_candidate, uu, zeta, rho, w, qq, ss) - logpost_cc(upsilon, kappaa, cc, cc(ss), uu, zeta, rho, w, qq, ss);

	  		if((log(R::runif(0, 1)) <= std::min(0.0, accept_prob_cc))){
			
				cc(ss) = cc_candidate;
				acceptance_cc += 1;
			}

			// Update uu
			arma::colvec opts_uu_aux(12, arma::fill::zeros);
			l = 0;
			for(int j = std::max(0.0, uu(ss) - 5.0); j <= std::min(cc(ss), uu(ss) + 5.0); j++){
				opts_uu_aux(l) = j;
				l += 1;
			} 
			arma::colvec opts_uu(l, arma::fill::zeros);
			opts_uu = opts_uu_aux.rows(0, (l - 1));
			
			int uu_candidate = C_sample(opts_uu);
			double accept_prob_uu = logpost_uu(upsilon, kappaa, cc, uu, uu_candidate, rho, w, qq, ss) - logpost_uu(upsilon, kappaa, cc, uu, uu(ss), rho, w, qq, ss);
      
	  		if((log(R::runif(0, 1)) <= std::min(0.0, accept_prob_uu))){
			
				uu(ss) = uu_candidate;
				acceptance_uu += 1;
			}
		
			// Update rho
			double shape1_rho_post = ncl_aux + upsilon - 1 + sum(uu.rows(ss - qq, ss));
			double shape2_rho_post = narea - ncl_aux + kappaa + sum(cc.rows(ss - qq, ss) - uu.rows(ss - qq, ss));
			rho(ss) = Rf_rbeta(shape1_rho_post, shape2_rho_post);
		}	
 
		arma::colvec beta_candidate = rmvnorm(1, betaa, sigma_beta_cand).t(); 
		double accept_prob_beta = logpost_beta(Y, O, X, theta_t, z_t, beta_candidate, mu_beta, sigma_beta, narea, ntime) - logpost_beta(Y, O, X, theta_t, z_t, betaa, mu_beta, sigma_beta, narea, ntime);
    
		if((log(R::runif(0, 1)) <= std::min(0.0, accept_prob_beta))){
			
			betaa = beta_candidate;
			acceptance_beta += 1;
			
			OXB_s.zeros();
			for(int i = 0; i < narea; i++){
				
				int s = 0, t = 0;
				while(s < nseason){
					while(t < ntime && season(t) == season_id(s)){ 
						
						arma::rowvec X1 = X(arma::span(i), arma::span(t), arma::span::all);
						XB(i,t) = arma::as_scalar(X1 * betaa);	
						OXB_s(i,s) += O(i,t) * exp(XB(i,t));
						
						t += 1;
					}
					s += 1;	
				}
			}
		}

		arma::colvec delta_candidate = rmvnorm(1, deltaa, sigma_delta_cand).t(); 
		double accept_prob_delta = logpost_delta(V, delta_candidate, z, mu_delta, sigma_delta, narea, nseason) - logpost_delta(V, deltaa, z, mu_delta, sigma_delta, narea, nseason);					
		
		if((log(R::runif(0, 1)) <= std::min(0.0, accept_prob_delta))){
			
			deltaa = delta_candidate;
			acceptance_delta += 1;
			
			for(int i = 0; i < narea; i++){
				for(int s = 0; s < nseason; s++){
					arma::rowvec V1 = V(arma::span(i), arma::span(s), arma::span::all);
					VD(i,s) = arma::as_scalar(V1 * deltaa);
				}
			}
			psi = exp(VD);
		}
 
		// ******************************************************************
		// Save sample
		if(iter >= nburn){ 
			if(belong(iter, seq_thin)){
				
				// Rprintf("iter = %i \n", iter);
				
				mat_ncl_post.row(jj) = ncl.cols(0, nseason - 1);
				cube_partition_post(arma::span::all, arma::span::all, arma::span(jj)) = cl_id(arma::span::all, arma::span(0, nseason - 1));

				cube_z_post(arma::span::all, arma::span::all, arma::span(jj)) = z(arma::span::all, arma::span::all);
				cube_theta_post(arma::span::all, arma::span::all, arma::span(jj)) = theta(arma::span::all, arma::span::all);

				mat_beta_post(arma::span(jj), arma::span::all) = betaa.t();
				mat_delta_post(arma::span(jj), arma::span::all) = deltaa.t();

				mat_rho_post(arma::span::all, arma::span(jj)) = rho.rows(0, nseason - 1);	
				mat_uu_post(arma::span::all, arma::span(jj)) = uu.rows(0, nseason - 1);
				mat_cc_post(arma::span::all, arma::span(jj)) = cc.rows(0, nseason - 1);
			
				vec_w_post(jj) = w;
				vec_zeta_post(jj) = zeta;
				vec_upsilon_post(jj) = upsilon; 
				vec_kappa_post(jj) = kappaa;
				
				jj += 1;
			}
		}	
	} 
	
	List L = List::create( _["ncl"] = mat_ncl_post, _["partition"] = cube_partition_post, _["z"] = cube_z_post, _["theta"] = cube_theta_post, 
						   _["beta"] = mat_beta_post, _["acceptance_beta"] = acceptance_beta/niter, _["delta"] = mat_delta_post, _["acceptance_delta"] = acceptance_delta/niter,
						   _["rho"] = mat_rho_post.t(), _["u"] = mat_uu_post.t(), _["acceptance_u"] = acceptance_uu/(niter*nseason), 
						   _["c"] = mat_cc_post.t(), _["acceptance_c"] = acceptance_cc/(niter*nseason), _["w"] = vec_w_post, _["zeta"] = vec_zeta_post,
						   _["upsilon"] = vec_upsilon_post, _["acceptance_upsilon"] = acceptance_upsilon/niter,
						   _["kappa"] = vec_kappa_post, _["acceptance_kappa"] = acceptance_kappa/niter );
	
	return L;	
}

// Main function - Poisson
// [[Rcpp::export]]
List STPPM_Poi(arma::mat Y, arma::mat O, arma::cube X, arma::colvec season, arma::mat adj_mat, int qq, double shape_theta, double rate_theta, 
			   double shape_upsilon, double rate_upsilon, double shape_kappa, double rate_kappa, double shape_zeta, double rate_zeta,
			   arma::colvec mu_beta, arma::mat sigma_beta, int niter, int nburn, int nthin){
   
	RNGScope scope;
	
	// *************** 
	// Sample size
	// *************** 
	
	arma::colvec season_id = unique(season);
	int narea = Y.n_rows, ntime = Y.n_cols, nseason = season_id.n_rows, nbeta = X.n_slices;
	
	// *************** 
	// Objects to save samples
	// *************** 
	
	arma::uvec seq_thin = arma::regspace<arma::uvec>(nburn, nthin, niter); // Generate a seq to save results according to nthin - (start, delta, end)
	
	int nsample = ((niter - nburn)/nthin), jj = 0; 
	double acceptance_beta = 0.0, acceptance_uu = 0.0, acceptance_cc = 0.0, acceptance_upsilon = 0.0, acceptance_kappa = 0.0;
	arma::colvec vec_w_post(nsample), vec_zeta_post(nsample), vec_upsilon_post(nsample), vec_kappa_post(nsample);
	arma::mat mat_ncl_post(nsample, nseason), mat_beta_post(nsample, nbeta), mat_rho_post(nseason, nsample), mat_uu_post(nseason, nsample), mat_cc_post(nseason, nsample);
	arma::cube cube_partition_post(narea, nseason, nsample), cube_theta_post(narea, nseason, nsample);
	
	// *************** 
	// Initialize partition (with only one cluster)
	// *************** 
  
	arma::rowvec ncl(nseason + qq, arma::fill::ones);
	arma::mat cl_id(narea, nseason + qq, arma::fill::ones);
  
	// *************** 
	// Initialize parameters according to their prior distribution
	// *************** 
	
	arma::colvec rho(nseason + 2*qq, arma::fill::value(0.05)), cc(nseason + 2*qq, arma::fill::zeros), uu(nseason + 2*qq, arma::fill::zeros);
	double w = 0.05, zeta = 1.0, upsilon = 5.0, kappaa = 100.0;
  
	arma::mat theta(narea, nseason, arma::fill::ones), z(narea, nseason, arma::fill::ones), theta_t(narea, ntime, arma::fill::ones), z_t(narea, ntime, arma::fill::ones);
	arma::colvec betaa(nbeta, arma::fill::zeros);
	arma::mat i_sigma_beta = arma::inv(sigma_beta);
 
	// *************** Auxiliary structures
	arma::mat XB(narea, ntime, arma::fill::zeros), OXB_s(narea, nseason, arma::fill::zeros), Y_s(narea, nseason, arma::fill::zeros);
	for(int i = 0; i < narea; i++){
		for(int s = 0; s < nseason; s++){
			for(int t = 0; t < ntime; t++){
			
				arma::rowvec X1 = X(arma::span(i), arma::span(t), arma::span::all);
				XB(i,t) = arma::as_scalar(X1 * betaa);	
			
				if(season(t) == season_id(s)){
					Y_s(i,s) += Y(i,t);
					OXB_s(i,s) += O(i,t) * exp(XB(i,t));
				}
			}
		}
	}
	
	arma::mat sigma_beta_cand(nbeta, nbeta, arma::fill::eye);
	sigma_beta_cand = sigma_beta_cand * 1e-5;
	double sd_ups = 2.0, sd_kap = 20.0;
	
	// ***************
	// Gibbs Sampling
	// ***************
	for(int iter = 1; iter < niter; iter++){
		
		double upsilon_candidate = r_truncnorm(upsilon, sd_ups, 0.0, datum::inf); 
		double accept_prob_upsilon = logpost_upsilon(upsilon_candidate, kappaa, cc.rows(0, nseason - 1), uu.rows(0, nseason - 1), rho.rows(0, nseason - 1), w, shape_upsilon, rate_upsilon, qq, nseason) 
									+ d_truncnorm(upsilon, upsilon_candidate, sd_ups, 0.0, datum::inf, 1)
									- logpost_upsilon(upsilon, kappaa, cc.rows(0, nseason - 1), uu.rows(0, nseason - 1), rho.rows(0, nseason - 1), w, shape_upsilon, rate_upsilon, qq, nseason)
									- d_truncnorm(upsilon_candidate, upsilon, sd_ups, 0.0, datum::inf, 1);				
									
		if((log(R::runif(0, 1)) <= std::min(0.0, accept_prob_upsilon))){
			
			upsilon = upsilon_candidate;
			acceptance_upsilon += 1;
		}

		double kappa_candidate = r_truncnorm(kappaa, sd_kap, 0.0, datum::inf);
		double accept_prob_kappa = logpost_kappa(upsilon, kappa_candidate, cc.rows(0, nseason - 1), uu.rows(0, nseason - 1), rho.rows(0, nseason - 1), w, shape_kappa, rate_kappa, qq, nseason)
								   + d_truncnorm(kappaa, kappa_candidate, sd_kap, 0, datum::inf, 1)
								   - logpost_kappa(upsilon, kappaa, cc.rows(0, nseason - 1), uu.rows(0, nseason - 1), rho.rows(0, nseason - 1), w, shape_kappa, rate_kappa, qq, nseason)
								   - d_truncnorm(kappa_candidate, kappaa, sd_kap, 0, datum::inf, 1);			 
    
		if((log(R::runif(0, 1)) <= std::min(0.0, accept_prob_kappa))){
			
			kappaa = kappa_candidate;
			acceptance_kappa += 1;
		}

		double shape1_w_post = upsilon + sum(uu.rows(0, nseason - 1));
		double shape2_w_post = kappaa + sum(cc.rows(0, nseason - 1) - uu.rows(0, nseason - 1));
		w = Rf_rbeta(shape1_w_post, shape2_w_post);

		double shape_zeta_post = shape_zeta + sum(cc.rows(0, nseason - 1));
		double scale_zeta_post = 1/(rate_zeta + nseason);
		zeta = Rf_rgamma(shape_zeta_post, scale_zeta_post);

		for(int ss = 0; ss < nseason; ss++){
      
			arma::uvec ts = arma::find(season == season_id(ss));
			
			// Update partition
			List new_part = C_partition_up(adj_mat, narea, ncl(ss), cl_id.col(ss), Y_s.col(ss), OXB_s.col(ss), z.col(ss), shape_theta, rate_theta, upsilon, kappaa, uu.rows(0, nseason - 1), cc.rows(0, nseason - 1), ss, qq);
			arma::colvec cl_aux = new_part["cl_id"];
			int ncl_aux = new_part["ncl"];
			cl_id.col(ss) = cl_aux; 
			ncl(ss) = ncl_aux;

			// Update cc
			arma::colvec opts_cc_aux(12, arma::fill::zeros); 
			int l = 0;
			for(int j = std::max(uu(ss) + 0.0, cc(ss) - 5.0); j <= (cc(ss) + 5.0); j++){
				opts_cc_aux(l) = j;
				l += 1;
			}
			arma::colvec opts_cc(l, arma::fill::zeros);
			opts_cc = opts_cc_aux.rows(0, (l - 1)); 

			int cc_candidate = C_sample(opts_cc);
			double accept_prob_cc = logpost_cc(upsilon, kappaa, cc, cc_candidate, uu, zeta, rho, w, qq, ss) - logpost_cc(upsilon, kappaa, cc, cc(ss), uu, zeta, rho, w, qq, ss);

	  		if((log(R::runif(0, 1)) <= std::min(0.0, accept_prob_cc))){
			
				cc(ss) = cc_candidate;
				acceptance_cc += 1;
			}
	
			// Update uu
			arma::colvec opts_uu_aux(12, arma::fill::zeros);
			l = 0;
			for(int j = std::max(0.0, uu(ss) - 5.0); j <= std::min(cc(ss), uu(ss) + 5.0); j++){
				opts_uu_aux(l) = j;
				l += 1;
			} 
			arma::colvec opts_uu(l, arma::fill::zeros);
			opts_uu = opts_uu_aux.rows(0, (l - 1));
			
			int uu_candidate = C_sample(opts_uu);
			double accept_prob_uu = logpost_uu(upsilon, kappaa, cc, uu, uu_candidate, rho, w, qq, ss) - logpost_uu(upsilon, kappaa, cc, uu, uu(ss), rho, w, qq, ss);
      
	  		if((log(R::runif(0, 1)) <= std::min(0.0, accept_prob_uu))){
			
				uu(ss) = uu_candidate;
				acceptance_uu += 1;
			}

			// Update rho
			int i_start = ss - qq; if(i_start < 0){ i_start = 0; }
			double shape1_rho_post = ncl_aux + upsilon - 1 + sum(uu.rows(i_start, ss));
			double shape2_rho_post = narea - ncl_aux + kappaa + sum(cc.rows(i_start, ss) - uu.rows(i_start, ss));
			rho(ss) = Rf_rbeta(shape1_rho_post, shape2_rho_post);

			// Update cluster-specific parameters
			arma::colvec Ys_aux = Y_s.col(ss), z_aux = z.col(ss), OXBs_aux = OXB_s.col(ss), theta_aux(narea, arma::fill::zeros);
			for(int j = 1; j < (ncl_aux + 1); j++){
        
				arma::uvec Sj = arma::find(cl_id.col(ss) == j);
				int size_cl = Sj.size();
				
				double Y_sj = arma::accu( Ys_aux.elem( Sj ) );
				double zOXB_sj = arma::accu( z_aux.elem( Sj ) % OXBs_aux.elem( Sj ));
	
				double shape_theta_post = Y_sj + shape_theta;
				double scale_theta_post = 1/(zOXB_sj + rate_theta);
        
				double theta_star = Rf_rgamma(shape_theta_post, scale_theta_post);
				
				arma::colvec theta_aux1( size_cl, arma::fill::value(theta_star) );
				theta_aux.elem(Sj) = theta_aux1;
			}
			theta.col(ss) = theta_aux;
			theta_t.each_col(ts) = theta_aux;
		}
		for(int ss = nseason; ss < (nseason + qq); ss++){
       
			// Update partition
			List new_part = C_partition_up2(adj_mat, narea, ncl(ss), cl_id.col(ss), rho(ss));
			arma::colvec cl_aux = new_part["cl_id"];
			int ncl_aux = new_part["ncl"];
			cl_id(arma::span::all, arma::span(ss)) = cl_aux; 
			ncl(ss) = ncl_aux;
 
			// Update cc
			arma::colvec opts_cc_aux(12, arma::fill::zeros); 
			int l = 0;
			for(int j = std::max(uu(ss) + 0.0, cc(ss) - 5.0); j <= (cc(ss) + 5.0); j++){
				opts_cc_aux(l) = j;
				l += 1;
			}
			arma::colvec opts_cc(l, arma::fill::zeros);
			opts_cc = opts_cc_aux.rows(0, (l - 1));
			
			int cc_candidate = C_sample(opts_cc);
			double accept_prob_cc = logpost_cc(upsilon, kappaa, cc, cc_candidate, uu, zeta, rho, w, qq, ss) - logpost_cc(upsilon, kappaa, cc, cc(ss), uu, zeta, rho, w, qq, ss);

	  		if((log(R::runif(0, 1)) <= std::min(0.0, accept_prob_cc))){
			
				cc(ss) = cc_candidate;
				acceptance_cc += 1;
			}
	
			// Update uu
			arma::colvec opts_uu_aux(12, arma::fill::zeros);
			l = 0;
			for(int j = std::max(0.0, uu(ss) - 5.0); j <= std::min(cc(ss), uu(ss) + 5.0); j++){
				opts_uu_aux(l) = j;
				l += 1;
			} 
			arma::colvec opts_uu(l, arma::fill::zeros);
			opts_uu = opts_uu_aux.rows(0, (l - 1));
			
			int uu_candidate = C_sample(opts_uu);
			double accept_prob_uu = logpost_uu(upsilon, kappaa, cc, uu, uu_candidate, rho, w, qq, ss) - logpost_uu(upsilon, kappaa, cc, uu, uu(ss), rho, w, qq, ss);
      
	  		if((log(R::runif(0, 1)) <= std::min(0.0, accept_prob_uu))){
			
				uu(ss) = uu_candidate;
				acceptance_uu += 1;
			}
     
			// Update rho
			double shape1_rho_post = ncl(ss) + upsilon - 1 + sum(uu.rows(ss - qq, ss));
			double shape2_rho_post = narea - ncl(ss) + kappaa + sum(cc.rows(ss - qq, ss) - uu.rows(ss - qq, ss));
			rho(ss) = Rf_rbeta(shape1_rho_post, shape2_rho_post);
		}	
   
		arma::colvec beta_candidate = rmvnorm(1, betaa, sigma_beta_cand).t(); 
		double accept_prob_beta = logpost_beta(Y, O, X, theta_t, z_t, beta_candidate, mu_beta, sigma_beta, narea, ntime) - logpost_beta(Y, O, X, theta_t, z_t, betaa, mu_beta, sigma_beta, narea, ntime);
    
		if((log(R::runif(0, 1)) <= std::min(0.0, accept_prob_beta))){
			
			betaa = beta_candidate;
			acceptance_beta += 1;
			
			OXB_s.zeros();
			for(int i = 0; i < narea; i++){
				
				int s = 0, t = 0;
				while(s < nseason){
					while(t < ntime && season(t) == season_id(s)){ 
						
						arma::rowvec X1 = X(arma::span(i), arma::span(t), arma::span::all);
						XB(i,t) = arma::as_scalar(X1 * betaa);	
						OXB_s(i,s) += O(i,t) * exp(XB(i,t));
						
						t += 1;
					}
					s += 1;	
				}
			}
		}
  
		// ******************************************************************
		// Save sample
		if(iter >= nburn){ 
			if(belong(iter, seq_thin)){
				
				// Rprintf("iter = %i \n", iter);
				
				mat_ncl_post.row(jj) = ncl.cols(0, nseason - 1);
				cube_partition_post(arma::span::all, arma::span::all, arma::span(jj)) = cl_id(arma::span::all, arma::span(0, nseason - 1));

				cube_theta_post(arma::span::all, arma::span::all, arma::span(jj)) = theta(arma::span::all, arma::span::all);
				mat_beta_post(arma::span(jj), arma::span::all) = betaa.t();

				mat_rho_post(arma::span::all, arma::span(jj)) = rho.rows(0, nseason - 1);	
				mat_uu_post(arma::span::all, arma::span(jj)) = uu.rows(0, nseason - 1);
				mat_cc_post(arma::span::all, arma::span(jj)) = cc.rows(0, nseason - 1);
			
				vec_w_post(jj) = w;
				vec_zeta_post(jj) = zeta;
				vec_upsilon_post(jj) = upsilon; 
				vec_kappa_post(jj) = kappaa;
				
				jj += 1;
			}
		}	
	} 
	
	List L = List::create( _["ncl"] = mat_ncl_post, _["partition"] = cube_partition_post, _["theta"] = cube_theta_post, 
						   _["beta"] = mat_beta_post, _["acceptance_beta"] = acceptance_beta/niter,
						   _["rho"] = mat_rho_post.t(), _["u"] = mat_uu_post.t(), _["acceptance_u"] = acceptance_uu/(niter*nseason), 
						   _["c"] = mat_cc_post.t(), _["acceptance_c"] = acceptance_cc/(niter*nseason), _["w"] = vec_w_post, _["zeta"] = vec_zeta_post,
						   _["upsilon"] = vec_upsilon_post, _["acceptance_upsilon"] = acceptance_upsilon/niter,
						   _["kappa"] = vec_kappa_post, _["acceptance_kappa"] = acceptance_kappa/niter );
	
	return L;	
}

// Logposterior distribution of upsilon (iid partitions)
double logpost_upsilon_iid(double upsilon, double kappaa, arma::colvec rho, double shape_upsilon, double rate_upsilon, int nseason){
  
	double post_value = nseason * (lgamma(upsilon + kappaa) - lgamma(upsilon)) + upsilon * (sum(log(rho)) - rate_upsilon) + (shape_upsilon - 1) * log(upsilon);
  
	return post_value; 
}

// Logposterior distribution of kappa (iid partitions)
double logpost_kappa_iid(double upsilon, double kappaa, arma::colvec rho, double shape_kappa, double rate_kappa, int nseason){
  
	double post_value = nseason * (lgamma(upsilon + kappaa) - lgamma(kappaa)) + kappaa * (sum(log(1 - rho)) - rate_kappa) + (shape_kappa - 1) * log(kappaa);
  
	return post_value; 
}

// Main function - PIG with iid partitions
// [[Rcpp::export]]
List STPPM_PIG_iid(arma::mat Y, arma::mat O, arma::cube X, arma::cube V, arma::colvec season, arma::mat adj_mat, double shape_theta, double rate_theta, 
				   double shape_upsilon, double rate_upsilon, double shape_kappa, double rate_kappa, arma::colvec mu_beta, arma::mat sigma_beta, 
				   arma::colvec mu_delta, arma::mat sigma_delta, int niter, int nburn, int nthin){

	RNGScope scope;
	
	// *************** 
	// Sample size
	// *************** 
	
	arma::colvec season_id = unique(season);
	int narea = Y.n_rows, ntime = Y.n_cols, nseason = season_id.n_rows, nbeta = X.n_slices, ndelta = V.n_slices;
	
	// *************** 
	// Objects to save samples
	// *************** 
	
	arma::uvec seq_thin = arma::regspace<arma::uvec>(nburn, nthin, niter); // Generate a seq to save results according to nthin - (start, delta, end)
	
	int nsample = ((niter - nburn)/nthin), jj = 0; 
	double acceptance_beta = 0.0, acceptance_delta = 0.0, acceptance_upsilon = 0.0, acceptance_kappa = 0.0;
	arma::colvec vec_upsilon_post(nsample), vec_kappa_post(nsample);
	arma::mat mat_ncl_post(nsample, nseason), mat_beta_post(nsample, nbeta), mat_delta_post(nsample, ndelta), mat_rho_post(nseason, nsample);
	arma::cube cube_partition_post(narea, nseason, nsample), cube_z_post(narea, nseason, nsample), cube_theta_post(narea, nseason, nsample);
	
	// *************** 
	// Initialize partition (with only one cluster)
	// *************** 
  
	arma::rowvec ncl(nseason, arma::fill::ones);
	arma::mat cl_id(narea, nseason, arma::fill::ones);
  
	// *************** 
	// Initialize parameters according to their prior distribution
	// *************** 
	
	arma::colvec rho(nseason, arma::fill::value(0.05)), cc(nseason, arma::fill::zeros), uu(nseason, arma::fill::zeros);
	double upsilon = 5.0, kappaa = 100.0;
  
	arma::mat theta(narea, nseason, arma::fill::ones), z(narea, nseason, arma::fill::ones), theta_t(narea, ntime, arma::fill::ones), z_t(narea, ntime, arma::fill::ones);
	arma::colvec betaa(nbeta, arma::fill::zeros), deltaa(ndelta, arma::fill::zeros);
	arma::mat i_sigma_beta = arma::inv(sigma_beta), i_sigma_delta = arma::inv(sigma_delta);
 
	// *************** Auxiliary structures
	arma::mat XB(narea, ntime, arma::fill::zeros), VD(narea, nseason, arma::fill::zeros), OXB_s(narea, nseason, arma::fill::zeros), Y_s(narea, nseason, arma::fill::zeros);
	for(int i = 0; i < narea; i++){
		for(int s = 0; s < nseason; s++){
			for(int t = 0; t < ntime; t++){
			
				arma::rowvec X1 = X(arma::span(i), arma::span(t), arma::span::all);
				XB(i,t) = arma::as_scalar(X1 * betaa);	
			
				if(season(t) == season_id(s)){
					Y_s(i,s) += Y(i,t);
					OXB_s(i,s) += O(i,t) * exp(XB(i,t));
				}
			}
			arma::rowvec V1 = V(arma::span(i), arma::span(s), arma::span::all);
			VD(i,s) = arma::as_scalar(V1 * deltaa);
		}
	}
	arma::mat psi(narea, nseason);
	psi = exp(VD);
	
	arma::mat sigma_beta_cand(nbeta, nbeta, arma::fill::eye), sigma_delta_cand(ndelta, ndelta, arma::fill::eye);
	sigma_beta_cand = sigma_beta_cand * 1e-5, sigma_delta_cand = sigma_delta_cand * 5e-3;
	double sd_ups = 2.0, sd_kap = 20.0;
	
	// ***************
	// Gibbs Sampling
	// ***************
	for(int iter = 1; iter < niter; iter++){
		
		double upsilon_candidate = r_truncnorm(upsilon, sd_ups, 0.0, datum::inf); 
		double accept_prob_upsilon = logpost_upsilon_iid(upsilon_candidate, kappaa, rho.rows(0, nseason - 1), shape_upsilon, rate_upsilon, nseason) 
									+ d_truncnorm(upsilon, upsilon_candidate, sd_ups, 0.0, datum::inf, 1)
									- logpost_upsilon_iid(upsilon, kappaa, rho.rows(0, nseason - 1), shape_upsilon, rate_upsilon, nseason)
									- d_truncnorm(upsilon_candidate, upsilon, sd_ups, 0.0, datum::inf, 1);				
									
		if((log(R::runif(0, 1)) <= std::min(0.0, accept_prob_upsilon))){
			
			upsilon = upsilon_candidate;
			acceptance_upsilon += 1;
		}

		double kappa_candidate = r_truncnorm(kappaa, sd_kap, 0.0, datum::inf);
		double accept_prob_kappa = logpost_kappa_iid(upsilon, kappa_candidate, rho.rows(0, nseason - 1), shape_kappa, rate_kappa, nseason)
								   + d_truncnorm(kappaa, kappa_candidate, sd_kap, 0, datum::inf, 1)
								   - logpost_kappa_iid(upsilon, kappaa, rho.rows(0, nseason - 1), shape_kappa, rate_kappa, nseason)
								   - d_truncnorm(kappa_candidate, kappaa, sd_kap, 0, datum::inf, 1);			 
    
		if((log(R::runif(0, 1)) <= std::min(0.0, accept_prob_kappa))){
			
			kappaa = kappa_candidate;
			acceptance_kappa += 1;
		}

		for(int ss = 0; ss < nseason; ss++){
    
			arma::uvec ts = arma::find(season == season_id(ss));
			
			// Update partition
			List new_part = C_partition_up(adj_mat, narea, ncl(ss), cl_id.col(ss), Y_s.col(ss), OXB_s.col(ss), z.col(ss), shape_theta, rate_theta, upsilon, kappaa, uu, cc, ss, 1);
			arma::colvec cl_aux = new_part["cl_id"];
			int ncl_aux = new_part["ncl"];
			cl_id.col(ss) = cl_aux; 
			ncl(ss) = ncl_aux;

			// Update rho
			double shape1_rho_post = ncl_aux + upsilon - 1;
			double shape2_rho_post = narea - ncl_aux + kappaa;
			rho(ss) = Rf_rbeta(shape1_rho_post, shape2_rho_post);

			// Update cluster-specific parameters
			arma::colvec Ys_aux = Y_s.col(ss), z_aux = z.col(ss), OXBs_aux = OXB_s.col(ss), theta_aux(narea, arma::fill::zeros);
			for(int j = 1; j < (ncl_aux + 1); j++){
        
				arma::uvec Sj = arma::find(cl_id.col(ss) == j);
				int size_cl = Sj.size();
				
				double Y_sj = arma::accu( Ys_aux.elem( Sj ) );
				double zOXB_sj = arma::accu( z_aux.elem( Sj ) % OXBs_aux.elem( Sj ));
	
				double shape_theta_post = Y_sj + shape_theta;
				double scale_theta_post = 1/(zOXB_sj + rate_theta);
        
				double theta_star = Rf_rgamma(shape_theta_post, scale_theta_post);
				
				arma::colvec theta_aux1( size_cl, arma::fill::value(theta_star) );
				theta_aux.elem(Sj) = theta_aux1;
			}
			theta.col(ss) = theta_aux;
			theta_t.each_col(ts) = theta_aux;

			// Update z
			arma::colvec z_aux1(narea, arma::fill::zeros);
			for(int i = 0; i < narea; i++){
        
				double par1 = Y_s(i,ss) - 0.5;
				double par2 = 2 * theta(i,ss) * OXB_s(i,ss) + psi(i,ss);

				z(i,ss) = C_rgig(par1, psi(i,ss), par2);
				z_aux1(i) = z(i,ss);
			}
			z_t.each_col(ts) = z_aux1; 
		}

		arma::colvec beta_candidate = rmvnorm(1, betaa, sigma_beta_cand).t(); 
		double accept_prob_beta = logpost_beta(Y, O, X, theta_t, z_t, beta_candidate, mu_beta, sigma_beta, narea, ntime) - logpost_beta(Y, O, X, theta_t, z_t, betaa, mu_beta, sigma_beta, narea, ntime);
    
		if((log(R::runif(0, 1)) <= std::min(0.0, accept_prob_beta))){
			
			betaa = beta_candidate;
			acceptance_beta += 1;
			
			OXB_s.zeros();
			for(int i = 0; i < narea; i++){
				
				int s = 0, t = 0;
				while(s < nseason){
					while(t < ntime && season(t) == season_id(s)){ 
						
						arma::rowvec X1 = X(arma::span(i), arma::span(t), arma::span::all);
						XB(i,t) = arma::as_scalar(X1 * betaa);	
						OXB_s(i,s) += O(i,t) * exp(XB(i,t));
						
						t += 1;
					}
					s += 1;	
				}
			}
		}

		arma::colvec delta_candidate = rmvnorm(1, deltaa, sigma_delta_cand).t(); 
		double accept_prob_delta = logpost_delta(V, delta_candidate, z, mu_delta, sigma_delta, narea, nseason) - logpost_delta(V, deltaa, z, mu_delta, sigma_delta, narea, nseason);					
		
		if((log(R::runif(0, 1)) <= std::min(0.0, accept_prob_delta))){
			
			deltaa = delta_candidate;
			acceptance_delta += 1;
			
			for(int i = 0; i < narea; i++){
				for(int s = 0; s < nseason; s++){
					arma::rowvec V1 = V(arma::span(i), arma::span(s), arma::span::all);
					VD(i,s) = arma::as_scalar(V1 * deltaa);
				}
			}
			psi = exp(VD);
		}
  
		// ******************************************************************
		// Save sample
		if(iter >= nburn){ 
			if(belong(iter, seq_thin)){
				
				// Rprintf("iter = %i \n", iter);
				
				mat_ncl_post.row(jj) = ncl.cols(0, nseason - 1);
				cube_partition_post(arma::span::all, arma::span::all, arma::span(jj)) = cl_id(arma::span::all, arma::span(0, nseason - 1));

				cube_z_post(arma::span::all, arma::span::all, arma::span(jj)) = z(arma::span::all, arma::span::all);
				cube_theta_post(arma::span::all, arma::span::all, arma::span(jj)) = theta(arma::span::all, arma::span::all);

				mat_beta_post(arma::span(jj), arma::span::all) = betaa.t();
				mat_delta_post(arma::span(jj), arma::span::all) = deltaa.t();

				mat_rho_post(arma::span::all, arma::span(jj)) = rho.rows(0, nseason - 1);	
			
				vec_upsilon_post(jj) = upsilon; 
				vec_kappa_post(jj) = kappaa;
				
				jj += 1;
			}
		}	
	} 
	
	List L = List::create( _["ncl"] = mat_ncl_post, _["partition"] = cube_partition_post, _["z"] = cube_z_post, _["theta"] = cube_theta_post, 
						   _["beta"] = mat_beta_post, _["acceptance_beta"] = acceptance_beta/niter,
						   _["delta"] = mat_delta_post, _["acceptance_delta"] = acceptance_delta/niter,
						   _["rho"] = mat_rho_post.t(),_["upsilon"] = vec_upsilon_post, _["acceptance_upsilon"] = acceptance_upsilon/niter,
						   _["kappa"] = vec_kappa_post, _["acceptance_kappa"] = acceptance_kappa/niter );
	
	return L;	
}