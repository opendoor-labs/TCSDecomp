#include <RcppArmadillo.h>
#include <Rcpp.h>

// Correctly setup the build environment
// [[Rcpp::depends(RcppArmadillo)]]

// Protect against compilers without OpenMP
#ifdef _OPENMP
#include <omp.h>
#endif

// R's implementation of the Moore-Penrose pseudo matrix inverse
// [[Rcpp::export]]
arma::mat Rginv(const arma::mat& m){
  arma::mat U, V;
  arma::vec S;
  arma::svd(U, S, V, m, "dc");
  arma::uvec Positive = arma::find(S > 1E-06 * S(1));
  if(all(Positive)){
    arma::mat D = diagmat(S);
    return V * (1/D * U.t());
  }else if(!any(Positive)){
    return arma::zeros(m.n_rows, m.n_cols);
  }else{
    S.elem(Positive) = 1/S.elem(Positive);
    arma::mat D = diagmat(S);
    return V * D * U.t();
  }
}

// st = Sys.time()
// B_tt = B_tl = matrix(0, ncol = ncol(yt), nrow = nrow(sp$Ft))
// P_tt = P_tl = as.list(NA, ncol(yti))
// 
// #Initialize the filter
// B_LL = sp$B0
// P_LL = sp$P0
// 
// for(j in 1:ncol(yt)){
//   ################## Kalman filter routine ##################
//   B_tl[, j] = sp$Dt + sp$Ft %*% B_LL  #Initial estimate of unobserved values conditional on t-1
//   P_tl[[j]] = sp$Ft %*% P_LL %*% t(sp$Ft) + sp$Qt #Initial estimate of the covariance matrix conditional on t-1
//   N_tl = yt[, j] - sp$At - sp$Ht %*% B_tl[, j] #Prediction error conditoinal on t-1
//   F_tl = sp$Ht %*% P_tl[[j]] %*% t(sp$Ht) + sp$Rt #Variance of the predictoin error conditional on t-1
//   K_t = P_tl[[j]] %*% t(sp$Ht) %*% MASS::ginv(F_tl) #Kalman gain conditional on t-1
//   B_tt[, j] = B_tl[, j] + K_t %*% N_tl #Final estimate of the unobserved values
//   P_tt[[j]] = P_tl[[j]] - K_t %*% sp$Ht %*% P_tl[[j]] #Final estiamte of the covariance matrix
// 
//   #Reinitialize for the next iteration
//   B_LL = B_tt[, j]
//   P_LL = P_tt[[j]]
// }
// Sys.time() - st
// 
// #Kalman Smoother
// for(j in (ncol(yt) - 1):1){
//   B_tt[, j] = B_tt[, j] + P_tt[[j]] %*% t(sp$Ft) %*% ginv(P_tl[[j + 1]]) %*% (B_tt[, j + 1] - B_tl[, j + 1])
//   P_tt[[j]] = P_tt[[j]] + P_tt[[j]] %*% t(sp$Ft) %*% ginv(P_tl[[j + 1]]) %*% (P_tt[[j + 1]] - P_tl[[j + 1]]) %*% t(P_tt[[j]] %*% t(sp$Ft) %*% ginv(P_tl[[j + 1]]))
// }

// [[Rcpp::export]]
Rcpp::List kalman_filter(const arma::mat B0, const arma::mat P0, const arma::mat Dt, 
                         const arma::mat At, const arma::mat Ft, const arma::mat Ht, 
                         const arma::mat Qt, const arma::mat Rt, const arma::mat yt, 
                         const arma::mat X, const arma::mat beta){
  
  //Define the storage matrices
  arma::mat B_tt(Ft.n_rows, yt.n_cols);
  arma::mat B_tl(B_tt.n_rows, yt.n_cols);
  arma::cube P_tt(Ft.n_rows, Ft.n_rows, yt.n_cols);
  arma::cube P_tl(Ft.n_rows, Ft.n_rows, yt.n_cols);
  arma::mat N_t(yt.n_rows, yt.n_cols);
  arma::cube F_t(yt.n_rows, yt.n_rows, yt.n_cols);
  arma::cube K_t(Ft.n_rows, yt.n_rows, yt.n_cols);
  arma::uvec nonna_idx;
  arma::uvec na_idx;
  
  //Initialize the filter
  arma::mat B_LL = B0;
  arma::mat P_LL = P0;
  double lnl = 0.0;
  
  //Kalm filter routine
  for(int i = 0; i < yt.n_cols; i++){
    B_tl.col(i) = Dt + Ft * B_LL; //Initial estimate of unobserved values conditional on t-1
    P_tl.slice(i) = Ft * P_LL * Ft.t() + Qt; //Initial estimate of the covariance matrix conditional on t-1
    N_t.col(i) = yt.col(i) - At - Ht * B_tl.col(i) - beta * X.col(i); //Prediction error conditoinal on t-1
    nonna_idx = arma::find_finite(N_t.col(i));
    na_idx = arma::find_nonfinite(N_t.col(i));
    if(!na_idx.is_empty()){
      arma::uvec cols;
      cols = i;
      N_t(na_idx, cols) = arma::vec(na_idx.n_elem, arma::fill::zeros);
    }
    
    F_t.slice(i) = Ht * P_tl.slice(i) * Ht.t() + Rt; //Variance of the predictoin error conditional on t-1
    K_t.slice(i) = P_tl.slice(i) * Ht.t() * inv(F_t.slice(i)); //Kalman gain conditional on t-1
    B_tt.col(i) = B_tl.col(i) + K_t.slice(i) * N_t.col(i); //Final estimate of the unobserved values
    P_tt.slice(i) = P_tl.slice(i) - K_t.slice(i) * Ht * P_tl.slice(i); //Final estiamte of the covariance matrix
    lnl = lnl + 0.5*arma::as_scalar((log(det(F_t.slice(i))) +  N_t.col(i).t() * inv(F_t.slice(i)) * N_t.col(i)));
    
    //Reinitialize for the next iteration
    B_LL = B_tt.col(i);
    P_LL = P_tt.slice(i);
  }
  
  return Rcpp::List::create(Rcpp::Named("loglik") = -lnl,
                            Rcpp::Named("B_tl") = B_tl,
                            Rcpp::Named("B_tt") = B_tt,
                            Rcpp::Named("P_tl") = P_tl,
                            Rcpp::Named("P_tt") = P_tt,
                            Rcpp::Named("F_t") = N_t,
                            Rcpp::Named("N_t") = N_t,
                            Rcpp::Named("K_t") = K_t);
}

// [[Rcpp::export]]
Rcpp::List kalman_smoother(const arma::mat B_tl, arma::mat B_tt, const arma::cube P_tl, 
                           arma::cube P_tt, const arma::mat Ft){
  int t = B_tt.n_cols - 1;
  arma::mat Ptt_x_Ft_x_PtInv = P_tt.slice(t - 1) * Ft.t() * Rginv(P_tl.slice(t));
  
  for(int i = t - 1; i >= 0; i--){
    Ptt_x_Ft_x_PtInv = P_tt.slice(i) * Ft.t() * Rginv(P_tl.slice(i + 1));
    B_tt.col(i) = B_tt.col(i) + Ptt_x_Ft_x_PtInv * (B_tt.col(i + 1) - B_tl.col(i + 1));
    P_tt.slice(i) = P_tt.slice(i) + Ptt_x_Ft_x_PtInv * (P_tt.slice(i + 1) - P_tl.slice(i + 1)) * Ptt_x_Ft_x_PtInv.t();
  }
  return Rcpp::List::create(Rcpp::Named("B_tt") = B_tt,
                            Rcpp::Named("P_tt") = P_tt);
}


//RcppArmadillo.package.skeleton(name = "TCSDecomp", path = "/Users/alexhubbard/Dropbox (Opendoor)/R Codes/Packages")
//compileAttributes(verbose=TRUE)
//library(tools)
//package_native_routine_registration_skeleton("/Users/alexhubbard/Dropbox (Opendoor)/R Codes/Packages/TCSDecomp")
//git config remote.origin.url git@github.com:opendoor-labs/TCSDecomp.git

//sourceCpp("/Users/alexhubbard/Dropbox (Opendoor)/R Codes/Packages/kfdecomp/src/kalmanfilter.cpp")