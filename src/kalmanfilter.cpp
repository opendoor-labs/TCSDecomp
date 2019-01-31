#include <RcppArmadillo.h>
#include <Rcpp.h>

// Correctly setup the build environment
// [[Rcpp::depends(RcppArmadillo)]]

// Add a flag to enable OpenMP at compile time
// [[Rcpp::plugins(openmp)]]

// Protect against compilers without OpenMP
#ifdef _OPENMP
#include <omp.h>
#endif

// B_TT = B_TL = matrix(0, ncol = ncol(yt), nrow = nrow(sp$Ft))
// P_TT = P_TL = as.list(NA, ncol(yti))
//
// #Initialize the filter
// B_LL = sp$a0
// P_LL = sp$P0
//
// for(j in 1:ncol(yt)){
//   ################## Kalman filter routine ##################
//   B_TL[, j] = sp$Dt + sp$Ft %*% B_LL  #Initial estimate of unobserved values conditional on t-1
//   P_TL[[j]] = sp$Ft %*% P_LL %*% t(sp$Ft) + sp$Qt #Initial estimate of the covariance matrix conditional on t-1
//   N_TL = yt[, j] - sp$At - sp$Ht %*% B_TL[, j] #Prediction error conditoinal on t-1
//   F_TL = sp$Ht %*% P_TL[[j]] %*% t(sp$Ht) + sp$Rt #Variance of the predictoin error conditional on t-1
//   K_T = P_TL[[j]] %*% t(sp$Ht) %*% ginv(F_TL) #Kalman gain conditional on t-1
//   B_TT[, j] = B_TL[, j] + K_T %*% N_TL #Final estimate of the unobserved values
//   P_TT[[j]] = P_TL[[j]] - K_T %*% sp$Ht %*% P_TL[[j]] #Final estiamte of the covariance matrix
//
//   #Reinitialize for the next iteration
//   B_LL = B_TT[, j]
//   P_LL = P_TT[[j]]
// }
//
// #Kalman Smoother
// for(j in (ncol(yt) - 1):1){
//   B_TT[, j] = B_TT[, j] + P_TT[[j]] %*% t(sp$Ft) %*% ginv(P_TL[[j + 1]]) %*% (B_TT[, j + 1] - B_TL[, j + 1])
//   P_TT[[j]] = P_TT[[j]] + P_TT[[j]] %*% t(sp$Ft) %*% ginv(P_TL[[j + 1]]) %*% (P_TT[[j + 1]] - P_TL[[j + 1]]) %*% t(P_TT[[j]] %*% t(sp$Ft) %*% ginv(P_TL[[j + 1]]))
// }

// [[Rcpp::export]]
Rcpp::List kf(arma::mat B0, arma::mat P0, arma::mat Dt, arma::mat At,
              arma::mat Ft, arma::mat Ht, arma::mat Qt, arma::mat Rt,
              arma::mat yt){

  //Define the storage matrices
  arma::mat B_TT(Ft.n_rows, yt.n_cols);
  arma::mat B_TL(B_TT.n_rows, yt.n_cols);
  arma::cube P_TT(Ft.n_rows, Ft.n_rows, yt.n_cols);
  arma::cube P_TL(Ft.n_rows, Ft.n_rows, yt.n_cols);
  arma::mat N_T(yt.n_rows, yt.n_cols);
  arma::cube F_T(yt.n_rows, yt.n_rows, yt.n_cols);
  arma::cube K_T(Ft.n_rows, yt.n_rows, yt.n_cols);

  //Initialize the filter
  arma::mat B_LL = B0;
  arma::mat P_LL = P0;
  double lnl = 0.0;

  //Kalm filter routine
  for(int i = 0; i < yt.n_cols; i++){
    B_TL.col(i) = Dt + Ft * B_LL; //Initial estimate of unobserved values conditional on t-1
    P_TL.slice(i) = Ft * P_LL * Ft.t() + Qt; //Initial estimate of the covariance matrix conditional on t-1
    N_T.col(i) = yt.col(i) - At - Ht * B_TL.col(i); //Prediction error conditoinal on t-1
    F_T.slice(i) = Ht * P_TL.slice(i) * Ht.t() + Rt; //Variance of the predictoin error conditional on t-1
    K_T.slice(i) = P_TL.slice(i) * Ht.t() * inv(F_T.slice(i)); //Kalman gain conditional on t-1
    B_TT.col(i) = B_TL.col(i) + K_T.slice(i) * N_T.col(i); //Final estimate of the unobserved values
    P_TT.slice(i) = P_TL.slice(i) - K_T.slice(i) * Ht * P_TL.slice(i); //Final estiamte of the covariance matrix
    lnl = lnl + 0.5*arma::as_scalar((log(det(F_T.slice(i))) +  N_T.col(i).t() * inv(F_T.slice(i)) * N_T.col(i)));

    //Reinitialize for the next iteration
    B_LL = B_TT.col(i);
    P_LL = P_TT.slice(i);
  }

  return Rcpp::List::create(Rcpp::Named("loglik") = -lnl,
                            Rcpp::Named("B_TL") = B_TL,
                            Rcpp::Named("B_TT") = B_TT,
                            Rcpp::Named("P_TL") = P_TL,
                            Rcpp::Named("P_TT") = P_TT,
                            Rcpp::Named("F_T") = N_T,
                            Rcpp::Named("N_T") = N_T,
                            Rcpp::Named("K_T") = K_T);
}

// [[Rcpp::export]]
Rcpp::List kf_smoother(arma::mat B_TL, arma::mat B_TT, arma::cube P_TL, arma::cube P_TT, arma::mat Ft){
  int t = B_TT.n_cols - 1;
  arma::mat Ptt_x_Ft_x_PtInv = P_TT.slice(t - 1) * Ft.t() * pinv(P_TL.slice(t));

  for(int i = t - 1; i >= 0; i--){
    Ptt_x_Ft_x_PtInv = P_TT.slice(i) * Ft.t() * pinv(P_TL.slice(i + 1));
    B_TT.col(i) = B_TT.col(i) + Ptt_x_Ft_x_PtInv * (B_TT.col(i + 1) - B_TL.col(i + 1));
    P_TT.slice(i) = P_TT.slice(i) + Ptt_x_Ft_x_PtInv * (P_TT.slice(i + 1) - P_TL.slice(i + 1)) * Ptt_x_Ft_x_PtInv.t();
  }
  return Rcpp::List::create(Rcpp::Named("B_TT") = B_TT,
                            Rcpp::Named("P_TT") = P_TT);
}


//RcppArmadillo.package.skeleton(name = "kfdecomp", path = "/Users/alexhubbard/Dropbox (Opendoor)/R Codes/Packages")
//compileAttributes(verbose=TRUE)
//library(tools)
//package_native_routine_registration_skeleton("/Users/alexhubbard/Dropbox (Opendoor)/R Codes/Packages/TCSDecomp")
//git config remote.origin.url git@github.com:opendoor-labs/TCSDecomp.git

//sourceCpp("/Users/alexhubbard/Dropbox (Opendoor)/R Codes/Packages/kfdecomp/src/kalmanfilter.cpp")