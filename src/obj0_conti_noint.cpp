#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
double obj0cpp_conti_noint(arma::vec w, arma::vec beta, arma::vec mu, double rho,
                           double alpha, double rho1,
                           arma::mat X, arma::mat Y, arma::mat white_coefs,
                           int idx_intcpt, int n_resps,
                           double n_obs, double n_nodes, int n_coefs,
                           arma::mat Sr, arma::mat Se) {

  arma::vec w1_temp = w;
  arma::mat W1_temp = reshape(w1_temp, idx_intcpt, n_resps);
  arma::mat W_temp = W1_temp;
  arma::mat W = W_temp % white_coefs;
  arma::mat W1 = W;

  arma::mat H = X * W;

  double loss = (0.5) * (1.0/double(n_obs)) * accu(sum(pow(Y - H, 2),0));

  arma::mat pp = Se * (W1 % W1) * trans(Sr);
  double h = trace(arma::expmat(pp)) - n_nodes;
  arma::vec vector = w - beta + mu;

  double result = loss + (rho1/2.0) * pow(h, 2) + alpha * h +
    (rho/2.0) * pow(norm(vector, 2), 2);

  return(result);
}
