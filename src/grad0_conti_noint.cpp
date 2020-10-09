#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::vec grad0cpp_conti_noint(arma::vec w, arma::vec beta, arma::vec mu, double rho,
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

  arma::mat grad_loss_mat = (1.0/double(n_obs)) * ((- trans(X) * (Y - H)) % white_coefs);

  // penalty part
  arma::mat pp = Se * (W1 % W1) * trans(Sr);
  arma::mat E = arma::expmat(pp);
  double h = trace(E) - n_nodes;

  // grad_h_mat
  arma::mat SetEtSr = trans(Se) * trans(E) * Sr;
  arma::mat grad_h_mat_coef = SetEtSr % (2.0 * W1);
  arma::mat grad_h_mat = grad_h_mat_coef;

  // grad_mat -> grad_vec
  arma::mat grad_mat = grad_loss_mat + grad_h_mat * (rho1 * h + alpha);
  arma::vec grad_vec = vectorise(grad_mat);

  arma::vec result = grad_vec + rho * (w - beta + mu);

  return(result);
}
