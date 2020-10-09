#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::vec grad0cpp_conti_int(arma::vec w, arma::vec beta, arma::vec mu,
                             double rho, double alpha, double rho1,
                             arma::mat X, arma::mat Y, arma::mat white_coefs,
                             int idx_intcpt, int n_resps,
                             double n_obs, double n_nodes, int n_coefs,
                             arma::mat Sr, arma::mat Se) {

  arma::vec w0_temp = w.tail(n_resps); arma::vec w1_temp = w.head(n_coefs);
  arma::mat W0_temp = reshape(w0_temp, 1, n_resps);
  arma::mat W1_temp = reshape(w1_temp, idx_intcpt-1, n_resps);
  arma::mat W_temp = join_cols(W1_temp, W0_temp);
  arma::mat W = W_temp % white_coefs;
  arma::mat W1 = W.head_rows(idx_intcpt-1);

  arma::mat H = X * W;

  // grad for continuous part

  arma::rowvec Deno = sum(pow(Y - H, 2),0);
  arma::mat grad_loss_conti_mat = - trans(X) * (Y - H);
  grad_loss_conti_mat.each_row() /= Deno;

  // combined grad loss
  arma::mat grad_loss_mat = grad_loss_conti_mat % white_coefs;

  // penalty part
  arma::mat pp = Se * (W1 % W1) * trans(Sr);
  arma::mat E = arma::expmat(pp);
  double h = trace(E) - n_nodes;

  // grad_h_mat
  arma::mat grad_h_mat(idx_intcpt, n_resps, arma::fill::zeros);
  arma::mat SetEtSr = trans(Se) * trans(E) * Sr;
  arma::mat grad_h_mat_coef = SetEtSr % (2.0 * W1);
  grad_h_mat.head_rows(idx_intcpt-1) = grad_h_mat_coef;

  // grad_mat -> grad_vec
  arma::mat grad_mat = grad_loss_mat + grad_h_mat * (rho1 * h + alpha);
  arma::mat grad_mat1 = grad_mat.head_rows(idx_intcpt-1), grad_mat0 = grad_mat.tail_rows(1);
  arma::vec grad_vec = join_cols(vectorise(grad_mat1), vectorise(grad_mat0));

  arma::vec result = grad_vec + rho * (w - beta + mu);

  return(result);

}
