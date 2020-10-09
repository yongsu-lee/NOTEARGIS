multi.W_update = function(W, Beta, Mu, rho, data_info, white_coefs, verbose,
                          rho1_max = 1e+20, h_tol = 1e-8){

  n_expls = data_info$n_expls
  idx_intcpt = data_info$n_expls + 1
  n_resps = data_info$n_resps
  X1 = data_info$X1

  Y = data_info$Y

  Sr = data_info$Sr
  Se = data_info$Se
  SSr = data_info$SSr

  n_obs = data_info$n_obs
  n_coefs = data_info$n_coefs
  n_nodes = data_info$n_nodes

  X = cbind(X1, rep(1,n_obs))
  colnames(X)[idx_intcpt] <- "(Intcpt)"

  ## Initialize parameters ####
  W = W * white_coefs;
  w = vectorizing(W * white_coefs, idx_intcpt)
  beta = vectorizing(Beta, idx_intcpt)
  mu = vectorizing(Mu, idx_intcpt)

  alpha = 0.0
  rho1 = 1.0
  h = Inf

  obj1 = function(w, alpha, rho1){
    obj0cpp_multi(w, beta, mu, rho, alpha, rho1, X, Y, white_coefs,
                  idx_intcpt, n_resps, n_obs, n_nodes, n_coefs,
                  Sr, Se)}
  grad1 = function(w, alpha, rho1){
    grad0cpp_multi(w, beta, mu, rho, alpha, rho1, X, Y, white_coefs,
                   idx_intcpt, n_resps, n_obs, n_nodes, n_coefs,
                   Sr, Se, SSr)}

  while (rho1 < rho1_max) {

    if (verbose == 1){ cat("rho 1 =", rho1, " and alpha = ", alpha, "\n") }

    obj2 = function(w){obj1(w, alpha, rho1)}
    grad2 = function(w){grad1(w, alpha, rho1)}

    w_new_temp = lbfgs::lbfgs(obj2, grad2, w, invisible = 1)$par

    W_new_temp = matricizing(w_new_temp, n_expls, n_resps)
    W_new = W_new_temp * white_coefs

    W1_new = W_new[-idx_intcpt, ] # coef part only
    h_new = sum(diag(expm::expm(Se %*% W1_new^2 %*% t(Sr)))) - n_nodes

    if (h_new > 0.25 * h) rho1 = rho1 * 10

    w = vectorizing(W_new, idx_intcpt) # w_new = vectorizing(W_new, idx_intcpt)
    h = h_new

    if (h <= h_tol) break

    alpha <- alpha + rho1 * h

  }

  return(W_new)
}


