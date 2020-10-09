conti.W_update = function(W, Beta, Mu, rho, data_info, white_coefs, verbose,
                          rho1_max = 1e+20, h_tol = 1e-8){

  ## Read data ####
  n_expls = data_info$n_expls
  intcpt = data_info$intcpt
  idx_intcpt = ifelse(intcpt=="none", data_info$n_expls, data_info$n_expls + 1)

  n_resps = data_info$n_resps
  X1 = data_info$X1

  Y = data_info$Y

  Sr = data_info$Sr
  Se = data_info$Se

  n_obs = data_info$n_obs
  n_coefs = data_info$n_coefs
  n_nodes = data_info$n_nodes

  if (intcpt == "none") {
    X = X1
  } else{
    X = cbind(X1, rep(1,n_obs))
    colnames(X)[idx_intcpt] <- "(Intcpt)"
  }

  ## Initialize parameters ####
  W = W * white_coefs;

  if(intcpt == "none") { w = c(W) } else { vectorizing(W, idx_intcpt) }
  if(intcpt == "none") { beta = c(Beta) } else { vectorizing(Beta, idx_intcpt)}
  if(intcpt == "none") { mu = c(Mu) } else {  vectorizing(Mu, idx_intcpt) }

  alpha = 0.0
  rho1 = 1.0
  h = Inf

  if (intcpt == "none"){

    obj1 = function(w, alpha, rho1){
      obj0cpp_conti_noint(w, beta, mu, rho, alpha, rho1, X, Y, white_coefs,
                          idx_intcpt, n_resps, n_obs, n_nodes, n_coefs,
                          Sr, Se)}
    grad1 = function(w, alpha, rho1){
      grad0cpp_conti_noint(w, beta, mu, rho, alpha, rho1, X, Y, white_coefs,
                           idx_intcpt, n_resps, n_obs, n_nodes, n_coefs,
                           Sr, Se)}

  } else { # for non-intercept case

    obj1 = function(w, alpha, rho1){
      obj0cpp_conti_int(w, beta, mu, rho, alpha, rho1, X, Y, white_coefs,
                        idx_intcpt, n_resps, n_obs, n_nodes, n_coefs,
                        Sr, Se)}
    grad1 = function(w, alpha, rho1){
      grad0cpp_conti_int(w, beta, mu, rho, alpha, rho1, X, Y, white_coefs,
                         idx_intcpt, n_resps, n_obs, n_nodes, n_coefs,
                         Sr, Se)}
  }


  while (rho1 < rho1_max) {

    if (verbose == TRUE){ cat("rho 1 =", rho1, " and alpha = ", alpha, "\n") }


    obj2 = function(w){obj1(w, alpha, rho1)}
    grad2 = function(w){grad1(w, alpha, rho1)}

    w_new_temp = lbfgs::lbfgs(obj2, grad2, w, invisible = 1)$par

    if (intcpt == "none"){
      W_new_temp = matricizing(w_new_temp, n_expls, n_resps, coef_only = T)
      W_new = W_new_temp * white_coefs
      W1_new = W_new
    } else {
      W_new_temp = matricizing(w_new_temp, n_expls, n_resps)
      W_new = W_new_temp * white_coefs
      W1_new = W_new[-idx_intcpt, ] # coef part only
    }

    h_new = sum(diag(expm::expm(Se %*% W1_new^2 %*% t(Sr)))) - n_nodes

    if (h_new > 0.25 * h) rho1 = rho1 * 10

    W = W_new;
    if (intcpt == "none") w = c(W) else w = vectorizing(W, idx_intcpt)

    h = h_new

    if (h <= h_tol) break

    alpha <- alpha + rho1 * h

  }

  return(W_new)

}
