admm_loop = function(W_update, Beta_update, W_init, Beta_init, white_coefs,
                     data_info, lambdas, admm_args,
                     add_stop_rule, fac_grp_lasso,
                     verbose, fit_hist){

  ## Set ADMM parameters ####
  max_iter = admm_args$max_iter; rho = admm_args$rho
  kappa = admm_args$kappa; abs_tol = admm_args$abs_tol
  rel_tol = admm_args$rel_tol; eps = admm_args$eps
  inner_verbose = admm_args$inner_verbose

  ## Read necessary arguments ####
  n_lams = length(lambdas)

  types_by_node = data_info$types_by_node
  n_resps_by_node = data_info$n_resps_by_node
  n_expls_by_node = data_info$n_expls_by_node
  n_resps = data_info$n_resps
  n_coefs = data_info$n_coefs
  n_nodes = data_info$n_nodes
  intcpt = data_info$intcpt
  idx_intcpt = ifelse(intcpt == "none", data_info$n_expls, data_info$n_expls + 1)

  ## Initialize parameters and set containers across lambdas ####
  W = W_init
  Beta = Beta_init
  Mu = Beta; Mu[] <- 0
  Beta_new_by_lam  = list(0)
  W_new_by_lam  = list(0)

  ## Set fit_hist container if specified ####
  if (fit_hist == TRUE) {history_Beta_by_lam = list(0); history_W_by_lam = list(0)}

  ## Path-wise loop starts ########
  for (ell in 1:n_lams){ # ell = 1

    if (verbose == T){ cat(paste0("Pathwise regularization for ", ell,
                                  ord(ell), " lambda starts ========\n")) }
    lambda = lambdas[ell]
    if (fit_hist == TRUE) history_Beta = list(0); history_W = list(0)

    r_norm = s_norm = eps_pri = eps_dual = c()
    k = 0
    conv = 0

    ## For the cold-start, uncomment this part
    # W = initialize_W(data_info)
    # Beta = Mu = matrix(0.0, idx_intcpt, n_resps)

    ### ADMM loop starts ========
    while (k < max_iter ){

      k = k + 1

      if (verbose == T){ cat("ADMM Iter =", k, "starts ======== \n") }

      #### W-update ----
      W_new <- W_update(W, Beta, Mu, rho, data_info, white_coefs, inner_verbose)
      W_hat <- kappa * W_new + (1-kappa) * Beta # over-relaxation

      #### Beta-update ----
      Beta_new <- Beta_update(W_hat, Mu, rho, lambda, fac_grp_lasso, data_info)

      #### Mu-update ----
      Mu_new <- Mu + (W_hat - Beta_new) # scaled version

      #### Termination Ready ----
      r_norm[k] = norm(W_new - Beta_new)
      s_norm[k] = norm(-rho * (Beta_new - Beta))
      eps_pri[k]= sqrt(n_coefs) * abs_tol + rel_tol * max(norm(W_new), norm(-Beta_new))
      eps_dual[k] = sqrt(n_coefs) * abs_tol + rel_tol * norm(rho * Mu_new)

      #### Termination check ----
      if ( (r_norm[k] < eps_pri[k] ) & (s_norm[k] < eps_dual[k]) ) {
        conv = 1
        if (verbose == T){ cat("ADMM converges \n") }
        break
      }

      W = W_new
      Beta = Beta_new
      Mu = Mu_new

      if (fit_hist == TRUE) {history_W[[k]] = W; history_Beta[[k]] = Beta}

    } ### end of ADMM loop ********

    ### Check convergence status ====
    if (verbose == T & conv != 1){
      cat("ADMM iteration hits the maximum number of iterations \n") }

    ### Save objects if history requested ====
    if (fit_hist == TRUE) {
      history_Beta_by_lam[[ell]] = history_Beta
      history_W_by_lam[[ell]] = history_W
    }

    ### Save objects as an pathwise solution ====
    colnames(W_new) <- colnames(Beta_new) <-
      Y_names(n_nodes, types_by_node, n_resps_by_node)

    if (intcpt == "none") {
      rownames(W_new) <- rownames(Beta_new) <-
        c(X_names(n_nodes, types_by_node, n_expls_by_node))
    } else {
      rownames(W_new) <- rownames(Beta_new) <-
        c(X_names(n_nodes, types_by_node, n_expls_by_node),"(Intcpt)")}

    Beta_new_by_lam[[ell]] = Beta_new
    W_new_by_lam[[ell]] = W_new

    ### Warm start ready ====
    W = W_new # may not necessary
    Beta = Beta_new
    Mu = matrix(0.0, idx_intcpt, n_resps)

    ### Additional stopping rule for pathwise solutions ===
    A_est = calc_A_est(Beta, data_info, intcpt, eps = eps)

    if (add_stop_rule == T){
      if (sum(A_est == 1) > 3 * n_nodes){
        cat("Current estimated n_edges is larger than 3*n_nodes.\n")
        cat("Pathwise regularization stops. \n")
        break # of current pathwise loop
      }
    }

  } # end of pathwise loop ********

  admm_results = c(list(Beta_est_by_lam = Beta_new_by_lam,
                        W_est_by_lam = W_new_by_lam, lambdas = lambdas),
                   if (fit_hist == T) list(history_W_by_lam = history_W_by_lam,
                                           history_Beta_by_lam = history_Beta_by_lam))

  return(admm_results)

}
