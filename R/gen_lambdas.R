gen_lambdas = function(data_info, n_lams = 30, eps = 0.05, ver = "1"){

  Sr = data_info$Sr
  Se = data_info$Se
  n_obs = data_info$n_obs
  X1 = data_info$X1
  Y = data_info$Y

  types_by_node = data_info$types_by_node
  conti_nodes = which(types_by_node == "c")
  multi_nodes = which(types_by_node == "m")
  idx_resps_by_node = data_info$idx_resps_by_node

  Ybar = apply(Y, 2, mean)
  Yc = sweep(Y, 2, Ybar, "-")
  X1tYc = t(X1) %*% Yc


  ## continuous part ####
  lambda_max_conti = NULL

  if (any(conti_nodes)) {

    conti_resp = unlist(idx_resps_by_node[(types_by_node == "c")])

    Yc_norm_sq = apply(Yc, 2, function(x) sum(x^2))
    scaled_X1tYc = sweep(X1tYc, 2, Yc_norm_sq, "/" )
    norm_scaled_X1tYc = sqrt(Se %*% (scaled_X1tYc^2) %*% t(Sr))
    diag(norm_scaled_X1tYc) <- 0

    lambda_max_conti = max(norm_scaled_X1tYc[,conti_nodes, drop = F])

  }


  ## multi-level part ####
  lambda_max_multi = NULL

  if (any(multi_nodes)){

    norm_X1tYc = sqrt(Se %*% (X1tYc^2) %*% t(Sr))
    diag(norm_X1tYc) <- 0

    multi_resp = unlist(idx_resps_by_node[multi_nodes])
    lambda_max_multi = max( (1/n_obs) * norm_X1tYc[,multi_nodes, drop = F] )

  }


  ## find the max ####
  lambda_max = max(lambda_max_conti, lambda_max_multi) * 1.01
  # ... *1.01 adjust lambda_max slightly larger to make it enough larger all the
  # ... all the coefs are zeros


  ## generate list of lambdas ####
  lambda_min = lambda_max * eps

  if (ver == "1"){
    lambdas = exp(seq(log(lambda_max), log(lambda_min), length.out = n_lams))
  } else if (ver == "2") {
    lambdas = seq(lambda_max, lambda_min, length.out = n_lams)
  }

  return(lambdas)

}
