initialize_W = function(data_info, intcpt = "always", black_list_nodes=c()){

  Y= data_info$Y
  types_by_node = data_info$types_by_node

  n_nodes = data_info$n_nodes
  idx_resps_by_node = data_info$idx_resps_by_node
  idx_expls_by_node = data_info$idx_expls_by_node

  idx_intcpt = ifelse(intcpt == "none", data_info$n_expls, data_info$n_expls + 1)
  n_resps = data_info$n_resps

  W_init = matrix(0, ncol = n_resps, nrow = idx_intcpt)
  white_coefs = matrix(1, ncol = n_resps, nrow = idx_intcpt)

  idx_resps_conti_only = unlist(idx_resps_by_node[(types_by_node == "c")])
  idx_resps_multi_only = unlist(idx_resps_by_node[(types_by_node == "m")])

  if (intcpt != "none"){
  W_init[idx_intcpt, idx_resps_conti_only] =
    apply(cbind(Y[ , idx_resps_conti_only]), 2, mean )
  W_init[idx_intcpt, idx_resps_multi_only] =
    apply(cbind(Y[ , idx_resps_multi_only]), 2,
                                     function(x) log(mean(x)/(1-mean(x))) )
  }

  for (j in 1:n_nodes){ # j = 5

    idx_resps_j = idx_resps_by_node[[j]];

    if (any(black_list_nodes %in% j)){
      white_coefs[, idx_resps_j] <- 0

    } else {
      idx_expls_j = idx_expls_by_node[[j]]
      white_coefs[idx_expls_j, idx_resps_j] <- 0
    }

  }

  init_W_info = list(W_init = W_init, white_coefs = white_coefs)

  return(init_W_info)

}
