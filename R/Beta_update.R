grp_lasso.Beta_update = function(W_hat, Mu, rho, lambda, fac_grp_lasso, data_info){

  idx_resps_by_node = data_info$idx_resps_by_node
  idx_expls_by_node = data_info$idx_expls_by_node
  n_nodes = data_info$n_nodes

  Beta_temp = W_hat + Mu
  Beta_new = Beta_temp

  for (j in 1:n_nodes){ # j = 3

    idx_resps_j = idx_resps_by_node[[j]]

    for (i in 1:n_nodes){ # i = 4
      
      if (i == j) next

      idx_expls_j = idx_expls_by_node[[i]]
      num_elmt = length(Beta_temp[idx_expls_j,idx_resps_j])
      fac = ifelse(fac_grp_lasso, sqrt(num_elmt), 1)
      Beta_new[idx_expls_j, idx_resps_j] = 
        grp_shrinkage(Beta_temp[idx_expls_j, idx_resps_j], (fac*lambda)/rho)

    }

  }

  return(Beta_new)

}
