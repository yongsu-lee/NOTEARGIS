multi.gen_data=function(n_obs, A_true, graph_true, W_true, seed=1, df = T){

  types_by_node = attr(W_true, "types_by_node")
  n_levels_by_node = attr(W_true, "n_resps_by_node")
  n_dummys_by_node = attr(W_true, "n_dummys_by_node")

  n_nodes = length(types_by_node)
  parents_by_node = lapply(as.data.frame(A_true), function(x) which(x==1))

  n_expls = nrow(W_true) - 1
  idx_intcpt = nrow(W_true)

  idx_resps_by_node = gen_indices_by_node(n_levels_by_node)
  idx_expls_by_node = gen_indices_by_node(n_dummys_by_node)

  X1_temp = matrix(NA, n_obs, n_expls)
  Z = matrix(NA, n_obs, n_nodes)

  topo_info = topo_sort(graph_true)

  set.seed(seed)

  for (j in 1:n_nodes){ # j = 2

    resp_node = topo_info[j]
    idx_resps_j = idx_resps_by_node[[resp_node]]
    n_resps_j = length(idx_resps_j)
    expl_nodes = parents_by_node[[resp_node]]

    if (length(expl_nodes) == 0L){

      q_j = n_levels_by_node[resp_node]

      try_gen_levels = 0
      while (try_gen_levels < 100){

        Z_j = factor(sample(0:(q_j-1), n_obs, replace = T ))

        if (n_levels_by_node[resp_node] == nlevels(Z_j)) break

        try_gen_levels = try_gen_levels + 1
        if (try_gen_levels > 100){
          stop("Something's wrong. Single level column has been generated.")}

      } ## end of while (try_gen_levels < 100)

      ## X matrix to generate other nodes data
      X1_temp[,idx_expls_by_node[[resp_node]]] = model.matrix(~factor(Z_j))[,-1]
      Z[, resp_node] = Z_j

    } else {

      idx_expls_j = unlist(idx_expls_by_node[expl_nodes])
      n_expls_j = length(idx_expls_j)

      q_j = n_levels_by_node[resp_node]

      etas1 = X1_temp[, idx_expls_j, drop = F] %*% (W_true[idx_expls_j, idx_resps_j])

      ## Add intercepts
      etas = sweep(etas1, 2, W_true[idx_intcpt, idx_resps_j], "+")
      probs = t(apply(exp(etas), 1, function(x) (x/sum(x))))

      try_gen_levels=0
      while (try_gen_levels < 100){

        Z_j = factor(apply(probs, 1, function(x)
          sample(0:(q_j-1), 1, replace=T, prob = x )) )

        if (n_levels_by_node[resp_node] == nlevels(Z_j) ) break

        try_gen_levels = try_gen_levels + 1
        if (try_gen_levels > 100){
          stop("Something's wrong. Single level column has been generated.")
        }
      }

      X1_temp[,idx_expls_by_node[[resp_node]]] = model.matrix(~factor(Z_j))[,-1]
      Z[, resp_node] = Z_j

    }

  } # end of for loop (j in 1:n_nodes)

  if (df == T){
    Z <- as.data.frame(Z)
    Z[, types_by_node=="m"] <- lapply(Z[,types_by_node=="m", drop = F], as.factor)
  }

  colnames(Z) <- paste0("Z", 1:n_nodes)
  return(Z)

}
