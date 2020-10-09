conti.gen_data = function(n_obs, A_true, graph_true, W_true, seed, df = T) {

  types_by_node = attr(W_true, "types_by_node")
  n_resps_by_node = attr(W_true, "n_resps_by_node")
  n_expls_by_node = attr(W_true, "n_dummys_by_node")

  n_nodes = length(types_by_node)
  parents_by_node = lapply(as.data.frame(A_true), function(x) which(x==1))

  n_resps = sum(n_resps_by_node)
  n_expls = sum(n_expls_by_node)
  idx_intcpt = nrow(W_true)
  if (n_expls == idx_intcpt) idx_intcpt = NA

  idx_resps_by_node = gen_indices_by_node(n_resps_by_node)
  idx_expls_by_node = idx_resps_by_node

  X1_temp = matrix(NA, n_obs, n_expls)
  Z = matrix(NA, n_obs, n_resps)

  topo_info = topo_sort(graph_true)

  set.seed(seed)

  for (j in 1:n_nodes){

    resp_node = topo_info[j]
    idx_resps_j = idx_resps_by_node[[resp_node]]
    n_resps_j = length(idx_resps_j)
    expl_nodes = parents_by_node[[resp_node]]

    if (length(expl_nodes) == 0L){

      Z_j = matrix(rnorm(n_obs * n_resps_j), ncol = n_resps_j )

      if (!is.na(idx_intcpt)) sweep(Z_j, 2, W_true[idx_intcpt, idx_resps_j], "+")

      X1_temp[, idx_expls_by_node[[resp_node]]] = Z_j
      Z[, idx_resps_j] = Z_j

    } else {

      idx_expls_j = unlist(idx_expls_by_node[expl_nodes])
      n_expls_j = length(idx_expls_j)

      Z_j = X1_temp[, idx_expls_j, drop=F] %*% (W_true[idx_expls_j, idx_resps_j]) +
        matrix(rnorm(n_obs * n_resps_j ), ncol = n_resps_j )
      if (!is.na(idx_intcpt)) sweep(Z_j, 2, W_true[idx_intcpt, idx_resps_j], "+")

      X1_temp[, idx_expls_by_node[[resp_node]]] = Z_j
      Z[, idx_resps_j] = Z_j
    }

  } # end of for loop (j in 1:n_nodes)

  if (df == T) Z <- as.data.frame(Z)

  var_name_temp = rep(paste0("Z",1:n_nodes), times= n_resps_by_node)
  elmts_temp = as.list(n_resps_by_node)
  elmts = lapply(elmts_temp, function(x) paste0(":",1:x))
  elmts[types_by_node=="c" & n_resps_by_node==1] <- ""
  elmts_names = unlist(elmts)
  var_names = paste0(var_name_temp, elmts_names)
  colnames(Z) <- var_names

  return(Z)

}

