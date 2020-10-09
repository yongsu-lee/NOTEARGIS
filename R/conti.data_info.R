conti.data_info = function(data_input, n_resps_by_node, intcpt = "none"){

  ## Read and set data ####
  if (is.data.frame(data_input)) {
    if (any(sapply(data_input, class) == "factor")){
      idx_multi_nodes = (sapply(data_input, class) == "factor")
      message("There exist multi-level variables." )
      message("These will be coerced to continuous data.")
      data_temp <- data_input
      data_input[, idx_multi_nodes] <-
        as.data.frame(lapply(data_temp[,idx_multi_nodes], as.numeric))

    } # otherwise, good to go.

  } else {
    message("Input data is not a dataframe. It will be converted to a dataframe.")
    data_input <- as.data.frame(data_input)
  }

  if (is.null(n_resps_by_node)){ # default: no group
    message("No 'n_resps_by_node' input. Data will be assumed to be a non-group case.")
    n_nodes = ncol(data_input)
    n_resps_by_node = rep(1, n_nodes)
  } else {
    n_nodes = length(n_resps_by_node)
  }

  n_obs = dim(data_input)[1]
  types_by_node = rep("c", n_nodes)

  n_expls_by_node = n_resps_by_node

  Y <- as.matrix(data_input)
  colnames(Y) =  Y_names(n_nodes, types_by_node, n_resps_by_node)
  X1 <- as.matrix(data_input)
  colnames(X1) =  X_names(n_nodes, types_by_node, n_expls_by_node)

  idx_resps_by_node = gen_indices_by_node(n_resps_by_node)
  idx_expls_by_node = idx_resps_by_node

  Sr = t(diag(n_nodes)[rep(1:n_nodes, times = n_resps_by_node),])
  Se = Sr

  n_resps = sum(n_resps_by_node)
  n_expls = n_resps

  n_coefs = n_expls * n_resps

  data_info = list( n_resps_by_node = n_resps_by_node, n_expls_by_node = n_expls_by_node,
                    types_by_node = types_by_node, intcpt = intcpt,
                    X1 = X1, Y = Y,
                    idx_resps_by_node = idx_resps_by_node,
                    idx_expls_by_node = idx_expls_by_node ,
                    n_obs = n_obs, n_nodes = n_nodes, n_coefs = n_coefs,
                    n_expls = n_expls, n_resps = n_resps,
                    Sr = Sr, Se = Se)

  return(data_info)

}
