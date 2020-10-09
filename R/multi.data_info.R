multi.data_info = function(data_input, n_resps_by_node, intcpt = "always"){

  ## Read and set data ####
  n_obs = dim(data_input)[1]
  n_nodes = ncol(data_input)

  if (is.data.frame(data_input)) {
    if (!all(sapply(data_input, class) == "factor")){
      idx_conti_nodes = (sapply(data_input, class) != "factor")
      message("There exist continuous variables.")
      message("These will be coerced to multi-level data.")
      data_temp <- data_input
      data_input[, idx_conti_nodes] <-
        as.data.frame(lapply(data_temp[,idx_conti_nodes], as.factor))

      } # otherwise, good to go.

  } else {
    message("Input data is not a dataframe. It will be converted to a dataframe.")
    data_temp <- as.data.frame(data_input)
    data_input <- as.data.frame(lapply(data_temp, as.factor))
  }

  if (is.null(n_resps_by_node)){
    message("Numbers of levels for each node is not provided.")
    message("The information will be obtained via 'data_input'")

    n_resps_by_node = sapply(data_input, nlevels)

  } else {

    if (!all(n_resps_by_node == sapply(data_input, nlevels)))
      stop("n_resp_by_node does not coincide with the info from data_input")

  }

  types_by_node = rep("m", n_nodes)
  n_expls_by_node = n_resps_by_node - 1

  X1 = gen_X1_design(data_input, types_by_node, n_expls_by_node)
  Y = gen_Y_design(data_input, types_by_node, n_resps_by_node)

  idx_resps_by_node = gen_indices_by_node(n_resps_by_node)
  idx_expls_by_node = gen_indices_by_node(n_expls_by_node)

  Sr = t(diag(n_nodes)[rep(1:n_nodes, times = n_resps_by_node),])
  Se = t(diag(n_nodes)[rep(1:n_nodes, times = n_expls_by_node),])
  SSr = Sr[rep(1:n_nodes, times = n_resps_by_node),]

  n_expls = sum(n_expls_by_node)
  n_resps = sum(n_resps_by_node)
  n_coefs = n_expls * n_resps

  data_info = list( n_resps_by_node = n_resps_by_node, n_expls_by_node = n_expls_by_node,
                    types_by_node = types_by_node, intcpt = intcpt,
                    X1 = X1, Y = Y,
                    idx_resps_by_node = idx_resps_by_node,
                    idx_expls_by_node = idx_expls_by_node,
                    n_obs = n_obs, n_nodes = n_nodes, n_coefs = n_coefs,
                    n_expls = n_expls, n_resps = n_resps,
                    Sr = Sr, Se = Se, SSr = SSr )

  return(data_info)

}
