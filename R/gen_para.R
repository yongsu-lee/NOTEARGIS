#' Generate true parameters using an adjacency matrix
#'
#' @param A_true binary matrix, a true adjacecny matrix from a desired graph.
#' @param types_by_node blah
#' @param n_resps_by_node blah
#' @param intcpt blah
#' @param prob_sparse double scalar, define a sparseness of parameters within each parameter sub-matrix with \eqn{0\le prob_sparse \le 1} - default is 0. If \eqn{prob_sparse = 0}, all the parameter for each submatrix will be non-zeros. Not recommnedable, but, if \eqn{prob_sparse = 1}, all the parameter for each submatrix will be zeros. See details below for the parameter sub-matrix.
#' @param range a possible range of (absolute value of) parameters.
#' @param seed blah
#'
#' @return double matrix, true parameter matrix, \code{W_true} is returend with \code{n_levels} or \code{n_group_elmts} as an attribute.
#' @export
#' @importFrom utils head tail
#' @importFrom stats runif model.matrix
#' @importFrom igraph graph_from_adjacency_matrix is_dag
#' @details The parameter matrix can be divided into submatrices. For example, a submatrix (i,j) consists of parameters related to (node i) -> (node j). The \code{prob_sparse} controls the sparseness of parameters within each submatrix.
gen_para=function(A_true, types_by_node=NULL, n_resps_by_node=NULL,
                  intcpt=c("always","for_child_only","none"),
                  range=c(0.5,1), prob_sparse=0, seed=1){


  if (missing(intcpt)) {

    if (any(types_by_node == "m")) intcpt = "always"
    if (all(types_by_node == "c")) intcpt = "none"
  } else {
    intcpt = match.arg(intcpt)
  }

  if (any(types_by_node == "m") & intcpt == "none"){
    message("Multinomial DAG does not support non-intercept model.")
    message("intcpt options will be coerced to 'always'")
    intcpt = "always"
  }


  graph_temp = igraph::graph_from_adjacency_matrix( A_true )

  if (!igraph::is_dag(graph_temp)) {
    stop("Provided adjacency matrix does not represent a DAG. Check the A_true.")
  }

  n_nodes = length(types_by_node)
  n_dummys_by_node = ifelse(types_by_node == "c",
                            n_resps_by_node, n_resps_by_node-1)

  idx_expls_by_node = gen_indices_by_node(n_dummys_by_node)
  idx_resps_by_node = gen_indices_by_node(n_resps_by_node)

  n_row = sum(n_dummys_by_node)
  # n_row = sum(sapply(idx_expls_by_node, length))
  n_col = sum(n_resps_by_node)
  # n_col = sum(sapply(idx_resps_by_node, length))

  W_temp = matrix(0, n_row, n_col)

  set.seed(seed)
  for (j in 1:n_nodes){ #j = 1

    idx_resps_j = idx_resps_by_node[[j]]

    for (i in 1:n_nodes){ # i = 2

      idx_expls_i = idx_expls_by_node[[i]]

      if (A_true[i,j] == 1){

        n_cell = length(idx_expls_i) * length(idx_resps_j)
        W_temp[idx_expls_i, idx_resps_j] = sample(c(0, 1), n_cell, replace = T,
                                                  prob = c(prob_sparse, 1-prob_sparse))
      }
    }
  }


  if (intcpt == "for_child_only") {

    idx_intcpts = which(apply(W_temp, 2, sum) != 0)
    W_true = rbind(W_temp, rep(0,n_col))

    for (j in idx_intcpts) W_true[n_row + 1, j] = 1

  } else if (intcpt == "always") {

    W_true = rbind(W_temp, rep(1,n_col))

  } else {

    W_true = W_temp
  }

  # Assign generated coefficients  ----
  min_para = min(range); max_para = max(range)
  idx_nnz = which(W_true  == 1)
  n_nnz = sum(W_true)
  coefs = round(runif(n_nnz, min_para, max_para), 4)

  # Flipping the signs ----
  signs = sample(c(-1,1), n_nnz, replace = T)
  W_true[idx_nnz] = coefs * signs

  # Attaching types of nodes to W_true ----
  attr(W_true, "types_by_node") = types_by_node
  attr(W_true, "n_resps_by_node") = n_resps_by_node
  attr(W_true, "n_dummys_by_node") = n_dummys_by_node

  # if(any(types_by_node == "c" & n_resps_by_node > 1)) grp = T


  colnames(W_true) <- Y_names(n_nodes, types_by_node, n_resps_by_node)

  if (intcpt == "none") {
    rownames(W_true) <- c(X_names(n_nodes, types_by_node, n_dummys_by_node))

  } else {
    rownames(W_true) <- c(X_names(n_nodes, types_by_node, n_dummys_by_node),
                          "(Intcpt)")}

  return(W_true)

}
