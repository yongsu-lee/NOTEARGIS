#' Generate a dataset using a true graph and corresponding adjacency matrix
#'
#' @param n_obs positive interge scalar, the desired number of observations.
#' @param A_true binary matrix, a true adjacecny matrix from a desired graph.
#' @param W_true double matrix, a true parameter matrix.
#' @param seed integer for seed control, default is 1.
#'
#' @return a positive integer or double matrix, a dataset with the desired number of observations.
#' @export
#' @importFrom igraph topo_sort
#' @importFrom stats model.matrix rnorm
gen_data = function(n_obs, A_true, W_true, seed = 1){

  if (missing(n_obs)) stop("Number of observation?")
  if (missing(A_true)) stop("True adjacency matrix?")
  if (missing(W_true)) stop("True paramter?")

  types_by_node = attr(W_true, "types_by_node")
  graph_true = graph_from_adjacency_matrix(A_true)

  ## grouped continuous case
  if (all(types_by_node == "c")){

    Z = conti.gen_data(n_obs, A_true, graph_true, W_true, seed, T)


  ## multinomial case
  } else if (all(types_by_node == "m")) {

    Z = multi.gen_data(n_obs, A_true, graph_true, W_true, seed, T)

  }

  return(Z)

}
