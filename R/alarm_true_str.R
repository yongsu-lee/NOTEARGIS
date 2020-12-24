#' The ALARM network true structure information
#'
#' True network information of ALARM network data from the
#' \href{https://www.bnlearn.com/bnrepository/discrete-medium.html#alarm}{Bayesian network repository}.
#'
#' @format A list object including:
#' \describe{
#' \item{A_true}{adjacency matrix of the true network}
#' \item{graph_true}{igraph object based on the adjacency matrix}
#' \item{n_edges}{the number of true edges}
#' }
#'
#' @usage
#' data(alarm_true_str)
#'
#' @references
#' [1] Beinlich IA, Suermondt HJ, Chavez RM, Cooper GF (1989). “The ALARM monitoring system: A case study with two probabilistic inference techniques for belief networks.” In AIME 89, pp. 247–256. Springer.
"alarm_true_str"
