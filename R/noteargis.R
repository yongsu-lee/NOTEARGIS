#' Learning a sparse DAG with Grouped Variables
#'
#' @param data_input a data.frame with attributes \code{data_type} and either \code{n_levels} or \code{n_grp_elmts}. Check \code{NOTEARGIS::gen_data} for details.
#' @param data_type a character vector, elements consist of either "c" or "m" indicating continuous variable and multinomial response, respectively.
#' @param n_resps_by_node positive integers vector, a list of numbers of elements for each group or numbers of levels for each multinomial variable.
#' @param intcpt intercept terms will be included if "always" and will be excluded if "none".
#' @param lambdas a double vector, user specified tuning parameter sequence. Typical usage is to have the function compute its own list of tuning parameters based on \code{n_lams} (and internal function \code{gen_lambdas}). Do not recommend to supply this argument. Default is \code{NA}.
#' @param n_lams a positive integer, the number of seuqences of tuning parameters. Default is \code{30}.
#' @param admm_args a list of value, to customize the several arguments for ADMM algorithm. See \code{?admm_arg_ctrl} for details.
#' @param add_stop_rule a logical scalar, if \code{add_stop_rule = TRUE}, the function will stop proceed fitting if the number of estimated edges are larger than \code{3*n_nodes}. Default is \code{TRUE}.
#' @param fac_grp_lasso a logical scalar, if \code{fac_grp_lasso = TRUE}, square-root of number of parameters will be considered as the multiplied factor in the group lasso update. Default is \code{FALSE}.
#' @param verbose a logical scalar, if \code{verbose = TRUE} you can see fitting status message in the console window. Default is \code{FALSE}.
#' @param fit_hist a logical scalar, if \code{verbose = TRUE} function returns all the parameters at each step of the ADMM loop. Do not recommend. Default is \code{FALSE}.
#'
#' @return An object with S3 class \code{"noteargis"} \item{A_est_by_lam}{A list of estimated adjacency matrices over the pathwise solutions. \code{graph_est} is saved as an attribute for each pathwise solution.} \item{Beta_new_by_lam}{A list of Beta estimates over the pathwise solutions.} \item{W_new_by_lam}{A list of W estimates over the pathwise solutions.} \item{lambdas}{A list of tuning parameters used for fitting all the pathwise solutions}\item{history_W_by_lam}{(if \code{history = TRUE}) A list of W estimates for each step of ADMM over the pathwise solutions.} \item{history_Beta_by_lam}{(if \code{history = TRUE}) A list of Beta estimates for each step of ADMM over the pathwise solutions.}
#' @export
#'
noteargis = function(data_input, data_type = c("c","m"),
                     n_resps_by_node=NULL, intcpt = NULL,
                     lambdas = NULL, n_lams=NULL, admm_args = admm_arg_ctrl(),
                     add_stop_rule=TRUE, fac_grp_lasso=FALSE,
                     verbose=FALSE, fit_hist=FALSE){

  n_obs = dim(data_input)[1]

  if (data_type == "m") {

    if (is.null(intcpt)) intcpt = "always"

    if (intcpt == "none"){
      message("Multinomial networks does not support 'intcpt = none' option.")
      message("'intcpt = none' will be coerced to 'always'.")
      intcpt = "always"
    }

    data_info = multi.data_info(data_input, n_resps_by_node)


  } else if (data_type == "c") {

    if (is.null(intcpt)) intcpt = "none"

    data_info = conti.data_info(data_input, n_resps_by_node, intcpt)

  }

  ## Set tuning parameters ####
  if (is.null(lambdas)){
    if (is.null(n_lams)){
      n_lams = 30
      lambdas = gen_lambdas(data_info, n_lams)
    } else {
      lambdas = gen_lambdas(data_info, n_lams)
    }
  } else {
    if (is.null(n_lams)) {
      n_lams = length(lambdas)
    } else {
      if (length(lambdas) != n_lams)
        n_lams = length(lambdas)
    }
  }


  ## Initialize parameters ####
  init_W_info = initialize_W(data_info, intcpt = intcpt)
  W_init = init_W_info$W_init
  white_coefs = init_W_info$white_coefs

  if (intcpt == "none"){
    Beta_init = Mu = matrix(0.0, data_info$n_expls, data_info$n_resps)
  } else {
    Beta_init = Mu = matrix(0.0, data_info$n_expls + 1, data_info$n_resps)
  }


  ## Set W_update according to types_by_node ####
  if (data_type == "m") {
    W_update = multi.W_update
  } else if (data_type == "c"){
    W_update = conti.W_update
  }


  ## Set penalizing method
  Beta_update = grp_lasso.Beta_update
  #... Later, can be modified by allowing various penalizing


  ## ADMM loop ####
  admm_results = admm_loop(W_update, Beta_update, W_init, Beta_init, white_coefs,
                           data_info, lambdas, admm_args,
                           add_stop_rule, fac_grp_lasso,
                           verbose, fit_hist)

  Beta_est_by_lam = admm_results$Beta_est_by_lam

  A_est_by_lam = push_dag(Beta_est_by_lam, data_info, lambdas)

  fit = c(list(A_est_by_lam = A_est_by_lam), admm_results)
  class(fit) <- "noteargis"
  return(fit)


} # end of function: noteargis()
