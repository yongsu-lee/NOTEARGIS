#' Title
#'
#' @param data_input blah
#' @param data_type blah
#' @param n_resps_by_node blah
#' @param intcpt blah
#' @param lambdas blah
#' @param n_lams blah
#' @param admm_args blah
#' @param add_stop_rule blah
#' @param fac_grp_lasso blah
#' @param verbose blah
#' @param fit_hist blah
#'
#' @return blah blah
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
