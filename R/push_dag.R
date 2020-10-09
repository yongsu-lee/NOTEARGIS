push_dag = function(Beta_est_by_lam, data_info, lambdas, eps = 1e-8,
                    method = c("no_path_constr", "path_constr", "hardthre")){


  if (missing(method)) method = "no_path_constr"
  method = switch(method, no_path_constr = "n", path_constr = "p", hardthre = "h")

  n_nodes = data_info$n_nodes
  intcpt = data_info$intcpt

  if (intcpt != "none") idx_intcpt = data_info$n_expls + 1

  Sr = data_info$Sr
  Se = data_info$Se

  n_lams = length(Beta_est_by_lam)

  A_est_by_lam = as.list(rep(0, n_lams))
  # graph_est_by_lam = as.list(rep(0, n_lams))

  n_edges_by_lam = c()
  modified_dag_by_lam = c()

  dagness_by_lam = c()
  init_dagness_by_lam = c()

  for (ell in 1:n_lams){ # ell = 1

    Beta_est =  Beta_est_by_lam[[ell]]

    if (intcpt == "none") {
      Beta1_est = Beta_est
    } else {
      Beta1_est = Beta_est[-idx_intcpt, ]
    }

    K = Se %*% Beta1_est^2 %*% t(Sr)
    nz_coef = sort(K[K > eps])
    A_temp =  (K > eps)*1

    graph_temp = graph_from_adjacency_matrix(A_temp)
    isdag = is_dag(graph_temp)
    init_dagness_by_lam[[ell]] = isdag

    if ( method != "h"){

      i = 0
      while (!isdag){

        A_temp_prev = A_temp

        i = i + 1
        A_temp = (K > nz_coef[i])*1

        graph_temp = graph_from_adjacency_matrix(A_temp)
        isdag = is_dag(graph_temp)

        if (ell > 1 & method == "p") {
          if (n_edges_by_lam[ell-1] >= sum(A_temp == 1)) {
            A_temp = A_temp_prev
            break
          }
        }

      } # end of while loop (!isdag)

    }

    n_edges_by_lam[ell] = sum(A_temp == 1) # to be attr
    # dagness_by_lam[ell] = isdag # to be attr
    # graph_est_by_lam[[ell]] = graph_temp
    colnames(A_temp) = paste0("Z",1:n_nodes)
    rownames(A_temp) = paste0("Z",1:n_nodes)
    A_est_by_lam[[ell]] = A_temp
    attr(A_est_by_lam[[ell]], "is_dag") = isdag
    attr(A_est_by_lam[[ell]], "n_edges") = sum(A_temp == 1)
    attr(A_est_by_lam[[ell]], "lambda") = lambdas[ell]
    attr(A_est_by_lam[[ell]], "graph_est") = graph_from_adjacency_matrix(A_temp)


  } # end of for loop ell = 1:n_lams

  # attr(A_est_by_lam,"is_dag") = dagness_by_lam
  attr(A_est_by_lam,"n_edges") = n_edges_by_lam

  return(A_est_by_lam)

} # end of function
