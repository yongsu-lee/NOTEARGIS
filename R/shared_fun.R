#' Control ADMM Arguments
#'
#' @param max_iter a positive integer, the number of loops in ADMM algorithm. Default is \code{100}.
#' @param rho a positive double scalar, a penalized parameter in the ADMM loop. Default is \code{1}.
#' @param kappa a postive double scalar, a relaxation pameter in the ADMM loop for the better convergnece. Recommended value is between 1.6 and 1.8(Default).
#' @param abs_tol a positive double scale, an absolute tolerance for stopping criterion in the ADMM loop.
#' @param rel_tol positive double scale, an relative tolerance for stopping criterion in the ADMM loop.
#' @param eps a positive scalar, a tolerance of estimates of parameter for the various purpose. Default is \code{1e-8}
#' @param inner_verbose blah
#'
#' @return a list object with arguments for customizing ADMM algorithm.
#' @export
#'
admm_arg_ctrl <- function(max_iter = 100, rho = 1.0, kappa = 1.8,
                          abs_tol = 1e-4, rel_tol = 1e-2, eps = 1e-8, inner_verbose = F){

  value = list(max_iter = max_iter, rho = rho, kappa = kappa, abs_tol = abs_tol,
               rel_tol = rel_tol, eps = eps, inner_verbose = inner_verbose)

  invisible(value)

}

#' Draw a DAG plot
#'
#' @param graph_true blah
#' @param data_type blah
#' @param font_size blah
#' @param main blah
#' @param cex.main blah
#' @param default_nodes_name blah
#' @param ... blah
#'
#' @return blah
#' @importFrom graphics title
#' @export
#'
dag_plot = function(graph_true, data_type, font_size = 1, main, cex.main = 2,
                    default_nodes_name = F, ...){

  if (missing(main)) main = "Graph"

  n_nodes = length(V(graph_true))
  V(graph_true)$color<- ifelse(data_type == "m", "black", "white")

  if (default_nodes_name == T) {
    nodes_name = V(graph_true)$name
  } else {
    nodes_name = 1:n_nodes
  }


  plot(graph_true,
       margin = c(0,0,0,0),
       vertex.label = nodes_name,
       vertex.label.color = ifelse(data_type == "m", "white", "black"),
       vertex.label.cex = font_size,
       vertex.label.font = 2,
       ...)

  title(main = main, cex.main = cex.main)

}

gen_node_types = function(n_conti_nodes, n_multi_nodes, seed = 1){
  n_nodes = n_conti_nodes + n_multi_nodes
  types_by_node = rep("m", n_nodes)

  set.seed(seed)
  sample(n_nodes, n_conti_nodes)
  types_by_node[sample(n_nodes, n_conti_nodes)] <- "c"
  names(types_by_node) = paste0("Z",1:n_nodes)
  return(types_by_node)

}

#' Random level number generate by node
#'
#' @param types_by_node blah
#' @param max_n_levels blah
#' @param seed blah
#'
#' @return blah
#' @export
#'
gen_node_levels = function(types_by_node, max_n_levels = 2, seed = 1){

  n_nodes = length(types_by_node)
  n_levels_by_node = rep(1, n_nodes)
  idx_multi_nodes = (types_by_node == "m")

  set.seed(seed)
  n_levels = sample(x = 1:(max_n_levels-1),
                    size = sum(idx_multi_nodes), replace = TRUE) + 1

  n_levels_by_node[idx_multi_nodes] <- n_levels
  names(n_levels_by_node) <- paste0("V",1:n_nodes)

  return(n_levels_by_node)

}

## For messages on algorithm iterations
ord = function(x) {
  ifelse(x == 1, "st", ifelse(x==2, "nd","th"))
}
#... appear in: admm_loop(); *.W_update()


gen_indices_by_node = function(sizes_by_node){

  idx = split(1:sum(sizes_by_node),
              rep(seq_along(sizes_by_node), times = sizes_by_node))

  return(idx)

}
#... (OLD version): gen_x_idx(); gen_y_idx()

X_names = function(n_nodes, types_by_node, n_expls_by_node){

  var_name_temp = rep(paste0("X",1:n_nodes), times= n_expls_by_node)
  elmts_temp = as.list(n_expls_by_node)

  dummys_temp = as.list(rep(0,n_nodes))

  for (j in 1:n_nodes){

    if (types_by_node[j] == "m"){
      dummys_temp[[j]] <- paste0(":D", 1:n_expls_by_node[j])
    } else if (types_by_node[j] == "c") {
      dummys_temp[[j]] <- paste0(":", 1:n_expls_by_node[j])
    }
  }

  dummys_temp[types_by_node=="c" & n_expls_by_node==1] <- ""
  dummys = unlist(dummys_temp)

  var_name = paste0(var_name_temp, dummys)

  return(var_name)
}

Y_names = function(n_nodes, types_by_node, n_resps_by_node){

  var_name_temp = rep(paste0("Y",1:n_nodes), times= n_resps_by_node)
  elmts_temp = as.list(n_resps_by_node)

  levels_temp = as.list(rep(0,n_nodes))

  for (j in 1:n_nodes){

    if (types_by_node[j] == "m"){
      levels_temp[[j]] <- paste0("=", 1:n_resps_by_node[j])
    } else if (types_by_node[j] == "c") {
      levels_temp[[j]] <- paste0(":", 1:n_resps_by_node[j])
    }
  }

  levels_temp[types_by_node=="c" & n_resps_by_node==1] <- ""
  levels = unlist(levels_temp)

  var_name = paste0(var_name_temp, levels)

  return(var_name)

}

gen_X1_design= function(data_input, types_by_node, n_expls_by_node){

  n = dim(data_input)[1]
  n_nodes = length(types_by_node)
  idx_expls_by_node = gen_indices_by_node(n_expls_by_node)

  X1 = matrix(NA, n, 0)

  for (j in 1:n_nodes){ # j = 1

    if (types_by_node[j] == "c" & n_expls_by_node[j] > 1){
      x_temp = data_input[, idx_expls_by_node[[j]]]
    } else {
      x_temp = data_input[,j]
    }

    if (types_by_node[j] == "m"){ # for multi-level node
      X1 = cbind(X1, model.matrix(~x_temp)[,-1])
    } else { # for gaussian node
      X1 = cbind(X1, x_temp)
    }
  }

  colnames(X1) <- X_names(n_nodes, types_by_node, n_expls_by_node)

  return(X1)
}

gen_Y_design= function(data_input, types_by_node, n_resps_by_node){

  n = dim(data_input)[1]
  n_nodes = length(types_by_node)
  idx_resps_by_node = gen_indices_by_node(n_resps_by_node)

  Y = matrix(NA, n, 0)

  for (j in 1:n_nodes){ # j = 1

    if (types_by_node[j] == "c" & n_resps_by_node[j] > 1){
      y_temp = data_input[, idx_resps_by_node[[j]]]
    } else {
      y_temp = data_input[,j]
    }


    if (types_by_node[j] == "m"){ # for multi-level node
      nc = nlevels(y_temp)
      Y = cbind(Y, diag(nc)[as.numeric(y_temp),])
    } else { # for gaussian node
      Y = cbind(Y, y_temp)
    }
  }

  colnames(Y) <- Y_names(n_nodes, types_by_node, n_resps_by_node)

  return(Y)

}

vectorizing = function(V, idx_intcpt){
  intcpts = cbind(V[idx_intcpt,])
  coefs = matrix(V[-idx_intcpt,], ncol = 1)
  v = rbind(coefs, intcpts)
  return(v)
}

matricizing = function(w, n_expls, n_resps, coef_only = F){
  n_coefs = n_expls * n_resps
  coefs = c(head(w, n = n_coefs))
  intcpts = c(tail(w, n = n_resps))

  W1 = matrix(coefs, nrow = n_expls, ncol = n_resps)

  if (coef_only == T) return(W1) else {
    W = rbind(W1, intcpts)
    rownames(W) <- NULL
    return(W)
  }
}


norm2 = function(x) {norm(as.matrix(x), "f")}

grp_shrinkage = function(x, a){
  z = max(0, 1- a/norm2(x)) * x
  return(z)
}

calc_A_est = function(Beta_new, data_info, intcpt = "always", eps = 1e-8){


  Se = data_info$Se
  Sr = data_info$Sr

  if (intcpt == "none"){
    Beta1_new = Beta_new
  } else {
    idx_intcpt = data_info$n_expls + 1
    Beta1_new = Beta_new[-idx_intcpt, ]
  }

  K = Se %*% Beta1_new^2 %*% t(Sr)
  A_est = (K > eps)*1

  return(A_est)
}

#' Print a \code{noteargis} object
#'
#' @param x fitted \code{noteargis} object
#' @param \dots additional arguments if available
#'
#' @return Display \code{graph_est_by_lam} and \code{A_est_by_lam}
#' @method print noteargis
#' @export
print.noteargis = function(x, ...){
  out = x
  cat("Number of Estimated Edges by Pathwise Solutions: \n")
  cat(attr(out$A_est_by_lam, "n_edges"))
  cat("\n\n")
  cat("Tuning Parameters: \n")
  cat(round(out$lambdas,4))
  cat("\n\n")
  cat("Estimated Adjacency Matrix by Pathwise Solutions: \n\n")
  print(out$A_est_by_lam)
}


skel = function(A){
  result_temp = A + t(A)
  A_skeleton = (result_temp != 0) * 1
  return(A_skeleton)
}


rev_edge = function(A_true, A_est){

  cpdag_A = pcalg::dag2cpdag(A_true)*1
  cpdag_est_A = pcalg::dag2cpdag(A_est)*1

  # find out undirected edges from true/estiamte graphs
  comp1 = cpdag_A == t(cpdag_A)
  comp2 = cpdag_est_A == t(cpdag_est_A)

  rev_edges_temp = ((A_est==1) & (t(A_true)==1))

  # Among the rev_edges_temp, we will exclude case that
  # ... corresponding edge is undirected in both
  # ... CPDAGs of estimated and true graph
  non_R = sum((rev_edges_temp & comp1) & comp2)
  return(sum(rev_edges_temp) - non_R)
}

