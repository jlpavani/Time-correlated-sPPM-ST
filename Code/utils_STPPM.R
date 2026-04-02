# *************************************************************************** #
# Packages
# *************************************************************************** #
list.of.packages <- c("igraph","GIGrvg")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
library("igraph","GIGrvg")

# *************************************************************************** #
# Functions
# *************************************************************************** #
tree_up = function(adj_mat, cl_lab){
  
  my_graph = igraph::graph_from_adjacency_matrix(adj_mat, mode = "undirected")
  adj_list = igraph::get.edgelist(my_graph) 
  nedge = nrow(adj_list) 
  
  omega = rep(NA, nedge) 
  
  for(i in 1:nedge){
    
    edge = adj_list[i, ]
    v1 = edge[1]; v2 = edge[2]
    
    if(cl_lab[v1] == cl_lab[v2]){
      omega[i]= runif(1, 0, 1)
    }else{
      omega[i] = runif(1, 10, 20)
    }
  }
  
  E(my_graph)$weight = omega
  
  return(igraph::mst(my_graph))
}

get_log_ratio = function(v1, v2, edge_status, part, ncl_t, narea, y_t, OXB_t, z_t, shape_theta, rate_theta, upsilon, kappaa,
                         uu, cc, tt, qq){ 
  
  # ******************************** STEP1: Identify all nodes that belongs to the same group of u and v
  left_component = igraph::bfs(part, root = v1, restricted = (1:narea)[-v2], unreachable = FALSE)$order
  left_component = left_component[!is.na(left_component)]
  
  right_component = igraph::bfs(part, root = v2, restricted = (1:narea)[-v1], unreachable = FALSE)$order
  right_component = right_component[!is.na(right_component)]
  
  # ******************************** STEP2: Calculate quantities for the two components and for the whole cluster 
  sum_dta1 = sum(y_t[left_component])
  sum_dta2 = sum(y_t[right_component])
  sum_dta = sum_dta1 + sum_dta2
  
  sum_off1 = sum(OXB_t[left_component] * z_t[left_component])
  sum_off2 = sum(OXB_t[right_component] * z_t[right_component])
  sum_off = sum_off1 + sum_off2
  
  # ******************************** STEP3: Calculate the log-ratio
  Cprime = ifelse(edge_status, ncl_t, ncl_t - 1)
  pos = tt:(tt-qq)
  
  log_const_beta = log(narea - Cprime + kappaa + sum((cc[pos[pos>0]] - uu[pos[pos>0]]), na.rm = T) - 1) - 
    log(Cprime + upsilon + sum(uu[pos[pos>0]], na.rm = T) - 1)
  
  log_const_gamma1 = lgamma(shape_theta) - shape_theta * log(rate_theta) + 
    lgamma(sum_dta + shape_theta) - lgamma(sum_dta1 + shape_theta) - lgamma(sum_dta2 + shape_theta)
  
  log_const_gamma2 = - (sum_dta + shape_theta) * log(sum_off + rate_theta) +
    (sum_dta1 + shape_theta) * log(sum_off1 + rate_theta) +
    (sum_dta2 + shape_theta) * log(sum_off2 + rate_theta) 
  
  return(log_const_beta + log_const_gamma1 + log_const_gamma2)
}

partition_up = function(adj_mat, narea, ncl_t, cl_lab, y_t, OXB_t, z_t,shape_theta, rate_theta, upsilon, kappaa, 
                        uu, cc, tt, qq){
  
  # ******************************** STEP1: Partition is equal to the tree
  part = tree_up(adj_mat, cl_lab)
  edges = igraph::get.edgelist(part)
  edges = edges[sample(narea - 1),]
  
  # ******************************** STEP2: Binary representation is equal to previous partition
  was_there = apply(edges, 1, function(edge){
    
    v1 = edge[1]; v2 = edge[2]
    
    if(cl_lab[v1] != cl_lab[v2]){
      e = paste(v1, "|", v2, sep = '')
      part <<- igraph::delete_edges(part, edges = e)
      FALSE
    }else{
      TRUE
    }
  })
  
  # ******************************** STEP3: Decision making - accept-reject criterion
  for(idx in 1:(narea - 1)){
    
    v1 = edges[idx, 1]
    v2 = edges[idx, 2]
    edge_status = was_there[idx]
    
    # ****** STEP3.1: log-ratio calculation
    log_ratio = get_log_ratio(v1, v2, edge_status, part, ncl_t, narea, y_t, OXB_t, z_t, shape_theta, rate_theta, 
                              upsilon, kappaa, uu, cc, tt, qq)
    
    # ****** STEP3.2: log-score calculation
    unif_value = runif(1, 0, 1)
    log_score = log(unif_value/(1 - unif_value))
    
    # ****** STEP3.3: decision
    if(log_ratio < log_score){
      if(edge_status == TRUE){ 
        e = paste (v1, "|" , v2, sep = '')
        part = igraph::delete_edges(part, edges = e)
      }
    }else{
      if(edge_status == FALSE){
        part = igraph::add_edges(part, c(v1, v2))
      }
    }
    ncl_t = igraph::clusters(part)$no
  }
  
  # ******************************** STEP4: Partition update
  res_cl_lab = igraph::clusters(part)
  
  return(list(cl_id = res_cl_lab$membership, ncl = res_cl_lab$no))
}

partition_up2 = function(adj_mat, narea, ncl_t, cl_lab, rho_s){
  
  # ******************************** STEP1: Partition is equal to the tree
  part = tree_up(adj_mat, cl_lab)
  edges = igraph::get.edgelist(part)
  edges = edges[sample(narea - 1),]
  
  # ******************************** STEP2: Binary representation is equal to previous partition
  was_there = apply(edges, 1, function(edge){
    
    v1 = edge[1]; v2 = edge[2]
    
    if(cl_lab[v1] != cl_lab[v2]){
      e = paste(v1, "|", v2, sep = '')
      part <<- igraph::delete_edges(part, edges = e)
      FALSE
    }else{
      TRUE
    }
  })
  
  # ******************************** STEP3: Decision making - accept-reject criterion
  for(idx in 1:(narea - 1)){
    
    v1 = edges[idx, 1]
    v2 = edges[idx, 2]
    edge_status = was_there[idx]
    
    # ****** STEP3.1: log-ratio calculation
    log_ratio = log(1 - rho_s) - log(rho_s)
    
    # ****** STEP3.2: log-score calculation
    unif_value = runif(1, 0, 1)
    log_score = log(unif_value/(1 - unif_value))
    
    # ****** STEP3.3: decision
    if(log_ratio < log_score){
      if(edge_status == TRUE){
        e = paste (v1, "|" , v2, sep = '')
        part = igraph::delete_edges(part, edges = e)
      }
    }else{
      if(edge_status == FALSE){
        part = igraph::add_edges(part, c(v1, v2))
      }
    }
    ncl_t = igraph::clusters(part)$no
  }
  
  # ******************************** STEP4: Partition update
  res_cl_lab = igraph::clusters(part)
  
  return(list(cl_id = res_cl_lab$membership, ncl = res_cl_lab$no))
}