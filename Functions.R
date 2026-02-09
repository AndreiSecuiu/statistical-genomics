########## Description ##########

#' This file contains all functions which implement structure learning for the 
#' project. A distinction is made between functions written by myself and those
#' provided by the professor

#' Required packages to run this file: stats, ROC, limma

########## My functions ##########

# ------ Analyze benefit of repeated restarts ------ 

restarts_benefit_analysis = function(scores){ 
  #' Based on an i.i.d. sample of "scores" (corresponding to local optima of some
  #' search method), this function estimates: 
  #'   ~ the expected value of the maximum of "k" draws
  #'   ~ the expected improvement of the new maximum over the old one
  #' after the search procedure has been employed "k" times. "k" is ranged from 
  #' 1 to N-1 (N being the number of elements in the sample)
  #' 
  #' *Note*: The function has been written with the help of Google Gemini, though
  #' its implementation, logic and commenting has been improved manually
  #' 
  #' --- Inputs:
  #' 
  #' scores = array. It contains a sample of the scores for the search method, 
  #' which must be of a length at least K
  #' 
  #' --- Outputs: 
  #' 
  #' A data frame containing the desired outputs on each column. For a row k:
  #'   ~ exp_max: the average value of the maximum of "k" samples
  #'   ~ exp_imp: the difference between the average maximum of "k+1" draws and
  #'   the average maximum of "k" draws
  
  N = length(scores)
  
  # - Pre-compute empirical CDF of scores
  
  # Sort for efficient computation
  scores = sort(scores)    
  
  # Get rid of ties; consecutively equal values are assigned the highest rank. 
  unique_idx = unique(rank(scores, ties.method = "max"))
  unique_vals = scores[unique_idx]
  
  F_N = unique_idx / N                # F_N = rank_value / N = empirical_CDF(value)
  F_prev = c(0, F_N[-length(F_N)])    # F_prev empirical_CDF(value-)
  pmf = F_N - F_prev                  # pmf = P(X = val) = empirical_CDF(value) - empirical_CDF(value-)
  
  # Initialize output
  OUT = data.frame(
    exp_max = numeric(N-1),
    exp_imp = numeric(N-1)
  )
  
  # - Compute metrics for each K
  
  for (k in 1:(N-1)) {
    
    # i) Expected maximum (E[M_k])
    
    pmf_max = F_N ^ k - F_prev ^ k
    OUT$exp_max[k] = sum(unique_vals * pmf_max)
  }
  
  # ii) Expected improvement
  exp_max_end = sum(unique_vals * (F_N ^ N - F_prev ^ N))
  OUT$exp_imp = diff(c(OUT$exp_max, exp_max_end))
  
  return(OUT)
}

# ------ Gradient Ascent ------

#' From a given starting point, do a greedy search after moving via single edge
#' operations. Adapted from Marco's strMCMC

strGA = function(data, inc_mat_start, max_iter){ 
  #' This function performs gradient ascent starting from the graph with incidence
  #' matrix "inc_mat_start". It stops either when a (local) maximum is found, 
  #' or when a number of iterations "max_iter" is reached
  #'
  #' --- Output:
  #' 
  #' A list containing the following:
  #' 
  #' ~ A vector tracking the log_BGe scores for the explored best neighbour graphs
  #' ~ The incidence matrix of the best graph
  #' ~ The status of the search (TRUE if a local optimum is reached, or FALSE if
  #' max_iter is reached first)
  
  n <- nrow(data)                               # number of nodes
  m <- ncol(data)                               # number of observations

  log_score_trace = array(dim = max_iter)   # log scores of intermediate graphs
  
  # - Initialize values (prior, T_m)
  
  fan.in = n-1
  v=1
  mu=numeric(n)
  a=n+2 
  
  T_0 = diag(1,n,n)
  T_0 <- (a-n-1) * T_0   # scaling T_0
  
  T_m = T_0 + (m-1)* cov(t(data)) + ((v*m)/(v+m))* (mu - rowMeans(data))%*%t(mu - rowMeans(data))
  
  # - Computation of the log BGe Score of the FIRST graph
  
  P_local_num = numeric(n) 
  P_local_den = numeric(n) 
  
  out = COMPUTE_BGE(T_m,m,P_local_num,P_local_den,inc_mat_start,a,v,T_0)
  
  bge_old     = out[[1]]
  P_local_num = out[[2]]
  P_local_den = out[[3]]
  
  log_score_trace[1] = bge_old                     
  
  # initialize the best (yet) DAG
  inc_mat_best = inc_mat_start
  anc_mat_best = ancestor(inc_mat_best)
  
  # - Main search
  
  MAX_FOUND = FALSE
  i = 1
  
  while((MAX_FOUND == FALSE) && (i <= max_iter)){
    # - Enumerate all valid neighbours via single edge operations
    
    OUT = GET_ALL_NEIGHBOURS(inc_mat_best, anc_mat_best, fan.in)
    all_neigh_scores = array(dim = length(OUT))
    
    for(j in 1:length(OUT)){
      out = OUT[[j]]
      
      # Compute the local score for each graph
      
      temp = COMPUTE_BGE_LOCAL(out$operation, out$child, out$parent, T_m, m, 
                               P_local_num, P_local_den, out$incidence, a, v, T_0)
      all_neigh_scores[j] = temp[[1]]
    }
    
    # - Find which neighbour has the maximum score. Compare to the previous best
    
    max_idx = which.max(all_neigh_scores)
    
    if(all_neigh_scores[max_idx] > log_score_trace[i]){      # Improvement possible
      out = OUT[[max_idx]]
      temp = COMPUTE_BGE_LOCAL(out$operation, out$child, out$parent, T_m, m, 
                               P_local_num, P_local_den, out$incidence, a, v, T_0)
      
      inc_mat_best = out$incidence      # Update incidence matrix
      anc_mat_best = out$ancestor       # Update the ancestor matrix
      
      log_score_trace[i+1] = temp[[1]]  # Update score
      P_local_num = temp[[2]]           # Update P_local (num & den)
      P_local_den = temp[[3]]
    } else{                                                 # Terminate the search (local max reached)
      MAX_FOUND = TRUE
    }
    
    # Update iteration idx
    i = i + 1
  }
  
  return(list(log_bge_trace = log_score_trace[1:i-1],
              incidence_best = inc_mat_best,
              optimum_found = MAX_FOUND)
         )
}

# --- Helper functions ---

#' Given an incidence matrix 'incidence' and its ancestor matrix ancestMat, provide
#' a list of all neighbour graphs that are accessible via single edge operations

# operation=1 <=> edge reversal 
# i.e. the edge "parent->child" is reversed to "child <- parent"

# operation=2 <=> edge deletion 
# i.e. the edge "parent->child" is deleted

# operation=3 <=> edge addition 
# i.e. the edge "parent->child" is added

GET_ALL_NEIGHBOURS <- function(incidence, ancestMat, fan.in) {
  
  n <- nrow(incidence)
  OUT <- list()
  
  # - Identify valid operations 

  # i) Deletions
  idx_delete <- which(incidence == 1)
  
  # ii) Additions        1- E(i,j) - I(i,j) - A(i,j)
  inter_add <- which(matrix(rep(1, n * n), nrow = n) - diag(1, n, n) - incidence - ancestMat > 0)
  add_mask <- matrix(numeric(n * n), nrow = n)
  add_mask[inter_add] <- 1
  
  # Enforce fan-in, i.e. refuse to add edges to nodes that already have fan.in parents 
  add_mask[, which(colSums(incidence) > fan.in - 1)] <- 0
  idx_add <- which(add_mask == 1)
  
  # iii) Reversals       I - (I^t * A)^t
  inter_rev <- which(incidence - t(t(incidence) %*% ancestMat) == 1)
  re_mask <- matrix(numeric(n * n), nrow = n)
  re_mask[inter_rev] <- 1
  
  # Enforce fan-in, i.e. refuse to reverse an edge of all parents that have fan.in
  # grandparents (otherwise the former child becomes a new parent which breaks the
  # fan.in restriction)
  re_mask[which(colSums(incidence) > fan.in - 1), ] <- 0 
  idx_rev <- which(re_mask == 1)
  
  # - Iterate and construct neighbours 
  
  # i) Deletions 
  
  if(length(idx_delete) > 0){
    for(new_edge in idx_delete){
      
      # Get nodes
      coords <- arrayInd(new_edge, .dim = c(n, n))
      parent <- coords[1]
      child  <- coords[2]
      
      # Modify Incidence
      inc_new <- incidence
      inc_new[new_edge] <- 0
      
      # Update Ancestor 
      anc_new <- ancestMat
      anc_new[c(child, which(ancestMat[, child] == 1)), ] <- 0    # delete all ancestors of the child and its descendants
      
      top_name <- des_top_order(inc_new, ancestMat, child)        # rebuild ancestor connections
      
      for (d in top_name) {
        for (g in which(inc_new[, d] == 1)) {
          anc_new[d, c(g, (which(anc_new[g, ] == 1)))] <- 1
        }
      }
      
      # Store result
      OUT[[length(OUT) + 1]] <- list(
        incidence = inc_new,
        ancestor = anc_new,
        operation = 2,
        parent = parent,
        child = child
      )
    }
  }
  
  # ii) Additions 
  
  if (length(idx_add) > 0) {
    for (new_edge in idx_add) {
      
      # Get nodes
      coords <- arrayInd(new_edge, .dim = c(n, n))
      parent <- coords[1]
      child  <- coords[2]
      
      # Modify incidence
      inc_new <- incidence
      inc_new[new_edge] <- 1
      
      # Update ancestor 
      anc_new <- ancestMat
      anc_parent <- which(ancestMat[parent, ] == 1) # ancestors of the new parent
      des_child  <- which(ancestMat[, child] == 1)  # descendants of the child
      anc_new[c(child, des_child), c(parent, anc_parent)] <- 1
      
      # Store result
      OUT[[length(OUT) + 1]] <- list(
        incidence = inc_new,
        ancestor = anc_new,
        operation = 3,
        parent = parent,
        child = child
      )
    }
  }
  
  # iii) Reversals
  
  if (length(idx_rev) > 0) {
    for (new_edge in idx_rev) {
      
      # Get nodes 
      coords <- arrayInd(new_edge, .dim = c(n, n))
      parent <- coords[1]
      child  <- coords[2]
      
      # Modify incidence 
      inc_new <- incidence
      inc_new[parent, child] <- 0
      inc_new[child, parent] <- 1
      
      # Update ancestor
      
      # a) Deletion
      anc_new <- ancestMat
      anc_new[c(child, which(ancestMat[, child] == 1)), ] <- 0    # delete all ancestors of the child and its descendants
      top_name <- des_top_order(inc_new, ancestMat, child)
      
      for (d in top_name) {
        for (g in which(inc_new[, d] == 1)) {
          anc_new[d, c(g, (which(anc_new[g, ] == 1)))] <- 1
        }
      }
      
      # b) Addition
      anc_parent_idx <- which(anc_new[child, ] == 1)           # child becomes parent
      des_child_idx  <- which(anc_new[, parent] == 1)          # parent becomes child
      anc_new[c(parent, des_child_idx), c(child, anc_parent_idx)] <- 1
      
      # Store result
      OUT[[length(OUT) + 1]] <- list(
        incidence = inc_new,
        ancestor = anc_new,
        operation = 1,
        parent = parent,
        child = child
      )
    }
  }
  
  return(OUT)
}

#' Given the number of nodes *n*, generate the incidence matrix for a random DAG
#' The user has the option of incorporating a fan.in constraint & tweak how dense
#' the resulting network becomes via 'prob'

#' The generation mechanism works as follows:
#'  
#'   ~ given a topological order, generate the edges
#'     ~ start with a matrix of n^2 random values in [0,1]
#'     ~ to enforce acyclicity, make it upper triangular
#'     ~ keep only the edges whose values are < prob (so an edge appears with prob 'prob')
#'   ~ generate a random topological order
#'   ~ assign the node labels to the generated edges

GENERATE_RANDOM_DAG = function(n, prob = 2/n, fan.in = n) {
  
  # - Sample the edges
  
  random_vals = matrix(runif(n * n), nrow = n)
  U = matrix(0, nrow = n, ncol = n)
  U[upper.tri(U)] = as.numeric(random_vals[upper.tri(U)] < prob)
  
  # - Enforce fan.in 
  
  # Check which nodes have too many parents, i.e. the indices where the column sums
  # exceed fan.in . In all those places, keep precisely fan.in nodes (randomly chosen)
  
  bad_nodes = which(colSums(U) > fan.in)
  
  if (length(bad_nodes) > 0){
    for (node in bad_nodes) {
      parents = which(U[, node] == 1)
      
      # Sample fan.in parents
      parents_kept = sample(parents, fan.in)
      
      # Update the column
      U[, node] = 0
      U[parents_kept, node] = 1
    }
  }
  
  # - Sample a topological order
  
  topo_order = sample(1:n)
  
  # - Assign the node labels to the generated edges
  
  incidence = matrix(0, nrow = n, ncol = n)
  incidence[topo_order, topo_order] = U
  
  return(incidence)
}


# ------ Simulated Annealing ------

#' Implement a simple SA algorithm. Adapted from Marco's strMCMC

strSA = function(data, incidence, iterations, temperature_schedule){
  #' This function performs simulated annealing starting from the graph with incidence
  #' matrix "incidence". It stops when the number of iterations "iterations" is reached
  #' 
  #' *Note*: "temperature_schedule" must be an array of length "iterations" containing
  #' the temperature values to be used for each Metropolis-Hastings step
  #' 
  #' --- Output:
  #' 
  #' A list containing the following:
  #' 
  #' ~ A vector tracking the log_BGe scores of the MCMC graphs
  #' ~ The incidence matrix of the final graph

  n <- nrow(data)                               # number of nodes
  m <- ncol(data)                               # number of observations
  
  log_score_trace = array(dim = iterations)   # log scores of intermediate graphs
  
  # - Initialize values (prior, T_m)
  
  fan.in = n-1
  v=1
  mu=numeric(n)
  a=n+2 
  
  T_0 = diag(1,n,n)
  T_0 <- (a-n-1) * T_0   # scaling T_0
  
  T_m = T_0 + (m-1)* cov(t(data)) + ((v*m)/(v+m))* (mu - rowMeans(data))%*%t(mu - rowMeans(data))
  
  # - Computation of the log BGe Score of the FIRST graph
  
  P_local_num = numeric(n) 
  P_local_den = numeric(n) 
  
  out = COMPUTE_BGE(T_m,m,P_local_num,P_local_den,incidence,a,v,T_0)
  
  bge_old     = out[[1]]
  P_local_num = out[[2]]
  P_local_den = out[[3]]
  
  log_score_trace[1] = bge_old
  
  # first ancestor matrix
  ancestMat <- ancestor(incidence)
  
  # - The Metropolis-Hastings routine
  
  for (i in 2:iterations){
      out = PROPOSE_NEW_DAG(incidence,ancestMat,fan.in)
      
      incidence_new  = out[[1]]
      ancestMat_new  = out[[2]]
      p_forward      = out[[3]]
      p_backward     = out[[4]]
      operation      = out[[5]]
      child_node     = out[[6]]
      parent_node    = out[[7]]
      
      out = COMPUTE_BGE_LOCAL(operation,child_node,parent_node,T_m,m,P_local_num,
                              P_local_den,incidence_new,a,v,T_0)
      
      bge_new         = out[[1]]
      P_local_num_new = out[[2]]
      P_local_den_new = out[[3]]
      
      # Decide if to accept the proposed change
      
      HR = log(p_backward) - log(p_forward)      # Hastings ratio (from proposal probabilities)
      Delta_BGe = bge_new - bge_old              # Change in posterior ratio
      
      log_p_accept <- Delta_BGe/temperature_schedule[i] + HR
      
      rand <- log(runif(1))
      
      if(log_p_accept > rand){
        incidence    <- incidence_new
        ancestMat    <- ancestMat_new
        bge_old      <- bge_new
        P_local_num  <- P_local_num_new
        P_local_den  <- P_local_den_new
      }
    
      log_score_trace[i] = bge_old
  }
  
  return(list(log_bge_trace = log_score_trace,
              incidence_final = incidence)
         )
}

########### Data generation (RAF pathway) & model evaluation ##########

#' The functions below have been written and provided by prof. Marco Grzegorczyk

# ------ Make test data

make_test_Data <- function(m_obs, var_noise){
  #' m_obs: number of samples
  #' var_noise: variance of Gaussian distributed noise terms
  
  edges <- 20
  a1 <- runif(edges,0.5,2)
  a2 <- sample(c(-1,1),edges, replace=TRUE)
  a <- a1*a2    # vector with regression coefficients for the 20 edges
  
  # 1. pip3
  x_pip3 <- rnorm(m_obs, sd=1)
  pip3 <- (x_pip3 - mean(x_pip3))/sd(x_pip3)
  
  # 2. plcg
  x_plcg <- a[1]* pip3 + rnorm(m_obs, sd=sqrt(var_noise))
  plcg <- (x_plcg - mean(x_plcg))/sd(x_plcg)
  
  # 3. pip2
  x_pip2 <- a[2]* pip3 + a[3]*plcg + rnorm(m_obs, sd=sqrt(var_noise))
  pip2 <- (x_pip2 - mean(x_pip2))/sd(x_pip2)
  
  # 4. pkc
  x_pkc <- a[4]* pip2 + a[5]*plcg + rnorm(m_obs, sd=sqrt(var_noise))
  pkc  <- (x_pkc - mean(x_pkc))/sd(x_pkc)
  
  # 5. pka
  x_pka <- a[6]* pkc + rnorm(m_obs, sd=sqrt(var_noise))
  pka  <- (x_pka - mean(x_pka))/sd(x_pka)
  
  # 6. jnk
  x_jnk <- a[7]* pkc + a[8]* pka + rnorm(m_obs, sd=sqrt(var_noise))
  jnk  <- (x_jnk - mean(x_jnk))/sd(x_jnk)
  
  # 7. p38
  x_p38 <- a[9]* pkc + a[10]* pka + rnorm(m_obs, sd=sqrt(var_noise))
  p38  <- (x_p38 - mean(x_p38))/sd(x_p38)
  
  # 8. raf
  x_raf <- a[11]* pkc + a[12]* pka + rnorm(m_obs, sd=sqrt(var_noise))
  raf  <- (x_raf - mean(x_raf))/sd(x_raf)
  
  # 9. mek
  x_mek <- a[13]* pkc + a[14]* pka + a[15]* raf + rnorm(m_obs, sd=sqrt(var_noise))        
  mek  <- (x_mek - mean(x_mek))/sd(x_mek)
  
  # 10. erk
  x_erk <- a[16]* pka + a[17]* mek + rnorm(m_obs, sd=sqrt(var_noise))
  erk  <- (x_erk - mean(x_erk))/sd(x_erk)
  
  # 11. akt
  x_akt <- a[18]* pip3 + a[19]* pka + a[20]* erk + rnorm(m_obs, sd=sqrt(var_noise))
  akt  <- (x_akt - mean(x_akt))/sd(x_akt)    
  
  daten <- cbind(pip3, plcg, pip2, pkc, pka, jnk, p38, raf, mek, erk, akt)
  
  return(t(daten))
}

# ------ Create the true RAF network incidence matrix

make_true_Net <- function(){
  
  NETWORK <- matrix(numeric(11*11),11,11)
  
  # 1. pip3
  
  # 2. plcg
  NETWORK[1,2] <- 1
  
  # 3. pip2
  NETWORK[1,3] <- 1
  NETWORK[2,3] <- 1
  
  # 4. pkc
  NETWORK[2,4] <- 1
  NETWORK[3,4] <- 1
  
  # 5. pka
  NETWORK[4,5] <- 1
  
  # 6. jnk
  NETWORK[4,6]  <- 1
  NETWORK[5,6]  <- 1
  
  # 7. p3B
  NETWORK[4,7]  <- 1
  NETWORK[5,7]  <- 1
  
  # 8. raf
  NETWORK[4,8]  <- 1
  NETWORK[5,8]  <- 1
  
  # 9. mek
  NETWORK[4,9]  <- 1
  NETWORK[5,9]  <- 1
  NETWORK[8,9]  <- 1
  
  # 10. erk
  NETWORK[5,10]  <- 1
  NETWORK[9,10]  <- 1
  
  # 11. akt
  NETWORK[1,11]  <- 1
  NETWORK[5,11]  <- 1
  NETWORK[10,11]  <- 1
  
  return(NETWORK)
}

# ------ Function for plotting the ROC curve, computing AUROC

draw_ROC <- function(postP, trueEdges, steps=0.01){
  
  n1 <- numeric(0)
  for(i in 1:dim(postP)[1]){
    for(j in 1:dim(postP)[2]){
      if (abs(i-j)>0){
        n1 <- c(n1, postP[i,j])
      }
    }
  }
  
  n2 <- numeric(0)
  for(i in 1:dim(trueEdges)[1]){
    for(j in 1:dim(trueEdges)[2]){
      if (abs(i-j)>0){     
        n2 <- c(n2, trueEdges[i,j])
      }
    }
  }
  
  sp <- sort(n1)         # posterior probabilities
  st <- n2[order(n1)]    # true ordered by posterior probabilities
  
  rocc.obj <- rocdemo.sca(n2, n1, rule=NULL, cutpts=NA)#seq(0,1,steps))
  
  return(list(plot((1-rocc.obj@"spec"), rocc.obj@"sens", type="l",las=1, xlab="1 - specificity", ylab="sensitivity", main='ROC curve'), abline(0,1, lty=2)))
  
  
}

compute_AUROC <- function(postP, trueEdges, steps=0.01){
  
  n1 <- numeric(0)
  for(i in 1:dim(postP)[1]){
    for(j in 1:dim(postP)[2]){
      if (abs(i-j)>0){
        n1 <- c(n1, postP[i,j])
      }
    }
  }
  
  n2 <- numeric(0)
  for(i in 1:dim(trueEdges)[1]){
    for(j in 1:dim(trueEdges)[2]){
      if (abs(i-j)>0){     
        n2 <- c(n2, trueEdges[i,j])
      }
    }
  }
  
  
  sp <- sort(n1)         # posterior probabilities
  st <- n2[order(n1)]    # true ordered by posterior probabilities
  
  rocc.obj <- rocdemo.sca(n2, n1, rule=NULL, cutpts=NA)#seq(0,1,steps))
  
  return(auROC(st,sp))
  
}

########## Marco's functions ##########

#' The functions below have been written and provided by prof. Marco Grzegorczyk

#' ------ Topological order of nodes in a DAG

top_order <- function(incidence){
  
  n <- nrow(incidence)
  Order <- numeric(n)
  fan_in <- numeric(n)
  no_fan_in <- numeric(0)
  m <- 1
  
  for (p in 1:n){                                       # number of parent nodes at the beginning
    fan_in[p] <- sum(incidence[,p])
  }
  
  no_fan_in <- which(fan_in==0)
  
  while (length(which(Order==0))>0){                    # as long as there is a node without an order
    fan_in[which(incidence[no_fan_in[1],]==1)] <- fan_in[which(incidence[no_fan_in[1],]==1)] - 1
    no_fan_in <- c(no_fan_in, c(which(incidence[no_fan_in[1],]==1),which(fan_in==0))[duplicated(c(which(incidence[no_fan_in[1],]==1),which(fan_in==0)))])
    Order[m] <- no_fan_in[1]
    no_fan_in <- no_fan_in[-1]
    m <- m+1
  }
  
  return(Order)
}

#' ------ Order the edges according to the algorithm required for extracting the CPDAG

#' The output is a matrix with 4 columns. Each row corresponds to an edge: 
#'  - "edges" (which entries in the incidence matrix has an edge (numbered from 1
#'  to n^2 from left to right)
#'  - "parents": the parent of the considered edge
#'  - "children": the child of the considered edge
#'  - "order": the order of the considered edge

order_edges <- function(incidence){  
  
  top.order <- top_order(incidence)
  n <- length(top.order)
  edges <- which(incidence!=0)
  children <- child(edges,n)
  parents <- parent(edges,n)
  m <- length(edges)
  ordered_edges  <- numeric(m)
  incidence_n <- incidence
  tog <- matrix(c(edges,parents,children,ordered_edges),ncol=4, byrow=FALSE)  
  k <- 1
  while(any(tog[,4]==0)){
    node1 <- top.order[which(colSums(incidence_n[,top.order])>0)][1]    # first node in top. order that has at least one parent
    par1<- tog[which(tog[,3]==node1),2]                # find the parents of  first child in the top. order that has an unordered edge incident into it
    g <- par1[which(par1>0)]
    f1 <- numeric(length(g))
    for (i in 1:length(g)){
      f1[i] <- which(top.order==g[i])
    }
    par2 <- g[which.max(f1)]                           # find the highest ordered node that has an edge leading into node1
    tog[which(tog[,2]==par2 & tog[,3]==node1),4] <- k     
    k <- k + 1
    incidence_n[tog[which(tog[,2]==par2 & tog[,3]==node1),1]] <- 0     # delete the edge in the "incidence" matrix
    tog[which(tog[,2]==par2 & tog[,3]==node1),2] <- 0
  }
  to <- matrix(c(edges,parents,children,tog[,4]),ncol=4,byrow=FALSE)
  return(to)                                          # return the whole matrix, the order is the fourth column
}

# --- Helper functions:

# input: the numbers of the edges in the incidence matrix and the number of nodes
child <- function(edges,n){ 
  p <- ceiling(edges/n)
  return(p)
}
parent <- function(edges,n){
  ch <- edges + n - child(edges,n)*n
  return(ch)
}

# ------  strMCMC and its helper functions

#' Take Marco's *corrected* structure MCMC function, then do one improvement
#' for speed: pre-allocate the arrays of final length in the memory

strMCMC = function(data,incidence,iterations,step_save){
  
  n <- nrow(data)    # number of nodes
  m <- ncol(data)    # number of observations
  
  fan.in = n-1
  
  v=1
  
  mu=numeric(n)
  
  a=n+2 
  
  T_0 = diag(1,n,n)
  
  T_0 <- (a-n-1) * T_0 # scaling T_0
  
  # Preallocate L1 & L2 for efficiency
  
  N_out = 1 + floor(iterations/step_save) # Number of final MCMC iterations
  
  L1 = array(dim = c(n, n, N_out))       # DAGs
  L2 = array(dim = N_out)                # log scores
  
  ################################################################################
  ### computation of the log BGe Score of the FIRST graph
  
  T_m = T_0 + (m-1)* cov(t(data)) + ((v*m)/(v+m))* (mu - rowMeans(data))%*%t(mu - rowMeans(data))
  
  P_local_num = numeric(n) 
  P_local_den = numeric(n) 
  
  out = COMPUTE_BGE(T_m,m,P_local_num,P_local_den,incidence,a,v,T_0)
  
  bge_old     = out[[1]]
  P_local_num = out[[2]]
  P_local_den = out[[3]]
  
  L1[,,1] = incidence                    
  L2[1] = bge_old                     
  
  # first ancestor matrix
  ancestMat <- ancestor(incidence)
  
  ####################################################################################################################################
  
  for (z in 2:N_out){
    
    for (count in 1:step_save){
      
      out = PROPOSE_NEW_DAG(incidence,ancestMat,fan.in)
      
      incidence_new  = out[[1]]
      ancestMat_new  = out[[2]]
      p_forward      = out[[3]]
      p_backward     = out[[4]]
      operation      = out[[5]]
      child_node     = out[[6]]
      parent_node    = out[[7]]
      
      out = COMPUTE_BGE_LOCAL(operation,child_node,parent_node,T_m,m,P_local_num,P_local_den,incidence_new,a,v,T_0)
      
      bge_new         = out[[1]]
      P_local_num_new = out[[2]]
      P_local_den_new = out[[3]]
      
      acceptance <- exp((bge_new +log(p_backward)) - (bge_old +log(p_forward)))
      
      rand <- runif(1)
      
      if(acceptance > rand){
        incidence    <- incidence_new
        ancestMat    <- ancestMat_new
        bge_old      <- bge_new
        P_local_num  <- P_local_num_new
        P_local_den  <- P_local_den_new
      }
    }
    
    L1[,,z] = incidence
    L2[z] = bge_old
  }
  return(list(L1,L2))
}

cpdag <- function(incidence){
  #' The output is a matrix with 5 columns. Each row corresponds to an edge: 
  #'  - "edges" (which entries in the incidence matrix has an edge (numbered from 1
  #'  to n^2 from left to right)
  #'  - "parents": the parent of the considered edge
  #'  - "children": the child of the considered edge
  #'  - "order": the order of the considered edge
  #'  - "status": 1 = compelled, -1 = reversible
  
  z <- order_edges(incidence)
  new_mat <- cbind(z,numeric(nrow(z)))    # edges, parents, children, order, zeros
  n_mat <- new_mat[order(new_mat[,4]),]   # sort the edges by its order
  vec <- numeric(nrow(z))
  while(any(vec==0)){                                  # while there are unlabeled edges            l.3
    if (length(vec)>1){                                  # if there are at least 2 edges
      first <- which(n_mat[,5]==0)[1]                    # first EDGE that ist labeled "unknown" (0)  l.4
      parent1 <- n_mat[first,2]                          # x   parent NODE
      child1 <- n_mat[first,3]                           # y   child NODE
      comp1 <- n_mat[which(n_mat[,3]==parent1 & n_mat[,5]==1),2]      # w NODES that have an edge incident into the parent labeled compelled)
    }
    if (length(vec)==1){
      first <- which(n_mat[5]==0)                      # first edge that ist labeled "unknown" (0)
      parent1 <- n_mat[2]                             # x   parent
      child1 <- n_mat[3]                              # y   child
      comp1 <- numeric(0)
    }
    for (j in comp1){                                   #                                            l.5
      if (incidence[j,child1]==0){                     # if w is not a parent of the child          l.6
        n_mat[first,5] <- 1                             # label x -> y compelled                     l.7
        n_mat[which(n_mat[,3]==child1),5] <- 1          # label every edge incident into y compelled l.7
        vec[first] <- 1
        vec[which(n_mat[,3]==child1)] <- 1
        break
      }
      if (incidence[j,child1]!=0)    {
        n_mat[which(n_mat[,2]==j & n_mat[,3]==child1),5] <- 1  # label w -> y compelled                l.10
        vec[which(n_mat[,2]==j & n_mat[,3]==child1)] <- 1
      }
    }
    if (length(vec)>1){
      if(n_mat[first,5]==0){
        
        moep <- n_mat[which(n_mat[,3]==child1 & n_mat[,2]!=parent1),2]      # other parents of the child
        if(length(moep)>0){                              #                     l.11
          for(o in moep){
            if(incidence[o,parent1]==0){
              vec[first] <- 1
              vec[which(n_mat[,3]==child1 & n_mat[,5]==0)] <- 1
              n_mat[first,5] <- 1                                     # label x -> y compelled
              n_mat[which(n_mat[,3]==child1 & n_mat[,5]==0),5] <- 1   # label all "unknown" edges incident into y compelled
              break
            }
            if(all(incidence[moep,parent1]!=0)){
              vec[first] <- -1
              vec[which(n_mat[,3]==child1 & n_mat[,5]==0)] <- -1
              n_mat[first,5] <- -1                                    # label x -> y reversible
              n_mat[which(n_mat[,3]==child1 & n_mat[,5]==0),5] <- -1  # label all "unknown" edges incident into y reversible
            }
          }
        }
        if(length(moep)==0){
          vec[first] <- -1
          vec[which(n_mat[,3]==child1 & n_mat[,5]==0)] <- -1
          n_mat[first,5] <- -1                                    # label x -> y reversible
          n_mat[which(n_mat[,3]==child1 & n_mat[,5]==0),5] <- -1  # label all "unknown" edges incident into y reversible
        }
      }
    }
    if (length(vec)==1){
      n_mat[5] <- -1                                    # label x -> y reversible
      vec <- -1
    }
  }
  return(n_mat)
}

ancestor <- function(incidence){
  #' calculation of the first ancestor matrix:
  
  incidence1 <- incidence
  incidence2 <- incidence
  k <- 1
  while (k < nrow(incidence)){
    incidence1 <- incidence1%*%incidence
    incidence2 <- incidence2 + incidence1
    k <-k+1
  }
  incidence2[which(incidence2[,]>0)] <- 1
  return(t(incidence2))}


c_function <- function(N,A){
  #' function for the computation of c(n, alpha)
  
  fact <- numeric(N)
  for (i in 1:N){
    fact[i] <- -lgamma((A+1-i)/2)
  }
  product <- sum(fact) -(A*N/2)*log(2)- (N*(N-1)/4)*log(pi)
  return(product)}


des_top_order <- function(incidence, ancest1,child){
  #' assign the topological order of the descendants of the child
  
  n <- nrow(incidence)
  top <- top_order(incidence)
  position_child <- which(top==child)
  top_all_after <- top[position_child:n]                # top. order without the "first" nodes
  desc <- which(ancest1[,child]==1)                     # descendants of the child
  inter_step <- c(child,desc,top_all_after)
  des_top <- inter_step[which(duplicated(inter_step))]
  return(des_top)
}

cpdag_list <- function(list.inc,E){    
  #' --- Input:
  #' 
  #' list.inc = an array of incidence matrices, as returned by strMCMC *in Andrei's*
  #' version
  #' E = index that marks the end of the burn-in phase
  
  # Initializations
  
  nodes <- dim(list.inc)[1]
  mat.sum <- matrix(numeric(nodes*nodes),nrow=nodes)
  
  L <- list()
  G <- list()
  
  # Compute the marginal edge posterior
  
  for (i in E:dim(list.inc)[3]){
    
    k <- cpdag(list.inc[,,i])    # cpdag of the i-th incidence matrix
    dummy <- matrix(numeric(nodes*nodes),nrow=nodes)
    
    if(length(nrow(k))!=0){
      dummy[k[,1]] <- k[,5]
      L[[i]] <- dummy
    }
    if(length(nrow(k))==0 && length(k)>0){
      dummy[k[1]] <- k[5]
      L[[i]] <- dummy
    }
    
    mat.com <-matrix(numeric(nodes*nodes),nrow=nodes)
    mat.re <- matrix(numeric(nodes*nodes),nrow=nodes)
    
    com <- which(L[[i]]>0)
    re <- which(L[[i]]<0)
    mat.com[com] <- 1
    mat.re[re] <- 1
    mat <- mat.com + mat.re + t(mat.re)
    G[[i]] <- mat
    mat.sum <- mat.sum + mat
  }
  
  return(list(L,G, (mat.sum/(length(list.inc)- E+1))))
}

# ---

extract_cpdag_from_dag <- function(true_incidence){    
  
  L <- list()
  
  nodes <- dim(true_incidence)[1]
  
  k <- cpdag(true_incidence)
  
  
  dummy <- matrix(numeric(nodes*nodes),nrow=nodes)
  
  if(length(nrow(k))!=0){
    dummy[k[,1]] <- k[,5]
    L <- dummy
  }
  
  if(length(nrow(k))==0 && length(k)>0){
    dummy[k[1]] <- k[5]
    L <- dummy
  }
  
  mat.com <- matrix(numeric(nodes*nodes),nrow=nodes)
  mat.re  <- matrix(numeric(nodes*nodes),nrow=nodes)
  
  com <- which(L>0)
  re  <- which(L<0)
  
  mat.com[com] <- 1
  mat.re[re] <- 1
  mat <- mat.com + mat.re + t(mat.re)
  
  return(mat)
}

PROPOSE_NEW_DAG <- function(incidence,ancestMat,fan.in){
  
  n = nrow(incidence)
  
  ### 1.) number of neighbour graphs by edge deletions
  num_deletion <- sum(incidence)
  
  ### 2.) number of neighbour graphs  by edge additions    1- E(i,j) - I(i,j) - A(i,j)
  inter_add <- which(matrix(rep(1,n*n),nrow=n) - diag(1,n,n) - incidence - ancestMat >0)
  add <- matrix(numeric(n*n),nrow=n)
  add[inter_add] <- 1
  add[,which(colSums(incidence)>fan.in-1)] <- 0
  num_addition <- sum(add)
  
  ### 3.) number of neighbour graphs by edge reversals    I - (I^t * A)^t
  inter_rev <- which(incidence - t(t(incidence)%*% ancestMat)==1)
  re <- matrix(numeric(n*n),nrow=n)
  re[inter_rev] <- 1
  re[which(colSums(incidence)>fan.in-1),] <- 0 
  num_reversal <- sum(re)
  
  ##### total number of neighbour graphs:
  total <- sum(num_deletion,num_addition,num_reversal)
  
  ### proposal probability:
  p_forward <- 1/total
  
  ############## sampling a new graph (or rather sampling an edge to shift)
  ### sample one of the three single edge operations
  random <- sample(1:total,1)
  
  operation <- 0 # if the single edge operation is (will be) an edge reversal
  if (random > total - num_reversal){operation <- 1}
  
  #### new incidence matrix
  incidence_new <- incidence
  
  if (random <= num_deletion){             # if edge deletion 
    if(length(which(incidence>0))>1){
      new_edge <- sample(which(incidence>0),1)} # sample one edge
    else{new_edge <- which(incidence>0)}
    incidence_new[new_edge] <- 0}            # and delete it
  
  if (random > (total - num_reversal)){      # if edge reversal was sampled
    if(num_reversal>1){
      new_edge <- sample(which(re==1),1)}  # sample one of the existing edges where a reversal leads to a valid graph
    else{
      new_edge <- which(re==1)}
    incidence_new[new_edge] <- 0             # delete it
    junk <- matrix(numeric(n*n),nrow=n)      # creating a matrix with all entries zero
    junk[new_edge] <- 1                      # an only a "1" at the entry of the new (reversed) edge
    incidence_new <- incidence_new + t(junk)}# sum the deleted matrix and the "junk-matrix"
  
  if (random <= (total - num_reversal) & random > num_deletion){     # if edge addition was sampled
    if(num_addition>1){
      new_edge <- sample(which(add==1),1)} # sample one of the existing edges where a addition leads to a valid graph
    else{
      new_edge <- which(add==1)}
    incidence_new[new_edge] <- 1             # and add it
  }
  
  
  #######################################################################################################################
  #################### Updating the ancestor matrix
  
  # creating a matrix with dimensions of the incidence matrix and all entries zero except for the entry of the chosen edge
  help_matrix <- matrix(numeric(n*n),nrow=n)
  help_matrix[new_edge] <- 1
  
  # numbers of the nodes that belong to the shifted edge
  parent <- which(rowSums(help_matrix)==1)
  child  <- which(colSums(help_matrix)==1)
  
  ### updating the ancestor matrix (after edge reversal)
  ## edge deletion
  ancestor_new <- ancestMat
  if (operation==1){
    ancestor_new[c(child,which(ancestMat[,child]==1)),] <- 0           # delete all ancestors of the child and its descendants                                           #
    top_name <- des_top_order(incidence_new, ancestMat, child)
    for (d in top_name){
      for(g in which(incidence_new[,d]==1)) {
        ancestor_new[d,c(g,(which(ancestor_new[g,]==1)))] <- 1
      }
    }
    ## edge addition
    anc_parent <- which(ancestor_new[child,]==1)                     # ancestors of the new parent
    des_child <- which(ancestor_new[,parent]==1)                     # descendants of the child
    ancestor_new[c(parent,des_child),c(child,anc_parent)] <- 1
  }
  
  ### updating the ancestor matrix (after edge deletion)
  if (random <= num_deletion){
    ancestor_new[c(child,which(ancestMat[,child]==1)),] <- 0           # delete all ancestors of the child and its descendants                                           #
    top_name <- des_top_order(incidence_new, ancestMat, child)
    for (d in top_name){
      for(g in which(incidence_new[,d]==1)) {
        ancestor_new[d,c(g,(which(ancestor_new[g,]==1)))] <- 1
      }
    }
  }
  
  # updating the ancestor matrix (after edge addition)
  if (random <= total - num_reversal & random > num_deletion){
    anc_parent <- which(ancestMat[parent,]==1)        # ancestors of the new parent
    des_child <- which(ancestMat[,child]==1)          # descendants of the child
    ancestor_new[c(child,des_child),c(parent,anc_parent)] <- 1
  }
  
  ######################################################################################################
  ####### ... the number of neighbour graphs/proposal probability for the proposed graph
  ### 1.) number of neighbour graphs obtained by edge deletions
  num_deletion_new <- sum(incidence_new)
  
  ### number of neighbour graphs obtained by edge additions    1- E(i,j) - I(i,j) - A(i,j)
  inter_add.new <- which(matrix(rep(1,n*n),nrow=n) - diag(1,n,n) - incidence_new - ancestor_new >0)
  add.new <- matrix(numeric(n*n),nrow=n)
  add.new[inter_add.new] <- 1
  add.new[,which(colSums(incidence_new)>fan.in-1)] <- 0
  num_addition_new <- sum(add.new)
  
  ### number of neighbour graphs obtained by edge reversals    I - (I^t * A)^t
  inter_rev.new <- which(incidence_new - t(t(incidence_new)%*% ancestor_new)==1)
  re.new <- matrix(numeric(n*n),nrow=n)
  re.new[inter_rev.new] <- 1
  re.new[which(colSums(incidence_new)>fan.in-1),] <- 0  # CORRECTED!!!???!!!
  num_reversal_new <- sum(re.new)
  
  ##### total number of neighbour graphs:
  total_new <- sum(num_deletion_new,num_addition_new,num_reversal_new)
  
  ### proposal probability:
  p_backward <- 1/total_new
  
  return(list(incidence_new,ancestor_new,p_forward,p_backward,operation,child,parent))
}

COMPUTE_BGE <- function(T_m,m,P_local_num,P_local_den,incidence,a,v,T_0){
  
  n = nrow(incidence)
  
  for (j in 1:n)  {
    n_nodes <- which(incidence[,j]==1)         # parents of j
    P_local_num[j] <- (-(length(n_nodes)+1)*m/2)*log(2*pi) + ((length(n_nodes)+1)/2)*log(v/(v+m)) + c_function((length(n_nodes)+1),a-n+length(n_nodes)+1)  - c_function((length(n_nodes)+1),a+m-n+length(n_nodes)+1)        + ((a-n+length(n_nodes)+1)/2)*log(det(as.matrix(T_0[sort(c(n_nodes,j)),sort(c(n_nodes,j))])))   + (-(a-n+length(n_nodes)+1+m)/2) *log(det(as.matrix(T_m[sort(c(n_nodes,j)),sort(c(n_nodes,j))])))
    if(sum(incidence[,j])>0){          # if j has at least one parent
      P_local_den[j] <- (-(length(n_nodes))*m/2)*log(2*pi) + (length(n_nodes)/2)*log(v/(v+m))     + c_function(length(n_nodes),    a-n+length(n_nodes)  )  - c_function(length(n_nodes),    a+m-n+length(n_nodes)  )        + ((a-n+length(n_nodes)  )/2)*log(det(as.matrix(T_0[n_nodes,n_nodes])))                         + (-(a-n+length(n_nodes)  +m)/2) *log(det(as.matrix(T_m[n_nodes,n_nodes])))
    }
    else{                              # if j has no parents
      P_local_den[j] <- 0
    }
  }
  
  
  bge_old = (sum(P_local_num))-(sum(P_local_den))
  
  return(list(bge_old,P_local_num,P_local_den))
  
}

COMPUTE_BGE_LOCAL <- function(operation,child,parent,T_m,m,P_local_num,P_local_den,incidence_new,a,v,T_0){
  
  n = nrow(incidence_new)
  
  P_local_num_new <- P_local_num
  P_local_den_new <- P_local_den
  
  n_nodes_new <- which(incidence_new[,child]==1)
  
  
  P_local_num_new[child] <- (-(length(n_nodes_new)+1)*m/2)*log(2*pi) + ((length(n_nodes_new)+1)/2)*log(v/(v+m)) + c_function((length(n_nodes_new)+1),a-n+length(n_nodes_new)+1) -c_function((length(n_nodes_new)+1),a+m-n+length(n_nodes_new)+1) + ((a-n+length(n_nodes_new)+1)/2)  *log(det(as.matrix(T_0[sort(c(n_nodes_new,child)),sort(c(n_nodes_new,child))])))+ (-(a-n+length(n_nodes_new)+1+m)/2)  *log(det(as.matrix(T_m[sort(c(n_nodes_new,child)),sort(c(n_nodes_new,child))])))
  
  if(sum(incidence_new[,child])>0){           # if at least one parent
    P_local_den_new[child] <- (-(length(n_nodes_new))*m/2)*log(2*pi) + (length(n_nodes_new)/2)*log(v/(v+m))     + c_function(length(n_nodes_new)    ,a-n+length(n_nodes_new)  ) - c_function(length(n_nodes_new),   a+m-n+length(n_nodes_new)  ) + ((a-n+length(n_nodes_new))  /2)  *log(det(as.matrix(T_0[n_nodes_new,n_nodes_new])))                              + (-(a-n+length(n_nodes_new)+m  )/2)  *log(det(as.matrix(T_m[n_nodes_new,n_nodes_new])))
  }
  else{                                       # if no parents
    P_local_den_new[child] <- 0
  }
  
  if (operation==1){                          # if single edge operation was an edge reversal
    
    n_nodesP <- which(incidence_new[,parent]==1)
    P_local_num_new[parent] <- (-(length(n_nodesP)+1)*m/2)*log(2*pi) + ((length(n_nodesP)+1)/2)*log(v/(v+m)) + c_function((length(n_nodesP)+1),a-n+length(n_nodesP)+1) - c_function((length(n_nodesP)+1),a+m-n+length(n_nodesP)+1)     + ((a-n+length(n_nodesP)+1)/2)*log(det(as.matrix(T_0[sort(c(n_nodesP,parent)),sort(c(n_nodesP,parent))])))+ (-(a+m-n+length(n_nodesP)+1)/2)  *log(det(as.matrix(T_m[sort(c(n_nodesP,parent)),sort(c(n_nodesP,parent))])))
    if(sum(incidence_new[,parent])>0){          # if parent at least one parent
      P_local_den_new[parent] <- (-(length(n_nodesP))*m/2)*log(2*pi) + (length(n_nodesP)/2)*log(v/(v+m))     + c_function(length(n_nodesP)    ,a-n+length(n_nodesP)  ) - c_function(length(n_nodesP)    ,a+m-n+length(n_nodesP)  )     + ((a-n+length(n_nodesP)) /2)*log(det(as.matrix(T_0[n_nodesP,n_nodesP])))                                + (-(a+m-n+length(n_nodesP))  /2)  *log(det(as.matrix(T_m[n_nodesP,n_nodesP])))
    }
    
    else{P_local_den_new[parent] <- 0}
  }
  
  
  bge_new <- (sum(P_local_num_new))-(sum(P_local_den_new))
  
  return(list(bge_new,P_local_num_new,P_local_den_new))
}
