

source("metric_functions.R")

############################
###supporting  functions#######
############################

#this function produces an update of the graph using the WWA algorithm
update_G <- function(adj, edge_prob, df_0, U, n) {
  p <- nrow(adj)
  
  return(update_G_Rcpp(
    adj = adj,
    edge_prob_mat = matrix(edge_prob, nrow = p, ncol = p),
    df = df_0 + n,
    df_0 = df_0,
    rate = diag(p) + U,
    n_edge = p,  # Number of single edge updates attempted in this MCMC update
    seed = sample.int(n = .Machine$integer.max, size = 1),
    loc_bal = FALSE
  ))
}

# this function samples from the G-Wishart distribution
rgwish <- function(adj, df = 3, rate = NULL) {
  p <- nrow(adj)
  if (is.null(rate)) rate <- diag(p)
  
  return(rgwish_Rcpp(
    adj = adj,
    df = df,
    rate = rate,
    seed = sample.int(n = .Machine$integer.max, size = 1)
  ))
}


############################
###RJ-WWA function#######
############################

#the function RJWWA_solve takes as input:
#- data from BDgraph. This "data" type contains data, a true graph 
#- iter_vec_thin. The iteration at which it needs to produce output

#It outputs:
#- the AUC at every iteration number indicated in iter_vec_thin
#- the p+ at every iteration number indicated in iter_vec_thin
#- the p- at every iteration number indicated in iter_vec_thin
#- the running time at every iteration number indicated in iter_vec_thin

RJWWA_solve = function(data,iter_vec_thin){
  
  
  #obtain true graph
  G_true = data$G
  response = G_true[upper.tri(G_true)]
  
  #obtain amount of observations and variables
  n = nrow(data$data)
  p = ncol(data$data)
  density = data$density
  graph = data$graph
  rep = data$rep
  algorithm_name = "RJWWA"
  
  #initialize adjaceny matrices
  adj <- matrix(0L, nrow = p, ncol = p) #initial graph of the Markov chain
  adj_cum <- matrix(0L, nrow = p, ncol = p)
  
  #intialize paramters
  edge_prob = 0.2 #prior edge inclusion prob
  df_0 = 3.0 #degree for G-wishart distribution
  U <- t(data$data) %*% data$data
  
  #run the algorithm, produce output at every iteration in the iter_vec_thin vector
  len = length(iter_vec_thin)
  AUC_vec = rep(0,len)
  pplus_vec = rep(0,len)
  pmin_vec = rep(0,len)
  time_vec = rep(0,len)
  
  for (j in 1:len){
    
    #determine amount of iterations
    if (j==1){
      maxiter = iter_vec_thin[1]
    }
    if (j>1){
      maxiter = iter_vec_thin[j] - iter_vec_thin[j-1] 
    }
    
    #run RJWWA for maxiter iterations
    start = proc.time()
    for (s in 1:maxiter) {
      adj = update_G(adj = adj, edge_prob = edge_prob, df_0 = df_0, U = U, n = n)
      adj_cum = adj_cum + adj
    }
    runtime = as.numeric( ( proc.time() - start)[ 3 ] ) #end and save time
    
    #save edge inclusion probabilities 
    plinks  = adj_cum / iter_vec_thin[j]
    predictor = plinks[upper.tri(plinks)]
    
    #calculate the AUC and the pplus and pmin
    AUC = calc_AUC_ROC(predictor=predictor,response=response)
    obj = calc_pplus_pmin(predictor=predictor,response=response)
    p_plus = obj$p_plus
    p_min = obj$p_min
    
    #save metrics
    AUC_vec[j] = AUC
    pplus_vec[j] = p_plus
    pmin_vec[j] = p_min
    time_vec[j] = runtime
    
    #save output per iteration in Rdata file (only for large scale problems)
    if (p > 250){
      output_per_iter = list()
      output_per_iter$AUC = AUC
      output_per_iter$pplus = p_plus
      output_per_iter$pmin = p_min
      output_per_iter$time = sum(time_vec)
      output_per_iter$iter = iter_vec_thin[j]
      filename = paste0("results_p",p,"_n",n,"_",graph,"_dens",density,"_",algorithm_name,"_rep",rep,"_iter",j,".Rdata")
      save(output_per_iter,file = filename)
    }
  }
  time_vec = cumsum(time_vec)
  
  return(list(AUC_vec=AUC_vec,pplus_vec=pplus_vec,pmin_vec=pmin_vec,time_vec=time_vec,iter_vec_thin=iter_vec_thin,predictor=predictor))
} 



