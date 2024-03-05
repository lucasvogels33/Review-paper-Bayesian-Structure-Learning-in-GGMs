
############################
###SSO function#######
############################

#the function SSO takes as input:
#- data from BDgraph. This "data" type contains data, a true graph 
#- iter_vec_thin. The iteration at which it needs to produce output

#It outputs:
#- the AUC at every iteration number indicated in iter_vec_thin
#- the p+ at every iteration number indicated in iter_vec_thin
#- the p- at every iteration number indicated in iter_vec_thin
#- the running time at every iteration number indicated in iter_vec_thin

SSO_solve = function(data,iter_vec_thin){
  
  #obtain true graph
  G_true = data$G
  response = G_true[upper.tri(G_true)]
  
  #obtain amount of observations and variables
  n = nrow(data$data)
  p = ncol(data$data)
  density = data$density
  graph = data$graph
  rep = data$rep
  algorithm_name = "SSO"
  
  #parameters
  burnin = 0
  var1 = 0.02
  var2 = 2
  lambda = 1
  save = FALSE
  cores = 1
  MCMCstart = "empty"
  sig.start = NULL
  if (p > n){sig.start = diag(1,p)}
  g.prior = 0.2
  plinks_previous = 0
  verbose = FALSE
  
  #run the algorithm, produce output at every iteration in the iter_vec_thin vector
  len = length(iter_vec_thin)
  AUC_vec = rep(0,len)
  pplus_vec = rep(0,len)
  pmin_vec = rep(0,len)
  time_vec = rep(0,len)
  
  for (j in 1:len){
    
    #determine amount of iterations
    olditer = 0
    if (j > 1){
      olditer = iter_vec_thin[j-1]
    }
    newiter = iter_vec_thin[j]
    iter = newiter-olditer
    
    #run bdgraph for maxiter iterations
    sample_ss  = ssgraph( data = data, iter = iter, burnin = burnin, var1 = var1, var2=var2,lambda=lambda,g.start=MCMCstart,sig.start=sig.start,save=save,cores=cores,g.prior=g.prior,verbose=verbose)
    
    #calculate new p matrix
    plinks_new = (olditer*plinks_previous + iter*sample_ss$p_links)/newiter
    predictor = plinks_new[upper.tri(plinks_new)]
    
    #calculate the AUC and the pplus and pmin
    AUC = calc_AUC_ROC(predictor=predictor,response=response)
    obj = calc_pplus_pmin(predictor=predictor,response=response)
    p_plus = obj$p_plus
    p_min = obj$p_min
    
    #save metrics
    AUC_vec[j] = AUC
    pplus_vec[j] = p_plus
    pmin_vec[j] = p_min
    
    #calculate the runtime
    time_init = 0
    if (j==1){time_init = sample_ss$time_init}
    runtime = time_init + sample_ss$time_all_iterations
    time_vec[j] = runtime
    
    #store data for next run
    plinks_previous = plinks_new
    MCMCstart = sample_ss
    
    #save output per iteration in Rdata file (only for large scale problems)
    if (p > 250){
      output_per_iter = list()
      output_per_iter$AUC = AUC
      output_per_iter$pplus = p_plus
      output_per_iter$pmin = p_min
      output_per_iter$time = sum(time_vec)/10^6
      output_per_iter$iter = newiter
      filename = paste0("results_p",p,"_n",n,"_",graph,"_dens",density,"_",algorithm_name,"_rep",rep,"_iter",j,".Rdata")
      save(output_per_iter,file = filename)
    }
  }
  time_vec = cumsum(time_vec)
  time_vec = time_vec/10^6 #ssgraph outputs the time in microseconds. Here we convert it to seconds
  
  return(list(AUC_vec=AUC_vec,pplus_vec=pplus_vec,pmin_vec=pmin_vec,time_vec=time_vec,iter_vec_thin=iter_vec_thin,predictor=predictor))
} 

