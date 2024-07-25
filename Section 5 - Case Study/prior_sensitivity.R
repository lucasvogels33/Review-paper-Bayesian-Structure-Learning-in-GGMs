library(BDgraph)
library(ssgraph)
library(huge)

#set parameters
algorithm = "BDA" #set algorithm to either BDA or SS
g.prior = 0.1 #set prior to any number between 0 and 1

##########################################
######load and normalize data ############
##########################################

data(geneExpression) #load data
gene_data = geneExpression 
new_data = huge.npn(x=geneExpression,npn.func="shrinkage") #new data is n x p matrix containing normalized data

##########################################
########define functions ###### ############
##########################################

#functions that takes as input:
#- the data 
#- the burnin iterations at which it needs to output the #edges in the Markov state, 
#- the iterations at which it needs to output the #edges in the Markov state, 
#- the initial graph of the Markov Chain

#as output:
#- edge_vec 
#- time_vec
#- iter_vec_thin
#- burnin_iter_vec_thin
#- plinks
#- plinks_iter_matrix 


SS_solve= function(data,burnin_iter_vec_thin,iter_vec_thin,MCMCstart,g.prior){
  
  #obtain p and n
  p = ncol(data)
  n = nrow(data)
  
  #parameters
  burnin = 0
  var1 = 0.02
  var2 = 2
  lambda = 1
  save = FALSE
  cores = 1
  sig.start = NULL
  if (p > n){sig.start = diag(1,p)}
  verbose = FALSE
  
  #intialize outputs
  len_iter = length(iter_vec_thin)
  len_burnin = length(burnin_iter_vec_thin)
  edge_vec = rep(0,len_iter+len_burnin)
  time_vec = rep(0,len_iter+len_burnin)
  
  #run burnin_iterations
  for (j in 1:len_burnin){
    print(burnin_iter_vec_thin[j]) #output progress
    
    #determine amount of iterations
    olditer = 0
    if (j > 1){
      olditer = burnin_iter_vec_thin[j-1]
    }
    newiter = burnin_iter_vec_thin[j]
    iter = newiter-olditer
    
    #run bdgraph for maxiter iterations
    sample_ss  = ssgraph( data = data,iter = iter, burnin = burnin, var1 = var1, var2=var2,lambda=lambda,g.start=MCMCstart,sig.start=sig.start,save=save,cores=cores,g.prior=g.prior,verbose=verbose)
    
    #calculate the runtime
    time_init = 0
    if (j==1){time_init = sample_ss$time_init}
    runtime = time_init + sample_ss$time_all_iterations
    
    #save metrics
    edge_vec[j] = sum(sample_ss$last_graph)/2
    time_vec[j] = runtime
    
    #save data for next run
    MCMCstart = sample_ss
    
  }
  
  #run MCMC iterations after burnin
  plinks_old = 0
  plinks_iter_matrix = matrix(0,p*(p-1)/2,len_iter)
  for (j in 1:len_iter){
    print(iter_vec_thin[j]) #output progress
    
    #determine amount of iterations
    olditer = 0
    if (j > 1){
      olditer = iter_vec_thin[j-1]
    }
    newiter = iter_vec_thin[j]
    iter = newiter-olditer
    
    #run bdgraph for maxiter iterations
    sample_ss  = ssgraph( data = data,iter = iter, burnin = burnin, var1 = var1, var2=var2,lambda=lambda,g.start=MCMCstart,sig.start=sig.start,save=save,cores=cores,g.prior=g.prior,verbose=verbose)
    
    #calculate new p matrix
    plinks_new = (olditer*plinks_old + iter*sample_ss$p_links)/(newiter)
    
    #calculate the runtime
    time_init = 0
    if (j==1){time_init = sample_ss$time_init}
    runtime = time_init + sample_ss$time_all_iterations
    
    #save metrics
    edge_vec[len_burnin + j] = sum(sample_ss$last_graph)/2
    time_vec[len_burnin + j] = runtime
    
    #save data for next run
    plinks_old = plinks_new
    MCMCstart = sample_ss

    #save plinks in matrix
    plinks_vec = plinks_new[upper.tri(plinks_new)]
    plinks_iter_matrix[,j] = plinks_vec
  }
  
  
  time_vec = cumsum(time_vec)
  time_vec = time_vec/(10^6)
  return(list(g.prior=g.prior,edge_vec=edge_vec,time_vec=time_vec,iter_vec_thin=iter_vec_thin,burnin_iter_vec_thin=burnin_iter_vec_thin,plinks=plinks_new,plinks_iter_matrix=plinks_iter_matrix))
}
BDA_solve = function(data,burnin_iter_vec_thin,iter_vec_thin,MCMCstart,g.prior){
  
  #obtain p and n
  p = ncol(data)
  n = nrow(data)

  #parameters
  burnin = 0
  jump = 1
  save = FALSE
  cores = 1
  verbose = FALSE
  
  #run the algorithm, produce output at every iteration in the iter_vec_thin vector
  len_iter = length(iter_vec_thin)
  len_burnin = length(burnin_iter_vec_thin)
  edge_vec = rep(0,len_iter+len_burnin)
  time_vec = rep(0,len_iter+len_burnin)
  
  for (j in 1:len_burnin){
    
    #determine amount of iterations
    olditer = 0
    if (j > 1){
      olditer = burnin_iter_vec_thin[j-1]
    }
    newiter = burnin_iter_vec_thin[j]
    iter = newiter-olditer
    
    #run bdgraph for maxiter iterations
    sample_bd  = bdgraph( data = data,algorithm = "bdmcmc", iter = iter, burnin = burnin, jump = jump, save = save,cores=cores,g.start=MCMCstart,g.prior=g.prior,verbose=verbose)  
    
    #calculate the runtime
    time_init = 0
    time_end = 0
    if (j==1){time_init = sample_bd$time_init}
    if (j==len_burnin){time_end = sample_bd$time_end}
    runtime = time_init + sample_bd$time_all_iterations + time_end
    
    #save metrics
    edge_vec[j] = sum(sample_bd$last_graph)/2
    time_vec[j] = runtime
    
    #save data for next run
    MCMCstart = sample_bd
    
  }
  
  #run MCMC iterations (after burnin)
  weights_old = 0
  plinks_old = 0
  plinks_iter_matrix = matrix(0,p*(p-1)/2,len_iter)
  for (j in 1:len_iter){
    
    #determine amount of iterations
    olditer = 0
    if (j > 1){
      olditer = iter_vec_thin[j-1]
    }
    newiter = iter_vec_thin[j]
    iter = newiter-olditer
    
    #run bdgraph for maxiter iterations
    sample_bd  = bdgraph( data = data,algorithm = "bdmcmc", iter = iter, burnin = burnin, jump = jump, save = save,cores=cores,g.start=MCMCstart,g.prior=g.prior,verbose=verbose)  
    
    #calculate new p matrix
    weights_new = sample_bd$sum_weights
    plinks_new = (weights_old*plinks_old + weights_new*sample_bd$p_links)/(weights_old+weights_new)
    
    #calculate the runtime
    time_init = 0
    time_end = 0
    if (j==1){time_init = sample_bd$time_init}
    if (j==len_iter){time_end = sample_bd$time_end}
    runtime = time_init + sample_bd$time_all_iterations + time_end
    
    #save metrics
    edge_vec[len_burnin+j] = sum(sample_bd$last_graph)/2
    time_vec[len_burnin+j] = runtime
    
    #save data for next run
    plinks_old = plinks_new
    weights_old = weights_old+weights_new
    MCMCstart = sample_bd

    #save plinks in matrix
    plinks_vec = plinks_new[upper.tri(plinks_new)]
    plinks_iter_matrix[,j] = plinks_vec
    
  }
  
  plinks = plinks_new
  time_vec = cumsum(time_vec)
  time_vec = time_vec/(10^6)
  return(list(g.prior=g.prior,edge_vec=edge_vec,time_vec=time_vec,iter_vec_thin=iter_vec_thin,burnin_iter_vec_thin=burnin_iter_vec_thin,plinks=plinks_new,plinks_iter_matrix=plinks_iter_matrix))
} 

##########################################
########run algorithms ###### ############
##########################################

data = new_data
MCMCstart = "empty"
if (algorithm == "SS"){
    burnin_iter_vec_thin = c(100)
    iter_vec_thin = c(5000)
    output = SS_solve(data=data,burnin_iter_vec_thin=burnin_iter_vec_thin,iter_vec_thin=iter_vec_thin,MCMCstart=MCMCstart,g.prior=g.prior)
}
if (algorithm == "BDA"){
    burnin_iter_vec_thin = c(25000)
    iter_vec_thin = c(250000,500000,750000)
    output = BDA_solve(data=data,burnin_iter_vec_thin=burnin_iter_vec_thin,iter_vec_thin=iter_vec_thin,MCMCstart=MCMCstart,g.prior=g.prior)
}
filename = paste0("prior_",g.prior,"_",algorithm,".Rdata")
save(output,file=filename)

