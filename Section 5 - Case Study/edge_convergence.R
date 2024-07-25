#set parameters
algorithm = "BDA" #set algorithm to either BDA or SS
density = 0.1 #set density to any number between 0 and 1

#download libraries
library(BDgraph)
library(ssgraph)
library(huge)

#######################################
#####load data ########################
######################################

data(geneExpression) #load data
gene_data = geneExpression 
new_data = huge.npn(x=geneExpression,npn.func="shrinkage") #new data is n x p matrix containing normalized data

#############################################################
################supporting functions #######################
#############################################################

#functions that takes as input:
#- the data 
#- the burnin iterations at which it needs to output the #edges in the Markov state, 
#- the iterations at which it needs to output the #edges in the Markov state, 
#- the initial graph of the Markov Chain

#It solves the data and outputs the #edges per Markov state
SS_solve= function(data,burnin_iter_vec_thin,iter_vec_thin,MCMCstart){
  
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
  g.prior = 0.2
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
    
  }
  
  
  time_vec = cumsum(time_vec)
  time_vec = time_vec/(10^6)
  return(list(edge_vec=edge_vec,time_vec=time_vec,iter_vec_thin=iter_vec_thin,burnin_iter_vec_thin=burnin_iter_vec_thin,plinks=plinks_new))
}
BDA_solve = function(data,burnin_iter_vec_thin,iter_vec_thin,MCMCstart){
  
  #parameters
  burnin = 0
  jump = 1
  save = FALSE
  cores = 1
  g.prior = 0.2
  verbose = TRUE
  
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
    
  }
  
  plinks = plinks_new
  time_vec = cumsum(time_vec)
  time_vec = time_vec/(10^6)
  return(list(edge_vec=edge_vec,time_vec=time_vec,iter_vec_thin=iter_vec_thin,burnin_iter_vec_thin=burnin_iter_vec_thin,plinks=plinks_new))
} 


#function that takes as input:
#- the algorithm
#- the data 
#- the densities of the Markov Chain initialization
#and outputs a matrix with #edges per iteration (column) per density (row)
edge_per_state =function(algorithm,data,burnin_iter_vec_thin,iter_vec_thin,densities){
  
  #calculate dimension and amount of links
  n = dim(data)[1]
  p = dim(data)[2]
  qp = p*(p-1)/2
  
  #create empty output matrix
  edge_matrix = matrix(0,length(densities),length(iter_vec_thin)+length(burnin_iter_vec_thin)+1)
  
  #fill output matrix
  q=0
  for (density in densities){
    q=q+1
    
    #create initial graph
    init_graph = array(0,c(p,p))
    upper_index = upper.tri(init_graph,diag=FALSE)
    init_graph[upper_index] = rbinom(qp,1,density)
    for (i in 2:p){
      for (j in 1:i){
        init_graph[i,j] = init_graph[j,i]
      }
    }
    
    #run algorithms 
    if (algorithm == "SS"){out = SS_solve(data=data,burnin_iter_vec_thin=burnin_iter_vec_thin,iter_vec_thin=iter_vec_thin,MCMCstart=init_graph)}
    if (algorithm == "BDA"){out = BDA_solve(data=data,burnin_iter_vec_thin=burnin_iter_vec_thin,iter_vec_thin=iter_vec_thin,MCMCstart=init_graph)}
    
    #save result
    edge_vec = out$edge_vec
    edge_vec = c(sum(init_graph)/2,edge_vec)
    edge_matrix[q,]=edge_vec
    
    burnin_total = burnin_iter_vec_thin[length(burnin_iter_vec_thin)]
    all_iter_vec = c(burnin_iter_vec_thin,burnin_total+iter_vec_thin)
    all_iter_vec = c(0,all_iter_vec)

  }
  
  return(list(edge_matrix=edge_matrix,all_iter_vec=all_iter_vec))
}

######################################
#######solve data and save results####
######################################

densities = c(density)

if (algorithm == "SS"){
  iter_vec_thin = c(20,40,60,80,100,125,150,175,200,250,300,400,500,750,1000)
  burnin_iter_vec_thin = c(1,2,3,4,5,10,20)
}
if (algorithm == "BDA"){
  iter_vec_thin = c(5000,10000,15000,20000,25000,30000)
  burnin_iter_vec_thin = c(10,20,30,40,50,100,200,300,400,500,750,1000,1500,2000,2500,3000,4000,5000,7500,10000,15000)
}

data = new_data
output = edge_per_state(algorithm=algorithm,data=data,burnin_iter_vec_thin=burnin_iter_vec_thin,iter_vec_thin=iter_vec_thin,densities=densities)

filename = paste0("edge_convergence_",algorithm,"_density",density,".Rdata")
save(output,file = filename)





