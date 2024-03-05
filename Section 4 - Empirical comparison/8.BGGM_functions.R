source("metric_functions.R")

#########################
###BGGM function#######
########################

#the function BGGM takes as input data from BDgraph. This "data" type contains data and a true graph\
#It outputs:
#- the AUC 
#- p+ and p- 
#- the time

BGGM_solve = function(data.sim){
  #save created data
  data = data.sim$data
  G_true = data.sim$G
  response= G_true[upper.tri(G_true)]
  
  ####solve data
  start = proc.time()
  fit = explore(data)
  E = select(fit, alternative = "exhaustive")
  runtime = as.numeric( ( proc.time() - start)[ 3 ] ) #end and save time
  
  #calculate metrics
  predictor = 1 - E$post_prob$prob_zero
  AUC = calc_AUC_ROC(predictor=predictor,response=response)
  obj = calc_pplus_pmin(predictor=predictor,response=response)
  p_plus = obj$p_plus
  p_min = obj$p_min
  
  return(list(AUC=AUC,pplus=p_plus,pmin=p_min,time=runtime,predictor=predictor))
}

