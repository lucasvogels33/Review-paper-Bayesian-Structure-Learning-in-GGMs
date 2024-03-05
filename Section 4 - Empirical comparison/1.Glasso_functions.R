

###################
####functions######
###################
source("metric_functions.R")

#########################
###glasso function#######
########################

#the function glasso takes as input data from BDgraph. This "data" type contains data and a true graph\
#It outputs:
#- the AUC based on a vector of lambda
#- the time to compute a path of graphs for each lambda
#- the time to select a labmda for ric, ebic and stars
#- final p+ and p- values for ric, ebic and stars


glasso_solve = function(data){

  #set parameters
  x = data$data
  nlambda = 30
  lambda.min.ratio = 0.01
  method = "glasso"
  verbose = TRUE
  
  ##############
  ####tuning####
  ##############
  
  #solve data for a vector of lambda
  start      = proc.time()
  huge_output = huge(
    x=x,
    lambda = NULL,
    nlambda = nlambda,
    lambda.min.ratio = lambda.min.ratio,
    method = method,
    verbose = verbose
  )
  time_tuning_prep = as.numeric( ( proc.time() - start)[ 3 ] ) #end and save tim
  
  time_tuning_ric = NULL
  time_tuning_stars = NULL
  time_tuning_ebic = NULL

  #select one lambda based on either ric, stars or ebic
  start      = proc.time()
  ric_out = huge.select(huge_output,criterion="ric") 
  time_tuning_ric = as.numeric( ( proc.time() - start)[ 3 ] ) #end and save tim
  
  if (p < 250){#the stars tuning method is too time consuming for large scale problems
    start      = proc.time()
    stars_out =huge.select(huge_output,criterion="stars",stars.thresh=0.05,rep.num=10) 
    time_tuning_stars = as.numeric( ( proc.time() - start)[ 3 ] ) #end and save time
  }

  start      = proc.time()
  ebic_out = huge.select(huge_output,criterion="ebic")
  time_tuning_ebic = as.numeric( ( proc.time() - start)[ 3 ] ) #end and save tim
  
  #compute AUC
  G_true = data$G
  response = G_true[upper.tri(G_true)]
  
  lambda_vec = huge_output$lambda
  lambda_max = max(lambda_vec)
  lambda_vec = lambda_vec/lambda_max
  predictor = rep(0,p*(p-1)/2)
  len = length(lambda_vec)
  G_path = huge_output$path
  for (i in 2:len) {
    G_new = matrix(unlist(G_path[i]),ncol = p)
    G_previous = matrix(unlist(G_path[i-1]),ncol = p)
    G_diff = G_new - G_previous
    G_diff = G_diff[upper.tri(G_diff)]
    predictor[which(G_diff==1)] = lambda_vec[i]
  }
  
  AUC = calc_AUC_ROC(response=response,predictor=predictor)
  
  ###########################
  #####model selection ######
  ###########################
  
  time_running_ric = NULL
  time_running_stars = NULL
  time_running_ebic = NULL
           
  pplus_ric = NULL
  pplus_stars = NULL
  pplus_ebic = NULL

  predictor_ric = NULL
  predictor_stars = NULL
  predictor_ebic = NULL
           
  pmin_ric = NULL
  pmin_stars = NULL
  pmin_ebic = NULL

  #ric
  start      = proc.time()
  huge_output_ric = huge(x=x,lambda = ric_out$opt.lambda,method=method)
  time_running_ric = as.numeric( ( proc.time() - start)[ 3 ] )
  G_est = huge_output_ric$path[[1]]
  predictor_ric = G_est[upper.tri(G_est)]
  TPR_FPR_obj = calc_TPR_FPR(predictor=predictor_ric,response=response)
  pplus_ric = TPR_FPR_obj$TPR
  pmin_ric = TPR_FPR_obj$FPR
  
  if (p < 250){#the stars tuning method is too time consuming for large scale problems
    #stars
    start      = proc.time()
    huge_output_stars = huge(x=x,lambda = stars_out$opt.lambda,method=method)
    time_running_stars = as.numeric( ( proc.time() - start)[ 3 ] )
    G_est = huge_output_stars$path[[1]]
    predictor_stars = G_est[upper.tri(G_est)]
    TPR_FPR_obj = calc_TPR_FPR(predictor=predictor_stars,response=response)
    pplus_stars = TPR_FPR_obj$TPR
    pmin_stars = TPR_FPR_obj$FPR
  }

  #ebic
  start      = proc.time()
  huge_output_ebic = huge(x=x,lambda = ebic_out$opt.lambda,method=method)
  time_running_ebic = as.numeric( ( proc.time() - start)[ 3 ] )
  G_est = huge_output_ebic$path[[1]]
  predictor_ebic = G_est[upper.tri(G_est)]
  TPR_FPR_obj = calc_TPR_FPR(predictor=predictor_ebic,response=response)
  pplus_ebic = TPR_FPR_obj$TPR
  pmin_ebic = TPR_FPR_obj$FPR
  
  
  
  #create list to output
  return(list(
           AUC = AUC,
           
           time_tuning_prep=time_tuning_prep,
           time_tuning_ric = time_tuning_ric,
           time_tuning_stars = time_tuning_stars,
           time_tuning_ebic = time_tuning_ebic,
           
           time_running_ric = time_running_ric,
           time_running_stars = time_running_stars,
           time_running_ebic = time_running_ebic, 
           
           pplus_ric = pplus_ric,
           pplus_stars = pplus_stars,
           pplus_ebic = pplus_ebic,

           predictor_ric = predictor_ric,
           predictor_stars = predictor_stars,
           predictor_ebic = predictor_ebic,
           
           pmin_ric = pmin_ric,
           pmin_stars = pmin_stars,
           pmin_ebic = pmin_ebic
         )
  )

}


