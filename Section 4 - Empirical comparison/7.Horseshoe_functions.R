
############################
###supporting functions####
############################

GHS = function(S,n,maxiter,Omega_start,Sigma_start,Lambda_sq_start,Nu_start,tau_sq_start,xi_start){
  p = nrow(S)
  omega_save = array(0,dim=c(p*(p-1)/2,maxiter))
  lambda_sq_save = array(0,dim=c(p*(p-1)/2,maxiter))
  tau_sq_save = rep(0,maxiter)
  
  
  ind_all = array(0,dim=c(p-1,p))
  
  for (i in 1:p){
    if (i == 1) {
      ind = 2:p
    } 
    if (i ==p) {
      ind = 1:(p-1)
    } 
    if ((i > 1)*(i<p)){
      ind = c(1:(i-1),(i+1):p)
    }
    
    ind_all[,i] = ind
  }
  
  #set initial values
  Omega = Omega_start
  Sigma = Sigma_start
  Lambda_sq = Lambda_sq_start
  Nu = Nu_start
  tau_sq = tau_sq_start
  xi = xi_start
  
  for (iter in 1: maxiter){
    for (i in 1:p){
      ind = ind_all[,i]     
      Sigma_11 = Sigma[ind,ind]
      sigma_12 = Sigma[ind,i]
      sigma_22 = Sigma[i,i]
      s_21 = S[ind,i]
      s_22 = S[i,i]
      lambda_sq_12 = Lambda_sq[ind,i]
      nu_12 = Nu[ind,i]
      
      #sample gamma and beta
      gamma = rgamma(1,n/2 + 1,s_22/2) #check if it samples the same value as in matlab
      inv_Omega_11 = Sigma_11 - sigma_12%*%t(sigma_12)/sigma_22
      inv_C = s_22*inv_Omega_11+diag(1/(lambda_sq_12*tau_sq),p-1)
      inv_C_chol = chol(inv_C)
      mu_i = solve(-inv_C, s_21,tol=1e-25)
      qq = rnorm(p-1,0,1)
      beta = mu_i+ solve(inv_C_chol,qq,tol=1e-25)
      omega_12 = beta
      omega_22 = gamma + beta %*% inv_Omega_11 %*% beta
      
      #sample lambda_sq and nu
      rate = (omega_12^2)/(2*tau_sq)+1/nu_12
      lambda_sq_12 = 1/rgamma(length(rate),1,rate)    
      nu_12 = 1/rgamma(length(lambda_sq_12),1,1+1/lambda_sq_12)
      
      #pdate Omega, Sigma, Lambda_sq, Nu
      Omega[i,ind] = omega_12
      Omega[ind,i] = omega_12
      Omega[i,i] = omega_22
      temp = inv_Omega_11%*%beta
      Sigma_11 = inv_Omega_11 + temp %*% t(temp) /gamma
      
      sigma_12 = -temp/gamma
      sigma_22 = 1/gamma
      
      Sigma[ind,ind] = Sigma_11
      Sigma[i,i] = sigma_22
      Sigma[i,ind] = sigma_12
      Sigma[ind,i] = sigma_12
      Lambda_sq[i,ind] = lambda_sq_12
      Lambda_sq[ind,i] = lambda_sq_12
      Nu[i,ind] = nu_12
      Nu[ind,i] = nu_12
    }   
    
    # sample tau_sq and xi
    omega_vector = Omega[upper.tri(Omega)]
    lambda_sq_vector = Lambda_sq[upper.tri(Lambda_sq)]
    rate = 1/xi + sum(  omega_vector^2 / (2*lambda_sq_vector))
    tau_sq = 1/rgamma(1,(p*(p-1)/2+1)/2, rate)
    xi = 1/rgamma(1,1,(1+1/tau_sq));    
    
    #save Omega, lambda_sq, tau_sq
    omega_save[,iter] = omega_vector
    lambda_sq_save[,iter] = lambda_sq_vector
    tau_sq_save[iter] = tau_sq
    
  }
  
  return(list(omega_save=omega_save,Omega=Omega,Sigma=Sigma,Lambda_sq=Lambda_sq,Nu=Nu,tau_sq=tau_sq,xi=xi))
  
}

############################
###horseshoe function#######
############################

#the function horseshoe_solve takes as input:
#- data from BDgraph. This "data" type contains data, a true graph 
#- iter_vec_thin. The iteration at which it needs to produce output

#It outputs:
#- the AUC at every iteration number indicated in iter_vec_thin
#- the p+ at every iteration number indicated in iter_vec_thin
#- the p- at every iteration number indicated in iter_vec_thin
#- the running time at every iteration number indicated in iter_vec_thin

horseshoe_solve = function(data,iter_vec_thin){
  
  
  #obtain true graph
  G_true = data$G
  response = G_true[upper.tri(G_true)]
  
  #obtain amount of observations and variables
  n = nrow(data$data)
  p = ncol(data$data)
  density = data$density
  graph = data$graph
  rep = data$rep
  algorithm_name = "Horseshoe"
  qp = p*(p-1)/2
  
  #obtain S = Y^t*Y
  S<-t(data.sim$data) %*% data.sim$data #this is how they calculate S in their original matlab algorithm. However, the scatter matrix S should be calculated using the sample means too.
  
  #initialize parameters
  Omega_start = diag(1,p)
  Sigma_start = diag(1,p)
  Lambda_sq_start = array(1,dim=c(p,p))
  Nu_start = array(1,dim=c(p,p))
  tau_sq_start = 1
  xi_start = 1
  
  #set the omega_save vector to be empty.
  omega_save = c()
  
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
    
    #run horseshoe for maxiter iterations
    start = proc.time()
    GHS_output = GHS(S=S,n=n,maxiter=maxiter,
                     Omega_start=Omega_start,
                     Sigma_start=Sigma_start,
                     Lambda_sq_start=Lambda_sq_start,
                     Nu_start=Nu_start,
                     tau_sq_start=tau_sq_start,
                     xi_start=xi_start)
    runtime = as.numeric( ( proc.time() - start)[ 3 ] ) #end and save time
    
    #save output for next run
    Omega_start = GHS_output$Omega 
    Sigma_start = GHS_output$Sigma
    Lambda_sq_start = GHS_output$Lambda_sq
    Nu_start = GHS_output$Nu
    tau_sq_start = GHS_output$tau_sq
    xi_start = GHS_output$xi
    
    #add the new omega path to the existing omega_path
    omega_save = matrix(c(omega_save,GHS_output$omega_save),nrow=p*(p-1)/2)
    
    #find the estimated graph based on the 50th percentile
    G_est = rep(0,qp)
    for (i in 1:qp){
      lb = quantile(omega_save[i,],0.25)
      ub = quantile(omega_save[i,],0.75)
      if ((lb >0) | (ub <0)){
        G_est[i] = 1
      }
    }
    
    #calculate the pplus and pmin of the estimated graph
    obj_pmin_plus = calc_pplus_pmin(predictor=G_est,response=response)
    pplus = obj_pmin_plus$p_plus
    pmin = obj_pmin_plus$p_min
    
    #calculate AUC
    nlambda = 30
    lambda_vec = seq(0,1,1/nlambda)
    G_old = rep(0,qp)
    G_new = rep(0,qp)
    predictor = rep(0,qp)
    
    for (lambda in lambda_vec) {
      lower_prct = lambda/2
      upper_prct = 1 - lower_prct
      for (i in 1:qp){
        lb = quantile(omega_save[i,],lower_prct)
        ub = quantile(omega_save[i,],upper_prct)
        if ((lb >0) | (ub <0)){
          G_new[i] = 1
        }
        
      }
      G_diff = G_new - G_old
      G_old = G_new
      predictor[which(G_diff==1)] = 1-lambda
    }
    AUC = calc_AUC_ROC(response=response,predictor=predictor)
    
    #save metrics
    AUC_vec[j] = AUC
    pplus_vec[j] = pplus
    pmin_vec[j] = pmin
    time_vec[j] = runtime
    
    if (p > 250){
      output_per_iter = list()
      output_per_iter$AUC = AUC
      output_per_iter$pplus = pplus
      output_per_iter$pmin = pmin
      output_per_iter$time = sum(time_vec)
      output_per_iter$iter = iter_vec_thin[j]
      filename = paste0("results_p",p,"_n",n,"_",graph,"_dens",density,"_",algorithm_name,"_rep",rep,"_iter",j,".Rdata")
      save(output_per_iter,file = filename)
    }
  }
  time_vec = cumsum(time_vec)
  
  return(list(AUC_vec=AUC_vec,pplus_vec=pplus_vec,pmin_vec=pmin_vec,time_vec=time_vec,iter_vec_thin=iter_vec_thin,predictor=G_est))
} 
  


