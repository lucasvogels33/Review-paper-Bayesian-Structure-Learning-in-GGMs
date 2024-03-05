

source("metric_functions.R")

###################################
######supporting  functions ############
###################################
Soft_thres<-function(a,b,lambda){
  return(sign(-b/a)*max(abs(-b/a)-lambda/a,0))
}
Cor_desc<-function(Theta_11_inverse,W_22,Theta_12,S_12,v0,v1,P_12,n){
  j=nrow(Theta_11_inverse)
  o=1
  Theta_12_1=10000
  Theta_12_2=Theta_12
  while(max(abs(Theta_12_1-Theta_12))>10^-3&o<=1000){
    o=o+1
    Theta_12_1=Theta_12
    for(i in 1:j){
      b=n/2*((Theta_11_inverse[i,]%*%Theta_12)*W_22-
               (Theta_11_inverse[i,i]%*%Theta_12[i])*W_22+S_12[i,])
      a=n/2*Theta_11_inverse[i,i]*W_22
      lambda=P_12[i]/v1+(1-P_12[i])/v0
      Theta_12[i]=Soft_thres(a=a,b=b,lambda=lambda)
    }
    if(is.nan(sum(Theta_12))|max(abs(Theta_12))>10){
      return(Theta_12_2)
    }
  }
  if(o==1001){
    return(Theta_12_2)
  }else{
    return(Theta_12) 
  }
}
Bagus<-function(S,n,v0,v1,maxiter,eta,tau,K_start,P_start){ #edit_LV
  
  #####Initialization######
  loss=NULL
  eigen=eigen1=NULL
  p_n=nrow(S)
  #W=S
  #diag(W)=diag(W)+2/n*tau
  #W=diag(diag(W))
  W=diag(1,p_n)
  #Theta=diag(1,p_n)
  Theta=K_start
  P = P_start
  o1=0
  
  ##iteration
  while(o1<maxiter){ 
    o1=o1+1
    
    ####E-step######
    ###Determine adaptive shrinkage parameter######
    for(i in 1:p_n){
      if(i==1){
        id=(i+1):p_n
      }else if(i!=p_n){
        id=c(1:(i-1),(i+1):p_n)
      }else{
        id=1:(i-1)
      }
      W_11=W[id,id]
      W_12=matrix(W[i,id])
      W_22=W[i,i]
      
      Theta_12=matrix(Theta[i,id])
      Theta_11=Theta[id,id]
      Theta_22=Theta[i,i]
      
      P_12=matrix(P[i,id])
      
      ####M-step#####
      S_12=matrix(S[i,id])
      
      #W_22=S[i,i]+2/n*tau
      
      Theta_11_inverse=W_11-W_12%*%t(W_12)/W_22
      #W_11-W_12%*%t(W_12)/W_22
      #solve(Theta_11)
      
      W_22=S[i,i]+2/n*tau
      ###Coordinate descent to update theta_12#####
      Theta_12=Cor_desc(Theta_11_inverse,W_22,Theta_12,S_12,v0,v1,P_12,n)
      ###Update theta_22, W#####
      
      Theta_22=1/(W_22)+t(Theta_12)%*%Theta_11_inverse%*%Theta_12
      
      Theta[i,id]=Theta[id,i]=Theta_12
      Theta[i,i]=Theta_22
      
      tmp=as.vector(Theta_22-t(Theta_12)%*%Theta_11_inverse%*%Theta_12)
      tmp1=Theta_11_inverse%*%Theta_12
      W_11=Theta_11_inverse+
        tcrossprod(tmp1)/tmp
      W_12=-tmp1/tmp
      
      W[id,id]=W_11
      W[i,id]=W[id,i]=W_12
      W[i,i]=W_22
    }
    
    P=1/(1+v1/v0*exp(-abs(Theta)/(v0)+abs(Theta)/(v1))*(1-eta)/eta)
  }
  
  return(list(Theta=Theta,P=P,tau=tau,p=p,W=W))
}
BIC_Bagus<-function(Theta,S,P,n){
  p_n=nrow(Theta)
  P1=diag(1,p_n)
  P1[P>0.5]=1
  if(min(eigen(Theta*P1)$values)>0){
    n*(sum(diag(S%*%(Theta*P1)))-log(det(Theta*P1)))+(log(n))*(sum(P>0.5)/2)
  }else{
    Inf
  }
}
Tune_Bagus=function(v0,tau,S,n,p,eta){
  pb <- txtProgressBar(min = 0, max = length(v0)*4, style = 3)
  bic=NULL
  for(m in 1:length(eta)){
    for(i in 1:length(v0)){
      v1=v0[i]*c(1.5,3,5,10)
      for(j in 1:length(v1)){
        # for(k in 1:length(tau)){
        ##Bayes EM####
        #v1=sqrt(n/(log(p_n)))/5
        #v0=sqrst(n/(p_n*log(p_n)))/5
        
        #tau=sqrt((log(p_n)/n))
        w=1
        l=1
        maxiter=30
        if (p > 250){maxiter = 5}
        
        #initialize precision matrix and edge inclusion prob matrix
        K_start = diag(1,p)
        P_start = P=1/(1+v1[j]/v0[i]*exp(-abs(K_start)/(v0[i])+abs(K_start)/(v1[j]))*(1-eta[m])/eta[m])
        
        result1<-Bagus(S,n,v0=v0[i],v1=v1[j],maxiter,eta=eta[m],tau=v0[i],K_start,P_start)
        
        bic=rbind(bic,list(v0=v0[i],v1=v1[j],tau=v0[i],eta=eta[m],BIC=BIC_Bagus(result1$Theta,S,result1$P,n)))
        
        setTxtProgressBar(pb, (m-1)*length(v0)*length(v1)+(i-1)*length(v1)+j)
        
      }
      #} 
    }
    
  }
  
  close(pb)
  return(bic[which.min(bic[,5]),])
}

###################################
######BAGUS  function ############
###################################

#the function bagus_solve takes as input:
#- data from BDgraph. This "data" type contains data, a true graph 
#- iter_vec_thin. The iteration at which it needs to produce output

#It outputs:
#- the AUC at every iteration number indicated in iter_vec_thin
#- the p+ at every iteration number indicated in iter_vec_thin
#- the p- at every iteration number indicated in iter_vec_thin
#- the running time at every iteration number indicated in iter_vec_thin
#- the time to compute the tuning parameters

BAGUS_solve = function(data,iter_vec_thin){
  
  #obtain true graph
  G_true = data$G
  response = G_true[upper.tri(G_true)]
  
  #obtain amount of observations and variables
  n = nrow(data$data)
  p = ncol(data$data)
  density = data$density
  graph = data$graph
  rep = data$rep
  algorithm_name = "BAGUS"
  
  #obtain sample covariance matrix
  S<-cov(data$data)  
  
  ##############
  ####tuning####
  ############## 
  v0=tau=sqrt(1/(n*log(p)))*c(0.4,2,4,20)
  eta = 0.5
  start = proc.time()
  Tune=Tune_Bagus(v0,tau,S,n,p,eta)
  time_tuning = as.numeric( ( proc.time() - start)[ 3 ] ) #end and save tim
  
  v0_t=Tune$v0
  v1_t=Tune$v1
  tau_t=Tune$tau
  
  ###########################
  #####model selection ######
  ###########################
  
  #initialize output vectors
  len = length(iter_vec_thin)
  AUC_vec = rep(0,len)
  pplus_vec = rep(0,len)
  pmin_vec = rep(0,len)
  time_vec = rep(0,len)
  
  #initialize precision matrix and edge inclusion prob matrix
  K_start = diag(1,p)
  P_start = P=1/(1+v1_t/v0_t*exp(-abs(K_start)/(v0_t)+abs(K_start)/(v1_t))*(1-eta)/eta)
  
  for (i in 1:len){
    if (i==1){
      maxiter = iter_vec_thin[1]
    }
    if (i>1){
      maxiter = iter_vec_thin[i] - iter_vec_thin[i-1] 
    }
    
    #run maxiter iterations of the BAGUS algorithm
    start = proc.time()
    output = Bagus(S,n,v0_t,v1_t,maxiter,eta,tau_t,K_start,P_start)
    runtime = as.numeric( ( proc.time() - start)[ 3 ] ) #end and save time
    
    #retrieve output algorithm  
    K_start = output$Theta
    P_start = output$P
    predictor = P_start[upper.tri(P_start)]
    
    #calculate metrics
    AUC = calc_AUC_ROC(predictor=predictor,response=response)
    obj = calc_pplus_pmin(predictor=predictor,response=response)
    pplus = obj$p_plus
    pmin = obj$p_min
    
    #save metrics
    AUC_vec[i] = AUC
    pplus_vec[i] = pplus
    pmin_vec[i] = pmin
    time_vec[i] = runtime
    
    #save output per iteration in Rdata file (only for large scale problems)
    if (p > 250){
      output_per_iter = list()
      output_per_iter$AUC = AUC
      output_per_iter$pplus = pplus
      output_per_iter$pmin = pmin
      output_per_iter$time = sum(time_vec)
      output_per_iter$iter = iter_vec_thin[i]
      filename = paste0("results_p",p,"_n",n,"_",graph,"_dens",density,"_",algorithm_name,"_rep",rep,"_iter",i,".Rdata")
      save(output_per_iter,file = filename)
    }
    
  }
  time_vec = cumsum(time_vec)
  
  return(list(AUC_vec=AUC_vec,pplus_vec=pplus_vec,pmin_vec=pmin_vec,time_vec=time_vec,iter_vec_thin=iter_vec_thin,time_tuning = time_tuning,predictor = predictor))
  
  
}



