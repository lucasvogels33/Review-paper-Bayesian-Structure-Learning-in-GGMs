#this file contains three functions to calculate the metrics:
#i) calc_pplus_pmin, ii)calc_AUC_ROC and iii)calc_TPR_FPR. 

# All three functions take as input a "response" and "predictor" vector.
# Every element of the response vector is either zero or one. The response 
# vector is the true value of the graph. A one for a present edge. A zero 
# for an absent edge. The predictor contains the edge inclusion probabilities.

#i) calc_pplus_min outputs the value of Pr+ and Pr- for the given predictor and response
#ii) calc_AUC_ROC outputs the AUC value for the given predictor and response
#iii) calc_TPR_FPR outputs the true positive rate and false positive rate for the given predictor and response



calc_pplus_pmin = function(predictor,response){
  if (length(response) != length(predictor)) {
    stop("response and predictor vector must be of same length")
  }
  
  #calculate the weighed CE and weighed MSE
  ones_index = which(response==1)
  zeroes_index = which(response==0)
  
  predictor_ones =  predictor[ones_index]
  predictor_zeroes = predictor[zeroes_index]
  
  p_plus = mean(predictor_ones)
  p_min = mean(predictor_zeroes)
  
  return(list(p_plus=p_plus,p_min=p_min))
}

calc_AUC_ROC = function(predictor,response){
  
  if (length(response) != length(predictor)) {
    stop("response and predictor vector must be of same length")
  }
  
  #order the vectors so that the predictor is increasing
  predictor.order = order(predictor,decreasing=FALSE)
  predictor.sorted = predictor[predictor.order]
  response.sorted = response[predictor.order]
  
  #determine amount of zeroes and ones
  ones = sum(response)
  zeroes = length(response)-ones
  
  #if there are duplicates
  if (sum(duplicated(predictor.sorted))>0){
    #create a vector with one index for every group of duplicates
    dup_index = cumsum(duplicated(predictor.sorted)==0)
    
    #create a vector sum_vec that sums the true positives in each group of duplicates   
    df <- data.frame(duplicates=dup_index,response.sorted=response.sorted)
    sum_vec = aggregate(response.sorted ~ duplicates, data=df, sum)[,2]
    
    #create a vector that averages the maximum amount of false positives of the current group with the previous group
    fp = cumsum(response.sorted==0)
    df <- data.frame(duplicates=dup_index,fp=fp)
    max_vec = aggregate(fp ~ duplicates, data=df, max)[,2]
    top = c(0,max_vec)
    bottom = c(max_vec,0)
    average_vec = head((top+bottom)/2,-1)
    
    #AUC is the dot product of the two vectors divided by the normalizing constant
    AUC = (sum_vec%*%average_vec)/(ones*zeroes)
    
  }
  
  #if there are no duplicates
  if (sum(duplicated(predictor.sorted))==0){
    fp = cumsum(response.sorted==0)
    AUC = sum(fp * response.sorted)
    AUC = AUC/(zeroes*ones)  
    
  }
  
  return(AUC)
}

calc_TPR_FPR = function(predictor,response){
  ones = sum(response)
  zeroes = length(response)-ones
  
  TPR = sum(predictor[which(response==1)])/ones
  FPR = sum(predictor[which(response==0)])/zeroes
  
  return(list(TPR=TPR,FPR=FPR))
  
}