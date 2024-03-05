#This Rfile contains the code to reproduce the results of Section 5 ("Case study") of the article:
#"Bayesian Structure Learning in Undirected Gaussian Graphical Models: Literature Review with Empirical Comparison"

#This file contains 9 chapters
#Chapter 1. code to create Figure 2
#Chapter 2. code to create Figure 3
#Chapter 3. code to create Figure 4 (left)
#Chapter 4. code to create Figure 4 (right)
#Chapter 5. code to create Figure 5 
#Chapter 6. code to create Figure 6 
#Chapter 7. code to create Figure 7
#Chapter 8. code to find similarities between our results and results in existing literature
#Chapter 9. code to create the table S3 in the supplementary material

#necessary libraries
library(BDgraph)
library(ssgraph)
library(matrixStats)
library(huge)
library(Matrix)

#Since, we will be saving and loading Rdata files, please make sure to set an appropirate working directory before running the code. 

#Chapter 1
#############################################################################
#####create histograms of univariate distributions of 6 random genes ########
#############################################################################

#load data
data(geneExpression) 
gene_data = geneExpression 

#select 6 random genes
set.seed(1)
random_select = sample.int(100,6,replace=FALSE)
random_select = sort(random_select,decreasing=FALSE)

#create titles for each of the 6 plots
gene_names = colnames(gene_data)
title_names = gene_names[random_select]
title_nr = c("Gene 1","Gene 34","Gene 39","Gene 43","Gene 68","Gene 87")
title = rep("",6)
for (i in 1:6){
  title[i] = paste0(title_nr[i],": ",title_names[i],sep="")
}

#set other parameters for the plots
par(cex=0.7,cex.lab=1.3, mai=c(0.1,0.6,0.4,0.1))
xlim = c(5,16)
ylim = c(0,25)
breaks = seq(5.25,15.75,0.5)

#plot the 6 univariate histograms
i = 1
par(fig = c(0.05,0.34,0.51,0.95))  
hist(gene_data[,random_select[i]],breaks=breaks,ylim=ylim,xlim=xlim,main=title[i],ylab="Frequency",xlab="")
i = 2
par(fig = c(0.36,0.64,0.51,0.95),new=TRUE)  
hist(gene_data[,random_select[i]],breaks=breaks,ylim=ylim,xlim=xlim,main=title[i],ylab="",xlab="")
i = 3
par(fig = c(0.66,0.95,0.51,0.95),new=TRUE)  
hist(gene_data[,random_select[i]],breaks=breaks,ylim=ylim,xlim=xlim,main=title[i],ylab="",xlab="")
i = 4
par(fig = c(0.05,0.34,0.05,0.49),new=TRUE)  
hist(gene_data[,random_select[i]],breaks=breaks,ylim=ylim,xlim=xlim,main=title[i],ylab="Frequency",xlab="") 
i = 5
par(fig = c(0.36,0.64,0.05,0.49),new=TRUE)  
hist(gene_data[,random_select[i]],breaks=breaks,ylim=ylim,xlim=xlim,main=title[i],ylab="",xlab="")
i = 6
par(fig = c(0.66,0.95,0.05,0.49),new=TRUE)  
hist(gene_data[,random_select[i]],breaks=breaks,ylim=ylim,xlim=xlim,main=title[i],ylab="",xlab="")

#Chapter 2
#############################################################
################create edge convergence statistics############
#############################################################

#This Chapter creates the edge convergence plots (Figure 3) for the SS and BD-A algorithms. It has three supporting functions:
#1. BDA_solve
#2. SS_solve
#3. edge_per_state
# We will first show these functions and then continue with the code to produce Figure 3

#Function 1. BDA_solve
#This function takes as input:
#- the covariance matrix produced by huge.npn
#- the iterations at which it needs to output the #edges in the Markov state, 
#- the initital graph
#- the amount of observations
#It solves the data using the BDA algorithm outputs the #edges per Markov state at every iteration specified.
BDA_solve = function(data,iter_vec_thin,MCMCstart,n){
  
  
  #obtain amount of variables
  p = dim(data)[1]
  
  #parameters
  burnin = 0
  jump = 1
  save = FALSE
  cores = 1
  g.prior = 0.2
  verbose = TRUE
  
  #run the algorithm, produce output at every iteration in the iter_vec_thin vector
  len = length(iter_vec_thin)
  edge_vec = rep(0,len)
  
  for (j in 1:len){
    
    #determine amount of iterations
    olditer = 0
    if (j > 1){
      olditer = iter_vec_thin[j-1]
    }
    newiter = iter_vec_thin[j]
    iter = newiter-olditer
    
    #run bdgraph for maxiter iterations
    sample_bd  = bdgraph( data = data, n =n,algorithm = "bdmcmc", iter = iter, burnin = burnin, jump = jump, save = save,cores=cores,g.start=MCMCstart,g.prior=g.prior,verbose=verbose)  
    
    #save metrics
    edge_vec[j] = sum(sample_bd$last_graph)/2
    MCMCstart = sample_bd
  }
  
  return(list(edge_vec=edge_vec))
}  

#Function 2. SS_solve
#This function takes as input:
#- the covariance matrix produced by huge.npn
#- the iterations at which it needs to output the #edges in the Markov state, 
#- the initital graph
#- the amount of observations
#It solves the data using the SS algorithm outputs the #edges per Markov state at every iteration specified.
SS_solve= function(data,iter_vec_thin,MCMCstart,n){
  
  
  #obtain amount of variables
  p = dim(data)[1]
  
  #parameters
  burnin = 0
  var1 = 0.02
  var2 = 2
  lambda = 1
  save = TRUE
  cores = 1
  sig.start = NULL
  if (p > n){sig.start = diag(1,p)}
  g.prior = 0.2
  verbose = FALSE
  
  #run the algorithm, produce output at every iteration in the iter_vec_thin vector
  len = length(iter_vec_thin)
  edge_vec = rep(0,len)
  
  #run bdgraph for maxiter iterations
  iter= max(iter_vec_thin)
  sample_ss  = ssgraph( data = data, n=n,iter = iter, burnin = burnin, var1 = var1, var2=var2,lambda=lambda,g.start=MCMCstart,sig.start=sig.start,save=save,cores=cores,g.prior=g.prior,verbose=verbose)
  
  #save output
  all_graphs = sample_ss$all_graphs[iter_vec_thin]
  sample_graphs = sample_ss$sample_graphs[all_graphs]
  
  for (j in 1:length(sample_graphs)){
    which_edge = which(unlist(strsplit(as.character(sample_graphs[j]),"")) == 1) 
    edge_vec[j]=length(which_edge)
  }

  return(list(edge_vec=edge_vec))
} 

#Function 3. edge_per_state
#This function that takes as input:
#- the algorithm name (either "SS" or "BDA")
#- the data as a covariance matrix 
#- the densities of the Markov Chain initialization
#It outputs a matrix with #edges per iteration (column) per density (row)
edge_per_state =function(algorithm,data,densities,n){
  
  #calculate dimension and amount of links
  p = dim(data)[1]
  qp = p*(p-1)/2
  
  #set iterations
  if (algorithm == "SS"){iter_vec_thin = c(1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100)}
  if (algorithm == "BDA"){iter_vec_thin = c(10,20,30,40,50,100,200,300,400,500,750,1000,1500,2000,2500,3000,4000,5000,7500,10000,15000,20000,25000,30000)}
  
  #create empty output matrix
  edge_matrix = matrix(0,length(densities),length(iter_vec_thin)+1)
  
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
  
    if (algorithm == "SS"){out = SS_solve(data=data,iter_vec_thin=iter_vec_thin,MCMCstart=init_graph,n)}
    if (algorithm == "BDA"){out = BDA_solve(data=data,iter_vec_thin=iter_vec_thin,MCMCstart=init_graph,n)}
    edge_vec = out$edge_vec
    edge_vec = c(sum(init_graph)/2,edge_vec)
    edge_matrix[q,]=edge_vec
  }
  
  iter_vec_thin = c(0,iter_vec_thin)
  
  return(list(edge_matrix=edge_matrix,iter_vec_thin=iter_vec_thin))
}

####the code to produce Figure 3####

#solve data and save in R.data file 
densities = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)
new_data = huge.npn(x=geneExpression,npn.func="skeptic") #new data is an estimate of the sample covariance matrix of the transformed sample
data = new_data
n = 60
for (algorithm in c("SS","BDA"){
  output = edge_per_state(algorithm=algorithm,data=data,densities=densities,n=n)
  filename = paste0("edge_convergence_",algorithm,".Rdata")
  save(output,file = filename)
}

#load data SS
algorithm = "SS"
filename = paste0("edge_convergence_",algorithm,".Rdata")
load(file=filename)
iter_vec_thin_ss = output$iter_vec_thin
edge_matrix_ss = output$edge_matrix

#load data (BD)
algorithm = "BDA"
filename = paste0("edge_convergence_",algorithm,".Rdata")
load(file=filename)
iter_vec_thin_bd = output$iter_vec_thin
edge_matrix_bd = output$edge_matrix

#########make plots ############
par(mfrow = c(1, 2))

par(cex.lab=1.1)
ylim = c(0,p*(p-1)/2)
xlab = "MCMC iterations"
ylab = "Amount of edges in state"

#plot ss
xlim = c(0,max(iter_vec_thin_ss))
plot(NA,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,main="SS algorithm")
len = dim(edge_matrix_ss)[1]
for (i in 1:len){
  points(x=iter_vec_thin_ss,y=edge_matrix_ss[i,],type="l")
}

#plot bd
xlim = c(0,max(iter_vec_thin_bd))
plot(NA,xlim=xlim,ylim=ylim,xlab=xlab,ylab="",main="BD-A algorithm")
len = dim(edge_matrix_bd)[1]
for (i in 1:len){
  points(x=iter_vec_thin_bd,y=edge_matrix_bd[i,],type="l")
}

#Chapter 3
#######################################################################
#############create links convergence plots ss algorithm###############
#######################################################################

#this Chapter contains the code to produce the link-convergence plot of the SS algorithm. Figure 4 (left).
#It contains one function (posterior_means) and the main code to create the plot.

#The posterior_means function takes as input an ssgraph or bdgraph object and outputs a (#edges x #iterations) matrix 
#where the ijth elements corresponds to the inclusion probability of edge i at iteration j
posterior_means = function(data){
  
  output = data
  all_graphs = output$all_graphs
  all_weights = output$all_weights
  sample_graphs = output$sample_graphs
  iterations = length(all_graphs)
  qp = length(unlist(strsplit(as.character(sample_graphs[1]),"")))
  
  mean_post = matrix(0, qp, iterations)
  sum_weights = 0
  sum_weights_vec = rep(0,qp)
  sum_weights_sq_vec = rep(0,qp)
  
  for (i in 1:iterations){
    if (i==1000){
      print(i)
    }
    if (i==5000){
      print(i)
    }
    which_edge = which(unlist(strsplit(as.character(sample_graphs[all_graphs[i]]),"")) == 1)
    sum_weights_vec[which_edge] = sum_weights_vec[which_edge]+all_weights[i]
    sum_weights = sum_weights + all_weights[i]
    mean_post[,i] = sum_weights_vec/sum_weights
  }
  
  return(list(mean_post=mean_post))
  
}

#Run the SS algorithm 

#Set the general parameters
jump = 1
save = TRUE
cores = 1
MCMCstart = "empty"
g.prior = 0.2
verbose = FALSE
data = new_data
n = 60
p = 100

#Set the SS specific parameters
var1 = 0.02
var2 = 2
lambda = 1 
iter = 2000
burnin = 50
sig.start = NULL
if (p > n){sig.start = diag(1,p)}

#run the SS algorithm
sample_ss  = ssgraph( data = data, n=n,iter = iter, burnin = burnin, var1 = var1, var2=var2,lambda=lambda,g.start=MCMCstart,sig.start=sig.start,save=save,cores=cores,g.prior=g.prior,verbose=verbose)

#calculate posterior means per iterations, and save in matrix
out_ss = posterior_means(sample_ss)
mean_post_ss = out_ss$mean_post

#############and select links to be plotted#########
links_needed = 10
real_iter_ss = iter-burnin
mean_post_ss = output_links_ss$mean_post_ss
plinks_final = mean_post_ss[,real_iter_ss]

#find the link with the highest edge inclusion prob
max_plinks = max(plinks_final)
index_max = which(plinks_final ==max_plinks)[1]

#find the link with the lowest edge inclusion prob
min_plinks = min(plinks_final)
index_min = which(plinks_final ==min_plinks)[1]

#find other links at random
links_needed_2 = links_needed - 2
qp = p*(p-1)/2
index_others = round(runif(links_needed_2,1,qp))

#all selected links
index_links = c(index_max,index_min,index_others)

##########plot convergence#################
plot(NA,xlim=c(1,real_iter_ss),xlab="MCMC iterations after burnin",ylab="Edge inclusion prob.",ylim=c(0,1),main="SS algorithm")
for (i in 1:links_needed){
  link = index_links[i]
  points(x=1:real_iter_ss,y=mean_post_ss[link,],type="l")
}

#Chapter 4
#######################################################################
#############create links convergence plots BDA algorithm###############
#######################################################################

#this Chapter contains the code to produce the link-convergence plot of the BDA algorithm. Figure 4 (right).
#It uses the posterior_means function given in Chapter 3 above. 

#Run the BD algorithm 
#set parameters
jump = 1
save = TRUE
cores = 1
MCMCstart = "empty"
g.prior = 0.2
verbose = FALSE
data = new_data
n = 60
p = 100
iter = 500000
burnin = 25000

#run code (these two lines can only be run on a cluster computer due to their memory requirement)
set.seed(3)
sample_bd  = bdgraph( data = data, n=n, algorithm = "bdmcmc", iter = iter, burnin = burnin, jump = jump, save = save,cores=cores,g.start=MCMCstart,g.prior=g.prior,verbose=verbose)  

#obtain posterior mean matrix
out = posterior_means(sample_bd)
mean_post = out$mean_post

#select links to be plotted ###########
links_needed = 10
iter_real = iter-burnin
plinks_final = mean_post[,iter_real]

#find the link with the highest edge inclusion prob
max_plinks = max(plinks_final)
index_max = which(plinks_final ==max_plinks)[1]

#find the link with the lowest edge inclusion prob
min_plinks = min(plinks_final)
index_min = which(plinks_final ==min_plinks)[1]

#find the other links (one in every interval)
index_links = c(index_max,index_min)
intervals = seq(0+1/links_needed,1-1/links_needed,1/links_needed)
for (i in 1:(links_needed-2)){
  lower = intervals[i]
  upper = intervals[i+1]
  index = which(plinks_final > lower & plinks_final < upper)[1]
  if (length(index)>0){ #if an index is found within this interval, then add it to the list
    index_links = c(index_links,index)
  }
}

#save only data for selected links to save memory
mean_post_links = mean_post[index_links,]

#select only one in 100 iterations to conserve memory
iter_vec_thin = c(1,seq(100,iter_real,100))
mean_post_links_iter = mean_post_links[,iter_vec_thin]

#save thinned matrix in .Rdata file for future reference
output = list(mean_post_links_iter=mean_post_links_iter)
filename = ("post_mean_links_iter_matrix.Rdata") 
save(output,file = filename)

#load .Rdata file
filename = ("post_mean_links_iter_matrix.Rdata")
load(file=filename)
mean_post_links_iter = output$mean_post_links_iter

#plot link convergence
plot(NA,xlim=c(1,iter_real),xlab="MCMC iterations after burnin",ylab="",ylim=c(0,1),main="BD algorithm")
links_needed = dim(mean_post_links_iter)[1]
for (i in 1:links_needed){
  points(x=iter_vec_thin,y=mean_post_links_iter[i,],type="l")
}

#Chapter 5
##########################################################
#############prior sensitivity ss############################
##########################################################

#this Chapter contains the code to produce the prior sensitivity plot of the SS algorithm (Figure 5).

#general parameters
jump = 1
save = FALSE
cores = 1
MCMCstart = "empty"
verbose = FALSE
data = new_data
n = 60
p = 100

#parameters ss algorithm 
var1 = 0.02
var2 = 2
lambda = 1 
iter = 500
burnin = 50
sig.start = NULL
if (p > n){sig.start = diag(1,p)}

#initialize input and output vectors
prior_vec = c(0.01,0.2,0.5)
nr_priors = length(prior_vec)
plinks_matrix = matrix(0,nr_priors,p*(p-1)/2)

#run the ss algorithm for every prior in prior_vec, and save edge incl. prob in plinks_matrix
count = 0
for (g.prior in prior_vec){
  count = count + 1
  sample_ss  = ssgraph( data = data, n=n,iter = iter, burnin = burnin, var1 = var1, var2=var2,lambda=lambda,g.start=MCMCstart,sig.start=sig.start,save=save,cores=cores,g.prior=g.prior,verbose=verbose)
  plinks = sample_ss$p_links
  plinks_vec = plinks[upper.tri(plinks)]
  plinks_matrix[count,] = plinks_vec
}

#save plinks_matrix for later use
output = list(plinks_matrix=plinks_matrix,prior_vec=prior_vec)
filename = ("ss_plinks_per_prior.Rdata") 
save(output,file = filename) 

#load plinks_matrix 
filename = ("ss_plinks_per_prior.Rdata") 
load(file=filename)
prior_vec = output$prior_vec
plinks_matrix = output$plinks_matrix

#plot results
par(mfrow = c(1, 3))
for (i in 1:nr_priors){
  boxplot(plinks_matrix[i,],ylim=c(0,1),main=paste0(100*prior_vec[i],"% prior sparsity"),boxwex=1.5)
}

#Chapter 6
##########################################################
#############prior sensitivity bd############################
##########################################################

#this Chapter contains the code to produce the prior sensitivity plot of the BDA algorithm (Figure 6).

#general parameters
jump = 1
save = FALSE
cores = 1
MCMCstart = "empty"
verbose = FALSE
data = new_data
n = 60
p = 100
iter = 100000
burnin = 25000
prior_vec = c(0.01,0.2,0.5)
nr_priors = length(prior_vec)

#run the bd algorithm for every prior in prior_vec, and save edge incl. prob in plinks_matrix in R.data files for future reference
count = 0
for (g.prior in prior_vec){
  count = count + 1
  sample_bd  = bdgraph( data = data, n=n, algorithm = "bdmcmc", iter = iter, burnin = burnin, jump = jump, save = save,cores=cores,g.start=MCMCstart,g.prior=g.prior,verbose=verbose)  
  plinks = sample_bd$p_links
  plinks_vec = plinks[upper.tri(plinks)]
  output = list(plinks_matrix=plinks_matrix,prior_vec=prior_vec)
  filename = paste0("bd_plinks_prior",g.prior,".Rdata") 
  save(output,file = filename)   
}

#load plinks_matrix in .Rdata files
plinks_matrix = matrix(0,nr_priors,p*(p-1)/2)
count = 0
for (prior in prior_vec){
  filename = paste0("bd_plinks_prior",prior,".Rdata") 
  load(file = filename)
  plinks_vec = output$plinks_vec
  count = count +1
  plinks_matrix[count,] = plinks_vec
}

#plot results
par(mfrow = c(1, 3))
for (i in 1:nr_priors){
  boxplot(plinks_matrix[i,],ylim=c(0,1),main=paste0(100*prior_vec[i],"% prior sparsity"),boxwex=1.5)
}

#Chapter 7
#######################################################################################
#############code to create gene network########
#######################################################################################

#this Chapter contains the code to produce the gene network (Figure 7).

#Run the BD algorithm 
#set parameters
jump = 1
save = FALSE
cores = 1
MCMCstart = "empty"
g.prior = 0.2
verbose = FALSE
data = new_data
n = 60
p = 100
iter = 500000
burnin = 25000

#run code (these two lines can only be run on a cluster computer due to their memory requirement)
set.seed(3)
sample_bd  = bdgraph( data = data, n=n, algorithm = "bdmcmc", iter = iter, burnin = burnin, jump = jump, save = save,cores=cores,g.start=MCMCstart,g.prior=g.prior,verbose=verbose)  

#obtain plinks matrix (the matrix containing the edge inclusion probabilities)
plinks = output$sample_bd$p_links

#fill the lower triangle of the matrix
p = 100
count = 0
for (i in 2:p){ #fill lower matrix
  for (j in 1:(i-1)){
    count = count + 1
    plinks[i,j] = plinks[j,i]
  } 
}
#determine the index of the 65 links with highest edge incl. prob.
library(igraph)
nr_links_show = 65
plinks_sort = sort(plinks,decreasing=TRUE)
cutoff = 0.5*plinks_sort[nr_links_show*2]+0.5*plinks_sort[nr_links_show*2+1]
index = which(plinks>cutoff)

#create adjacency matrix of an undirected graph containing only those 65 links
G = matrix(0,100,100)
G[index] = 1
G_ig = graph_from_adjacency_matrix(adjmatrix=G,mode="undirected")

#obtain vector of weights of these 65 links
new_matrix = G*plinks
weights = new_matrix[lower.tri(new_matrix)]
weights = weights[which(weights>0)]

#normalize this weight vector
weights_norm = (weights-min(weights))/(max(weights)-min(weights))

#create partial correlation matrix C
K = output$sample_bd$K_hat
K_diag_vec = 1/sqrt(diag(K))
K_diag = diag(K_diag_vec)
C = -K_diag%*%K%*%K_diag

#obtain vector of part. corr. of all 65 links
G = matrix(0,100,100)
G[index] = 1
new_matrix = G*C
part_cor_vec = new_matrix[lower.tri(new_matrix)]
part_cor_vec = part_cor_vec[which(part_cor_vec!=0)]

#make edges blue that have neg part. correlation, and red that have pos. partial correlation
edge.color = rep("red",length(part_cor_vec))
edge.color[which(part_cor_vec<0)] = "blue"

#save layout and labels 
V(G_ig)$label = 1:100
layout_bd = layout_nicely(G_ig,dim=2)

#remove vertices that are unconnected
isolated = which(degree(G_ig)==0)
layout_isolated = layout_bd[-isolated,]
G_ig_2 = delete.vertices(G_ig, isolated)

#plot network of edge inclusion probabilities, edge width is edge inclusion prob, edge color is sign of partial correlation
plot.igraph(G_ig_2,mode="undirected",layout=layout_isolated,diag=FALSE,edge.width=5*weights_norm,vertex.size=0.1,edge.color = edge.color,vertex.label.color ="black",vertex.label.dist=0.5)

#Chapter 8
#################################################################################################################################
###########find similarities between our network and networks in literature ####################################
#################################################################################################################################

#####find similarities with the K-Horseshoe paper ####
plinks = output$sample_bd$p_links

#fill lower part of matrix
count = 0
for (i in 2:p){ #fill lower matrix
  for (j in 1:(i-1)){
    count = count + 1
    plinks[i,j] = plinks[j,i]
  } 
}

#find index of all links that have edge incl. prob higher than 0.5
cutoff = 0.5
index = which(plinks>cutoff)

#create adjacency matrix of those links
G = matrix(0,100,100)
G[index] = 1

#create vector "connectivity" with every entry the degree of a node
connectivity = colSums(G)

#print most connected genes
max_conn = max(connectivity)
gene_names = rownames(plinks)
for (i in max_conn:5){
  print("The following genes have degree")
  print(i)
  index = which(connectivity==i)
  print(gene_names[index])
}

#####find similarities with the Bhadra and Mallick 2013 ####

#the below tree structure appear in Figure 7. We check here if they also appear in Bhadra and Mallick 2013
gene_names[c(86,3,5,23,57,27)]
gene_names[c(90,63,53,47,97,58,26,95,51)]
gene_names[c(73,40)]
gene_names[c(38,13,18,4,8,6,9,10,16)]
gene_names[c(25,93,1,2,12,17)]


#Chapter 9
#################################################################################################################################
###########Create table S3 in supplemtary material ####################################
##############################################################################################

install.packages("lazyWeave")
library(lazyWeave)
table_r = matrix(0,nrow=100,ncol=2)
table_r[,1] = 1:100
table_r[,2] = gene_names
table_latex = lazy.matrix(table_r,align="left")







  
