

# set working directory to folder that contains the files created by:
#- edge_convergence.R
#- links_convergence.R
#- prior_sensitivity.R
#- networks.R

############################################################
######################edge convergence######################
############################################################
dev.off()
par(mfrow = c(1, 2))

algorithm = "SS"
density_vec = seq(0,1,0.1)
filename = paste0("edge_convergence_",algorithm,"_density",density_vec[1],".Rdata")
load(file=filename)
total_iter = tail(output$all_iter_vec,1)
xlim = c(0,total_iter)
plot(NA,xlim=c(0,20),ylim=c(0,5000),main=paste0("SS-O algorithm"),xlab="MCMC iterations",ylab = "Number of edges in state")
for (density in density_vec){
  filename = paste0("edge_convergence_",algorithm,"_density",density,".Rdata")
  load(file=filename)
  points(x=output$all_iter_vec,y=output$edge_matrix,type="l",lw=1.5)
}

algorithm = "BDA"
density_vec = seq(0,1,0.1)
filename = paste0("edge_convergence_",algorithm,"_density",density_vec[1],".Rdata")
load(file=filename)
total_iter = tail(output$all_iter_vec,1)
xlim = c(0,total_iter)
plot(NA,xlim=c(0,30000),ylim=c(0,5000),main=paste0("BD-A algorithm"),xlab="MCMC iterations",ylab = "Number of edges in state")
for (density in density_vec){
  filename = paste0("edge_convergence_",algorithm,"_density",density,".Rdata")
  load(file=filename)
  points(x=output$all_iter_vec,y=output$edge_matrix,type="l",lw=1.5)
}

#######################################################################
#############create links convergence plots ###########################
#######################################################################
options(scipen=999)
dev.off()
par(mfrow = c(1, 2))
links_needed = 10

#create plot for SS algorithm
algorithm = "SS"
filename = paste0("plinks_convergence_",algorithm,".Rdata")
load(file=filename)
plinks_iter_matrix = output$plinks_iter_matrix
iter_vec_thin = output$iter_vec_thin
plinks_final = plinks_iter_matrix[,length(iter_vec_thin)]

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

#plot convergence
max_iter = iter_vec_thin[length(iter_vec_thin)]
plot(NA,xlim=c(1,max_iter),xlab="MCMC iterations after burnin",ylab="Edge inclusion prob.",ylim=c(0,1),main="SS-O algorithm")
for (i in 1:links_needed){
  index = index_links[i]
  points(x=iter_vec_thin,y=plinks_iter_matrix[index,],type="l")
}

#create plot for BDA algorithm
algorithm = "BDA"
filename = paste0("plinks_convergence_",algorithm,".Rdata")
load(file=filename)
plinks_iter_matrix = output$plinks_iter_matrix
iter_vec_thin = output$iter_vec_thin
plinks_final = plinks_iter_matrix[,length(iter_vec_thin)]

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

#plot convergence
max_iter = iter_vec_thin[length(iter_vec_thin)]
plot(NA,xlim=c(1,max_iter),xlab="MCMC iterations after burnin",ylab="Edge inclusion prob.",ylim=c(0,1),main="BD-A algorithm")
for (i in 1:links_needed){
  index = index_links[i]
  points(x=iter_vec_thin,y=plinks_iter_matrix[index,],type="l")
}


##########################################################
#############prior sensitivity############################
##########################################################

#general plotting settings
options(scipen=999)
dev.off()
par(mfrow = c(1, 3))

#plot prior sens. for SS algorithm
algorithm = "SS"
prior_vec = c(0.01,0.2,0.5)
nr_priors = length(prior_vec)
for (prior in prior_vec){
  filename = paste0("prior_",prior,"_",algorithm,".Rdata")
  load(file=filename)
  plinks = output$plinks
  plinks_vec = plinks[upper.tri(plinks)]
  boxplot(plinks_vec,ylim=c(0,1),main=paste0(100*prior,"% prior sparsity"),boxwex=1.5)
}  
 
#plot prior sens. for BD-A algorithm
algorithm = "BDA"
prior_vec = c(0.01,0.2,0.5)
nr_priors = length(prior_vec)
for (prior in prior_vec){
  filename = paste0("prior_",prior,"_",algorithm,".Rdata")
  load(file=filename)
  plinks = output$plinks
  plinks_vec = plinks[upper.tri(plinks)]
  boxplot(plinks_vec,ylim=c(0,1),main=paste0(100*prior,"% prior sparsity"),boxwex=1.5)
} 


#######################################################################################
#############networks with edge incl prob and part. correlations########
#######################################################################################


#load final edge inclusion probabilities bd
par(mfrow = c(1, 2))
algorithm_vec = c("SS","BDA")
nr_links_show = 60

#make layout 
G_sum = matrix(0,100,100)
for (algorithm in algorithm_vec){
  filename = paste0("network_",algorithm,".Rdata") 
  load(file=filename)
  plinks = output$plinks
  
  #fill the lower triangle of the matrix
  p = 100
  count = 0
  for (i in 2:p){ #fill lower matrix
    for (j in 1:(i-1)){
      count = count + 1
      plinks[i,j] = plinks[j,i]
    } 
  }
  
  #determine x links with highest edge incl. prob.
  library(igraph)
  plinks_sort = sort(plinks,decreasing=TRUE)
  cutoff = 0.5*plinks_sort[nr_links_show*2]+0.5*plinks_sort[nr_links_show*2+1]
  index = which(plinks>cutoff)
  
  #create undirected graph with those x links
  G = matrix(0,100,100)
  G[index] = 1
  G_sum = G_sum + G #save sum graph
  
}
index_incl = which(G_sum > 0)
G_all = matrix(0,100,100)
G_all[index_incl] = 1
G_ig_all = graph_from_adjacency_matrix(adjmatrix=G_all,mode="undirected")
V(G_ig_all)$label = 1:100
layout = layout_nicely(G_ig_all,dim=2)

#plot graphs
for (algorithm in algorithm_vec){
  filename = paste0("network_",algorithm,".Rdata") 
  load(file=filename)
  plinks = output$plinks
  
  #fill the lower triangle of the matrix
  p = 100
  count = 0
  for (i in 2:p){ #fill lower matrix
    for (j in 1:(i-1)){
      count = count + 1
      plinks[i,j] = plinks[j,i]
    } 
  }
  
  #determine x links with highest edge incl. prob.
  library(igraph)
  plinks_sort = sort(plinks,decreasing=TRUE)
  cutoff = 0.5*plinks_sort[nr_links_show*2]+0.5*plinks_sort[nr_links_show*2+1]
  index = which(plinks>cutoff)
  
  #create undirected graph with those x links
  G = matrix(0,100,100)
  G[index] = 1
  G_ig = graph_from_adjacency_matrix(adjmatrix=G,mode="undirected")
  
  #create precision matrix K
  final_iter = dim(output$K_iter_matrix)[2]
  K_vec = output$K_iter_matrix[,final_iter]
  K = matrix(0,p,p)
  K[upper.tri(K,diag=TRUE)] = K_vec
  count = 0
  for (i in 2:p){ #fill lower matrix
    for (j in 1:(i-1)){
      count = count + 1
      K[i,j] = K[j,i]
    } 
  }
  
  #create partial correlation matrix C
  K_diag_vec = 1/sqrt(diag(K))
  K_diag = diag(K_diag_vec)
  C = -K_diag%*%K%*%K_diag
  
  #obtain vector of part. corr. of all edges that have incl. prob > cutoff
  G = matrix(0,p,p)
  G[index] = 1
  new_matrix = G*C
  part_cor_vec = new_matrix[lower.tri(new_matrix)]
  part_cor_vec = part_cor_vec[which(part_cor_vec!=0)]
  
  #make edges red that have neg part. correlation, and blue that have pos. partial correlation
  edge.color = rep("blue",length(part_cor_vec))
  edge.color[which(part_cor_vec<0)] = "red"
  
  #make a weights vector (that will determine the width of the edge according to the size of the partial correlation)
  weights = abs(part_cor_vec)
  
  #save layout and labels and remove vertices that are unconnected
  V(G_ig)$label = 1:100
  
  #delete unconnected nodes
  isolated = which(degree(G_ig)==0)
  layout_isolated = layout[-isolated,]
  G_ig_2 = delete.vertices(G_ig, isolated)
  
  #plot network of edge inclusion probabilities, edge width is edge inclusion prob, edge color is sign of partial correlation
  if (algorithm=="SS"){main = "SS-O algorithm"}
  if (algorithm=="BDA"){main = "BD-A algorithm"}
  plot.igraph(G_ig_2,mode="undirected",layout=layout_isolated,diag=FALSE,edge.width=10*weights,vertex.size=0.1,edge.color = edge.color,vertex.label.color ="black",vertex.label.dist=0.5,main=main)
}

#################################################################################################################################
###########find similarities between our network and networks in literature ####################################
#################################################################################################################################

#####find similarities with the K-Horseshoe paper ####
#load file 
algorithm = "BDA"
filename = paste0("network_",algorithm,".Rdata") 
load(file=filename)
plinks = output$plinks


#fill lower part of matrix
count = 0
for (i in 2:p){ #fill lower matrix
  for (j in 1:(i-1)){
    count = count + 1
    plinks[i,j] = plinks[j,i]
  } 
}
connectivity = rowSums(plinks)
sort(connectivity,decreasing=TRUE) 

#The below table gives the ranking of the 5 most connected genes in the K-Horseshoe paper
#GI_27754767-I, our rank: 1
#GI_3422299-S, our rank: 8
#GI_22027487-S, our rank: 30 something
#GI_29894333-A, our rank: 15
#GI_27477086-S, our rank: 50 something

gene_names = colnames(plinks)
gene_names[c(4,8,18,16,15,38,10)]
gene_names[c(94,41,60,39)]








  
