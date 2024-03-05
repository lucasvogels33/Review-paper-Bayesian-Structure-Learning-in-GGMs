##############################################################
# this is the main Rfile to produce the results in Section 4 #
# For a selected p, n, graph, density and algorithm, it      #
# outputs two Rdata files: i) a file containing the true     #
# graphical structure and ii) a file containing the AUC, Pr+ #
# Pr-, and runtime at selected iterations during the running #
# of the algorithm                                           #
##############################################################

#install packages for all algorithms (except RJ-WWA algorithm)
install.packages("BDgraph") #for BDA, GMPLD algorithms and for data creation
install.packages("ssgraph")    #for SS algorithm
install.packages("huge")       #for glasso algorithm
install.packages("BGGM")       #for BGGM algorithm

#download the necessary libraries
library(BDgraph)    #for BDA, GMPLD algorithms and for data creation
library(ssgraph)    #for SS algorithm
library(huge)       #for glasso algorithm
library(BGGM)       #for BGGM algorithm

#select the instance
p = 10              #number of variables
n = 20              #number of observations
graph = "random"    #graph type, options are "random" or "cluster"
density = 0.1       #density of the created graph
rep = 1             #replication number, should be an integer between 1 and 16.

#select the algorithm 
algorithm = "BDA"   #options are "Glasso", "RJWWA", "BDA", "SSO", "BAGUS" "GMPLBD" "Horseshoe" "BGGM"

#install packages and libraries for RJ-WWA algorithm 
if (algorithm=="RJWWA"){
    pkgs <- c("BH", "Rcpp", "withr")
    install.packages(
        pkgs = pkgs[!pkgs %in% rownames(installed.packages())], dependencies = TRUE,
        repos = "https://cran.rstudio.com/"
    )
    library(BH)         #for RJWWA algorithm
    library(Rcpp)       #for RJWWA algorithm
    library(withr)      #for RJWWA algorithm
    withr::with_makevars(
        new = c(PKG_LIBS = "$(LAPACK_LIBS) $(BLAS_LIBS)"),
        code = Rcpp::sourceCpp(file = "ggm_new5.cpp")
    )
}

#--load all supporting functions
source("1.Glasso_functions.R")
source("2.RJWWA_functions.R")
source("3.BDA_functions.R")
source("4.SSO_functions.R")
source("5.BAGUS_functions.R")
source("6.GMPLBD_functions.R")
source("7.Horseshoe_functions.R")
source("8.BGGM_functions.R")
source("metric_functions.R")

#--set value of "prob" -- prob is an argument in the data creation function bdgraph.sim -- It needs to be set so that the density of the resulting graph is equal to "density"#
if (graph=="random"){prob = density}
if ((graph=="cluster")*(p==10)*(density==0.1)){prob=0.225}
if ((graph=="cluster")*(p==100)*(density==0.01)){prob=0.020204}
if ((graph=="cluster")*(p==100)*(density==0.1)){prob=0.202041}
if ((graph=="cluster")*(p==1000)*(density==0.001)){prob=0.0080565}
if ((graph=="cluster")*(p==1000)*(density==0.01)){prob=0.0805645}

#set remaining parameters
type = "Gaussian"
size = NULL
vis = FALSE

#Algorithms calculate the metrics AUC, p+, p- at set iterations. These iterations are specified here for every instance and algorithm
if ((p==10)*(algorithm=="BDA")){iter_vec_thin = c(10,25,50,75,100,250,500,750,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000,15000,20000,25000,30000) }
if ((p==10)*(algorithm=="GMPLBD")){iter_vec_thin = c(10,50,100,1000,2000,3000,4000,5000,7500,10000,15000,20000,25000,30000) }
if ((p==10)*(algorithm=="SSO")){iter_vec_thin =c(10,20,30,40,50,75,100,250,500,750,1000,1500,2000,2500,3000)}
if ((p==10)*(algorithm=="BAGUS")){iter_vec_thin = c(1,2,3,4,5,10,15,20,30,50)}
if ((p==10)*(algorithm=="Horseshoe")){iter_vec_thin = c(10,50,100,200,300,400,500,750,1000,1500,2000)}
if ((p==10)*(algorithm=="RJWWA")){iter_vec_thin =c(10,50,100,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000,15000,20000,25000,30000)}

if ((p==100)*(algorithm=="BDA")){iter_vec_thin = c(10,50,100,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000,15000,20000,25000,30000,50000)}
if ((p==100)*(algorithm=="GMPLBD")){iter_vec_thin = c(10,50,100,1000,2000,3000,4000,5000,7500,10000,15000,20000,25000,30000)}
if ((p==100)*(algorithm=="SSO")){iter_vec_thin = c(10,20,30,40,50,75,100,250,500,750,1000,1500,2000,2500,3000)}
if ((p==100)*(algorithm=="BAGUS")){iter_vec_thin = c(1,2,3,4,5,10,15,20,25,30)}
if ((p==100)*(algorithm=="Horseshoe")){iter_vec_thin = c(5,10,20,30,40,50,75,100,200,300,400,500,750,1000,1500,2000)}
if ((p==100)*(algorithm=="RJWWA")){iter_vec_thin = c(10,50,100,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000,15000,20000,25000,30000,50000)}

if ((p==1000)*(density==0.001)*(algorithm=="BDA")){iter_vec_thin = c(10)} #we won't run this for p=1000 (not feasible)
if ((p==1000)*(density==0.001)*(algorithm=="GMPLBD")){iter_vec_thin = c(50,100,500,1000,2500,5000,10000,25000,50000,75000,100000,150000,200000)} 
if ((p==1000)*(density==0.001)*(algorithm=="SSO")){iter_vec_thin = c(10,20,30,40,50,75,100)}
if ((p==1000)*(density==0.001)*(algorithm=="BAGUS")){iter_vec_thin = c(1,2,3,4,5,10,15,20)}
if ((p==1000)*(density==0.001)*(algorithm=="Horseshoe")){iter_vec_thin = c(2,5,10,50,100)}
if ((p==1000)*(density==0.001)*(algorithm=="RJWWA")){iter_vec_thin = c(10)} #we won't run this for p=1000 (not feasible)

#Create data using bdgraph.sim() function and save the created true graph "G" and other parameters in a list called "data_output"
set.seed(rep)
data.sim = bdgraph.sim( p = p, n = n, graph = graph, prob = prob, size = size, type = type, vis = vis )
data.sim$density = density
data.sim$rep = rep
data_output = list(G = data.sim$G,p=p,n=n,graph = graph,density=density,rep=rep)

#Save the list "data_output" in an R.data file
filename = paste0("truth_p",p,"_n",n,"_",graph,"_dens",density,"_rep",rep,".Rdata")
save(data_output,file = filename)

#Solve the data and save the solutions in a list called "output"
if (algorithm=="Glasso"){output = glasso_solve(data=data.sim)}
if (algorithm=="RJWWA"){output = RJWWA_solve(data=data.sim,iter_vec_thin=iter_vec_thin)}
if (algorithm=="BDA"){output = BDA_solve(data=data.sim,iter_vec_thin=iter_vec_thin)}
if (algorithm=="SSO"){output = SSO_solve(data=data.sim,iter_vec_thin=iter_vec_thin)}
if (algorithm=="BAGUS"){output = BAGUS_solve(data=data.sim,iter_vec_thin=iter_vec_thin)}
if (algorithm=="GMPLBD"){output = GMPLBD_solve(data=data.sim,iter_vec_thin=iter_vec_thin)}
if (algorithm=="Horseshoe"){output = horseshoe_solve(data=data.sim,iter_vec_thin=iter_vec_thin)}
if (algorithm=="BGGM"){output = BGGM_solve(data=data.sim)}

#Save the solutions in an R.data file
filename = paste0("results_p",p,"_n",n,"_",graph,"_dens",density,"_",algorithm,"_rep",rep,".Rdata")
save(output,file = filename)








