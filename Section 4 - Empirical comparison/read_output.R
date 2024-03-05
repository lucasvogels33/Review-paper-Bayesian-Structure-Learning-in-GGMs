#This file reads all the ".Rdata" files created by the file
#"run_and_save.R" and creates all the tables and plots of Section 4.

#this file has three sections
#1. Code to create convergence plots (Figure 1)
#2. Code to create AUC, Pr+, Pr- tables (Table 6,7 & 8)
#3. Code to create computation cost table (Table 9)

#before running this file, make sure to set the working directory to 
#the folder where all the .Rdata files are stored


#Section 1
##################################################
################convergence plots################
##################################################

p = 1000
n = 1050
graph = "random"
density = 0.001
algorithm_vec = c("Glasso","SSO","BAGUS" ,"GMPLBD" ,"Horseshoe")
yaxis = "p_min" #possible values are AUC, p_plus, or p_min
col_vec = c(1,4,5,6,7)
lty_vec = c(1,2,3,4,5,6)

title = ""
if (yaxis=="AUC"){ylim = c(0.5,1)}
if (yaxis=="p_plus"){ylim = c(0,1)}
if (yaxis=="p_min"){ylim = c(0,1)}
if (yaxis=="AUC"){ylab = "AUC"}
if (yaxis=="p_plus"){ylab = "Pr+"}
if (yaxis=="p_min"){ylab = "Pr-"}
plot(NA,xlim=c(0,15000),ylim=ylim,xlab="Time (in seconds)",ylab=ylab,main = title)

count = 0
count_lty = 0
for (method in algorithm_vec ){
  
  ##################
  ####load data#####
  ##################
  #initialize metrics
  AUC = 0
  p_plus = 0
  p_min = 0
  time = 0
  AUC_matrix = c()
  p_plus_matrix = c()
  p_min_matrix = c()
  time_matrix = c()
  
  rep_list = 1:16
  if ((p==100)*(n==40)*(graph=="random")*(density==0.01)){
    if(method=="RJWWA"){rep_list = c(1,2,7,9,11,12,13,14,15,16)} 
    if(method=="Horseshoe"){rep_list = c(2,3,4,6,10,11,13,14)} 
  }
  if ((p==100)*(n==40)*(graph=="random")*(density==0.1)){
    if(method=="RJWWA"){rep_list = c(1:8,10:16)} 
    if(method=="Horseshoe"){rep_list = c(1:8,10,11,13:16)} 
  }
  if ((p==100)*(n==700)*(graph=="random")*(density==0.01)){
    if(method=="RJWWA"){rep_list = c(1:6,8:16)} 
    if(method=="Horseshoe"){rep_list = c(2,3,15)} 
  }
  if ((p==100)*(n==700)*(graph=="random")*(density==0.1)){
    if(method=="RJWWA"){rep_list = c(1,2,4:16)} 
    if(method=="Horseshoe"){rep_list = c(1:4,6:15)} 
  }
  if ((p==100)*(n==40)*(graph=="cluster")*(density==0.01)){
    if(method=="RJWWA"){rep_list = c(1,2,5,6,7)} 
    if(method=="Horseshoe"){rep_list = c(1,2,4,6,10,13,14)} 
  }
  if ((p==100)*(n==40)*(graph=="cluster")*(density==0.1)){
    if(method=="RJWWA"){rep_list = c(1:4,6:16)} 
    if(method=="Horseshoe"){rep_list = c(1,2,3,4,5,6,8,10,11,13,14,15)} 
  }
  if ((p==100)*(n==700)*(graph=="cluster")*(density==0.01)){
    if(method=="RJWWA"){rep_list = 1:16} 
    if(method=="Horseshoe"){rep_list = c(2,4,13,15)} 
  }
  if ((p==100)*(n==700)*(graph=="cluster")*(density==0.1)){
    if(method=="RJWWA"){rep_list = c(1:6,8:16)} 
    if(method=="Horseshoe"){rep_list = c(1,2,3,4,7,8,10,11,12,13,15)} 
  }
  if ((p==1000)*(n==60)*(graph=="random")*(density==0.001)){
    if(method=="Horseshoe"){rep_list = c(1:4,6:16)} 
  }
  if ((p==1000)*(n==1050)*(graph=="random")*(density==0.001)){
    if(method=="Horseshoe"){rep_list = c(1:6,8:16)} 
  }
  replications = length(rep_list)
  
  
  for (rep in rep_list){
    
    
    #load file
    filename = paste0("results_p",p,"_n",n,"_",graph,"_dens",density,"_",method,"_rep",rep,".Rdata")
    load(file=filename)
    
    #save AUC and time for BGGM 
    if (method == "BGGM"){
      AUC = AUC + output$AUC
      p_plus = p_plus + output$pplus
      p_min = p_min + output$pmin
      time = time + output$time
    }
    #save AUC and time for Glasso
    else if (method=="Glasso"){
      AUC = AUC + output$AUC
      p_plus = p_plus + output$pplus_ric
      p_min = p_min + output$pmin_ric
      time = time + output$time_tuning_prep + output$time_tuning_ric + output$time_running_ric
    }
    #save AUC and time for BAGUS
    else if (method=="BAGUS"){
      time_matrix = rbind(time_matrix,c(output$time_tuning,output$time_tuning + output$time_vec))
      AUC_matrix = rbind(AUC_matrix,c(0.5,output$AUC_vec))
      p_plus_matrix = rbind(p_plus_matrix,c(0,output$pplus_vec))
      p_min_matrix = rbind(p_min_matrix,c(0,output$pmin_vec))
    }
    #save AUC and time for other methods
    else {
      time_matrix = rbind(time_matrix,output$time_vec)
      AUC_matrix = rbind(AUC_matrix,output$AUC_vec)
      p_plus_matrix = rbind(p_plus_matrix,output$pplus_vec)
      p_min_matrix = rbind(p_min_matrix,output$pmin_vec)
    }
    
  }
  
  
  #######################
  ####plot averages #####
  #######################
  count = count +1
  col = col_vec[count]
  
  #plot average point (in the case of BGGM or Glasso)
  if (method == "BGGM"){
    AUC = AUC/replications
    p_plus = p_plus/replications
    p_min = p_min/replications
    time = time/replications
    if (yaxis == "AUC"){points(x=time,y=AUC,col=col,bg = col,pch=21)}
    if (yaxis == "p_plus"){points(x=time,y=p_plus,col=col,bg = col,pch=21)}
    if (yaxis == "p_min"){points(x=time,y=p_min,col=col,bg = col,pch=21)}
  }
  else if (method=="Glasso"){
    AUC = AUC/replications
    p_plus = p_plus/replications
    p_min = p_min/replications
    time = time/replications
    if (yaxis == "AUC"){points(x=time,y=AUC,col=col,bg = col,pch=22)  }
    if (yaxis == "p_plus"){points(x=time,y=p_plus,col=col,bg = col,pch=22)}
    if (yaxis == "p_min"){points(x=time,y=p_min,col=col,bg = col,pch=22)}
  }
  #plot average lines (in the case of the other methods)
  else {
    count_lty = count_lty + 1
    lty = lty_vec[count_lty]
    time_vec = c(0,colMeans(time_matrix)) #we add the point AUC=0.5 at time t=0)
    AUC_vec = c(0.5,colMeans(AUC_matrix)) #we add the point AUC=0.5 at time t=0)
    p_plus_vec = c(0,colMeans(p_plus_matrix)) #we add the point p+ = 0 at time t=0)
    p_min_vec = c(0,colMeans(p_min_matrix)) #we add the point p- = 0 at time t=0)
    if (yaxis == "AUC"){points(x=time_vec,y=AUC_vec,type="l",col=col,lw=2,lty=lty) }
    if (yaxis == "p_plus"){points(x=time_vec,y=p_plus_vec,type="l",col=col,lw=2,lty=lty) }
    if (yaxis == "p_min"){points(x=time_vec,y=p_min_vec,type="l",col=col,lw=2,lty=lty) }
  }
  
}
if (p<250){
  legend( "topright", inset=.03, algorithm_vec, col = 1:8,pt.bg= c(1,NA,NA,NA,NA,NA,NA,8), pch = c(22,NA,NA,NA,NA,NA,NA,21), lw = c(NA,2,2,2,2,2,2,NA),lty = c(NA,1,1,1,1,1,1,NA), cex = 1, box.lty = 0 )
}
if (p==1000){
  legend( "topright", inset=.03, algorithm_vec, col = col_vec,pt.bg= c(1,NA,NA,NA,NA), pch = c(22,NA,NA,NA,NA), lw = c(NA,2,2,2,2),lty = c(NA,1,2,3,4), cex = 1, box.lty = 0 )
}



#Section 2
##################################################
################ AUC, Pr+, Pr- tables ############
##################################################

final_values = function(p=10,n=20,graph="cluster",density=0.1,rep_list=1:16,method="RJWWA"){
  
  #initiate vectors
  AUC_vec = c()
  pplus_vec = c()
  pmin_vec = c()
  
  #load data into the vectors
  for (rep in rep_list){
    
    #load file
    filename = paste0("results_p",p,"_n",n,"_",graph,"_dens",density,"_",method,"_rep",rep,".Rdata")  
    load(file = filename)
    
    #save metrics for BGGM
    if (method == "BGGM"){
      AUC_vec = c(AUC_vec,output$AUC)
      pplus_vec = c(pplus_vec,output$pplus)
      pmin_vec = c(pmin_vec,output$pmin)
    }
    #save metrics for GLASSO (currently only ric results are loaded)
    else if (method=="Glasso"){
      AUC_vec = c(AUC_vec,output$AUC)
      pplus_vec = c(pplus_vec,output$pplus_ric)
      pmin_vec = c(pmin_vec,output$pmin_ric)
    }
    #save AUC and time for other algorithms
    else {
      AUC_vec = c(AUC_vec,tail(output$AUC_vec, n=1))
      pplus_vec = c(pplus_vec,tail(output$pplus_vec, n=1))
      pmin_vec = c(pmin_vec,tail(output$pmin_vec, n=1))  
    }
    
  }
  
  output_list = list()
  output_list$AUC_vec = AUC_vec
  output_list$pplus_vec = pplus_vec
  output_list$pmin_vec = pmin_vec
  
  return(output_list)
}

make_latex_table = function(graph_list=c("cluster"),p_list = c(10),method_list=c("mpl_bd"),
                            file_p_plus="pplus.txt",file_p_min="pmin.txt",
                            file_AUC = "AUC.txt",round_nr=2){
  
  method_string = paste(method_list,collapse=" & ")
  
  for (paste_filename in c(file_p_plus,file_p_min,file_AUC)){
    cat("p & graph & density & n & ",method_string," \\\\",file=paste_filename,"\n", append = TRUE )
    cat("\\hline",file=paste_filename,"\n", append = TRUE )
  }
  
  for (p in p_list){
    print(p)
    for (graph in graph_list){
      print(graph)
      if (p==10){density_list=c(0.1)}
      if (p==100){density_list=c(0.01,0.1)}
      if (p==1000){density_list=c(0.001)}
      for (density in density_list){
      
        print(density)
      
        if (p==10){n_list=c(20,350)}
        if (p==100){n_list=c(40,700)}
        if (p==1000){n_list=c(60,1050)}
        
        for (n in n_list){
          for (paste_filename in c(file_p_plus,file_p_min,file_AUC)){
            cat(p,"&",graph," & ",density*100,"\\% &",n," &",  file=paste_filename, append = TRUE )
          }
          print(n)
          
          for (method in method_list){
            print(method)
            
            #BGGM does not work for n<p. BGGM, RJWWA and BDA do not work for p=1000
            if ((n < p)*(method=="BGGM")|(p==1000)*(method=="BGGM")|(p==1000)*(method=="RJWWA")|(p==1000)*(method=="BDA")){
              text = "-&"
              cat(text,file = file_p_plus, append = TRUE )
              cat(text,file = file_p_min, append = TRUE )
              cat(text,file = file_AUC, append = TRUE )
            }
            else {
              
              #there are 16 replications. However, for p=100 and for RJWWA and Horseshoe some replications failed
              rep_list = 1:16
              if ((p==100)*(n==40)*(graph=="random")*(density==0.01)){
                if(method=="RJWWA"){rep_list = c(1,2,7,9,11,12,13,14,15,16)} 
                if(method=="Horseshoe"){rep_list = c(2,3,4,6,10,11,13,14)} 
              }
              if ((p==100)*(n==40)*(graph=="random")*(density==0.1)){
                if(method=="RJWWA"){rep_list = c(1:8,10:16)} 
                if(method=="Horseshoe"){rep_list = c(1:8,10,11,13:16)} 
              }
              if ((p==100)*(n==700)*(graph=="random")*(density==0.01)){
                if(method=="RJWWA"){rep_list = c(1:6,8:16)} 
                if(method=="Horseshoe"){rep_list = c(2,3,15)} 
              }
              if ((p==100)*(n==700)*(graph=="random")*(density==0.1)){
                if(method=="RJWWA"){rep_list = c(1,2,4:16)} #RJWWA one replication failed
                if(method=="Horseshoe"){rep_list = c(1:4,6:15)} #Horseshoe one replication failed
              }
              
              if ((p==100)*(n==40)*(graph=="cluster")*(density==0.01)){
                if(method=="RJWWA"){rep_list = c(1,2,5,6,7)} 
                if(method=="Horseshoe"){rep_list = c(1,2,4,6,10,13,14)} 
              }
              if ((p==100)*(n==40)*(graph=="cluster")*(density==0.1)){
                if(method=="RJWWA"){rep_list = c(1:4,6:16)} 
                if(method=="Horseshoe"){rep_list = c(1,2,3,4,5,6,8,10,11,13,14,15)} 
              }
              if ((p==100)*(n==700)*(graph=="cluster")*(density==0.01)){
                if(method=="RJWWA"){rep_list = 1:16} 
                if(method=="Horseshoe"){rep_list = c(2,4,13,15)} 
              }
              if ((p==100)*(n==700)*(graph=="cluster")*(density==0.1)){
                if(method=="RJWWA"){rep_list = c(1:6,8:16)} 
                if(method=="Horseshoe"){rep_list = c(1,2,3,4,7,8,10,11,12,13,15)} 
              }
              if ((p==1000)*(n==60)*(graph=="random")*(density==0.001)){
                if(method=="Horseshoe"){rep_list = c(1:4,6:16)} 
              }
              if ((p==1000)*(n==1050)*(graph=="random")*(density==0.001)){
                if(method=="Horseshoe"){rep_list = c(1:6,8:16)} 
              }
              
              output_list = final_values(p=p,n=n,graph=graph,density=density,rep_list=rep_list,method=method)
              
              #p_plus
              mean = (round(mean(output_list$pplus_vec),round_nr))
              sd = (round(sd(output_list$pplus_vec),round_nr))
              text = paste0(mean," (",sd,") & ")
              cat(text,file = file_p_plus, append = TRUE )
              
              #p_min
              mean = (round(mean(output_list$pmin_vec),round_nr))
              sd = (round(sd(output_list$pmin_vec),round_nr))
              text = paste0(mean," (",sd,") & ")
              cat(text,file = file_p_min, append = TRUE )
              
              #AUC
              mean = (round(mean(output_list$AUC_vec),round_nr))
              sd = (round(sd(output_list$AUC_vec),round_nr))
              text = paste0(mean," (",sd,") & ")
              cat(text,file = file_AUC, append = TRUE )
            }
          
          } #closes method
          
          for (paste_filename in c(file_p_plus,file_p_min,file_AUC)){
            cat("\\\\",file = paste_filename,"\n", append = TRUE )
          }
          
        }#closes n
      }  #closes p
    }#closes density
  }#closes graph
  
}

p_list = c(10,100,1000)
graph_list = c("random","cluster")
method_list = c("Glasso", "RJWWA", "BDA" ,"SSO" ,"BAGUS" ,"GMPLBD" ,"Horseshoe" ,"BGGM")
make_latex_table(graph_list=graph_list,p_list = p_list,method_list=method_list)

#Section 3
####################################################
############## computational cost table  ###########
####################################################

#supporting functions
convergence_values = function(graph="random",p = 10,method="GMPLBD",n=100,density=0.1,rep_list=rep_list){
  
  cost_vec = c()
  for (rep in rep_list){
    
    #load AUC_vec and time_vec
    filename = paste0("results_p",p,"_n",n,"_",graph,"_dens",density,"_",method,"_rep",rep,".Rdata")
    load(filename)
    
    
    #save metrics for BGGM
    if (method == "Glasso"){
      conv_time = output$time_tuning_prep + output$time_tuning_ric + output$time_running_ric
      cost_vec = c(cost_vec,conv_time)
    }
    else if (method == "BGGM"){
      conv_time = output$time
      cost_vec = c(cost_vec,conv_time)
    }
    else {
      AUC_vec = output$AUC_vec
      time_vec = output$time_vec
      
      #calculate computational cost
      AUC_final = tail(AUC_vec, n=1)
      diff_AUC_vec = AUC_final - AUC_vec
      epsilon = 0.01
      conv_index = min(which(diff_AUC_vec < epsilon))
      conv_time = time_vec[conv_index]
      
      # add the tuning time for the BAGUS method
      if (method == "BAGUS"){
        conv_time = conv_time + output$time_tuning
      }
      
      #save computaional cost
      cost_vec = c(cost_vec,conv_time)
      
    }
    
  }
  
  return(cost_vec)
}

make_latex_table_comp_cost = function(graph_list=c("cluster"),p_list = c(10),method_list=c("mpl_bd"),
                            paste_filename = "cost.txt",round_nr=2){
  
  method_string = paste(method_list,collapse=" & ")
  
  
  cat("p & graph & density & n & ",method_string," \\\\",file=paste_filename,"\n", append = TRUE )
  cat("\\hline",file=paste_filename,"\n", append = TRUE )
  
  for (p in p_list){
    print(p)
    for (graph in graph_list){
      print(graph)
      if (p==10){density_list=c(0.1)}
      if (p==100){density_list=c(0.01,0.1)}
      if (p==1000){density_list=c(0.001)}
      for (density in density_list){
        
        print(density)
        
        if (p==10){n_list=c(20,350)}
        if (p==100){n_list=c(40,700)}
        if (p==1000){n_list=c(60,1050)}
        
        for (n in n_list){
          
          cat(p,"&",graph," & ",density*100,"\\% &",n," &",  file=paste_filename, append = TRUE )
          print(n)
          
          for (method in method_list){
            print(method)
            
            #BGGM does not work for n<p. BGGM, RJWWA and BDA do not work for p=1000
            if ((n < p)*(method=="BGGM")|(p==1000)*(method=="BGGM")|(p==1000)*(method=="RJWWA")|(p==1000)*(method=="BDA")){
              text = "-&"
              cat(text,file = paste_filename, append = TRUE )
            }
            else {
              
              #there are 16 replications. However, for p=100 and for RJWWA and Horseshoe some replications failed
              rep_list = 1:16
              if ((p==100)*(n==40)*(graph=="random")*(density==0.01)){
                if(method=="RJWWA"){rep_list = c(1,2,7,9,11,12,13,14,15,16)} 
                if(method=="Horseshoe"){rep_list = c(2,3,4,6,10,11,13,14)} 
              }
              if ((p==100)*(n==40)*(graph=="random")*(density==0.1)){
                if(method=="RJWWA"){rep_list = c(1:8,10:16)} 
                if(method=="Horseshoe"){rep_list = c(1:8,10,11,13:16)} 
              }
              if ((p==100)*(n==700)*(graph=="random")*(density==0.01)){
                if(method=="RJWWA"){rep_list = c(1:6,8:16)} 
                if(method=="Horseshoe"){rep_list = c(2,3,15)} 
              }
              if ((p==100)*(n==700)*(graph=="random")*(density==0.1)){
                if(method=="RJWWA"){rep_list = c(1,2,4:16)} #RJWWA one replication failed
                if(method=="Horseshoe"){rep_list = c(1:4,6:15)} #Horseshoe one replication failed
              }
              
              if ((p==100)*(n==40)*(graph=="cluster")*(density==0.01)){
                if(method=="RJWWA"){rep_list = c(1,2,5,6,7)} 
                if(method=="Horseshoe"){rep_list = c(1,2,4,6,10,13,14)} 
              }
              if ((p==100)*(n==40)*(graph=="cluster")*(density==0.1)){
                if(method=="RJWWA"){rep_list = c(1:4,6:16)} 
                if(method=="Horseshoe"){rep_list = c(1,2,3,4,5,6,8,10,11,13,14,15)} 
              }
              if ((p==100)*(n==700)*(graph=="cluster")*(density==0.01)){
                if(method=="RJWWA"){rep_list = 1:16} 
                if(method=="Horseshoe"){rep_list = c(2,4,13,15)} 
              }
              if ((p==100)*(n==700)*(graph=="cluster")*(density==0.1)){
                if(method=="RJWWA"){rep_list = c(1:6,8:16)} 
                if(method=="Horseshoe"){rep_list = c(1,2,3,4,7,8,10,11,12,13,15)} 
              }
              if ((p==1000)*(n==60)*(graph=="random")*(density==0.001)){
                if(method=="Horseshoe"){rep_list = c(1:4,6:16)} 
              }
              if ((p==1000)*(n==1050)*(graph=="random")*(density==0.001)){
                if(method=="Horseshoe"){rep_list = c(1:6,8:16)} 
              }
              
              
              cost_vec = convergence_values(graph=graph,p = p,method=method,n=n,density=density,rep_list=rep_list)
              
              #calculate mean and sd and paste in file
              mean = round(mean(cost_vec),round_nr)
              sd = round(sd(cost_vec),round_nr)
              text = paste0(mean," (",sd,") & ")
              cat(text,file = paste_filename, append = TRUE )
              
            }
            
          } #closes method
          
          
          cat("\\\\",file = paste_filename,"\n", append = TRUE )
          
          
        }#closes n
      }  #closes p
    }#closes density
  }#closes graph
  
}

#creating the table
graph_list = c("random","cluster")
p_list = c(10,100,1000)
method_list = c("Glasso", "RJWWA", "BDA" ,"SSO" ,"BAGUS" ,"GMPLBD" ,"Horseshoe" ,"BGGM")
make_latex_table_comp_cost(graph_list=graph_list,p_list = p_list,method_list=method_list,round_nr =0)


