This folder contains the code to reproduce the results of Section 4 ("Empirical comparison") of the article:
"Bayesian Structure Learning in Undirected Gaussian Graphical Models: Literature Review with Empirical Comparison"

This folder contains 5 types of files:
1. The main file: "run_and_save.R".
2. The main file: "read_output.R"
3. Supporting file: "metric_functions"
4. The eight supporting files: "..._functions"
5. The file "ggm_new5.cpp" 

Only the first two need to be ran. The others contain supporting code. The purpose of each file type is now explained.

1. Main file: "run_and_save.R".
    Running this file outputs two Rdata files: 
    i)  a file containing the true graphical structure and 
    ii) a file containing the estimated graphical structure 
        together with the AUC, Pr+ Pr-, and runtime at selected 
        iterations during the running of the algorithm
    Running "run_and_save.R" for different values of p, 
    n, graph, density algorithm and replication number, creates all 
    necessary ".Rdata" files. 

2. Main file: "read_output.R"
    This file reads all the ".Rdata" files created by the file
    "run_and_save.R" and creates all the tables and plots of Section 4.
    
3. Supporting file: "metric_functions"
    This file contains all supporting functions to calculate    
    the AUC, Pr+, Pr-, TPR and FPR                             

4. The eight supporting files: "..._functions"
    There is one file for each of the eight algorithms: 
    "1.Glasso_functions.R"
    "2.RJWWA_functions.R"
    ....
    "8.BGGM_functions.R" 
    All these files contain one main function and potentially one 
    or more supporting functions. The main function takes as input 
    the data of the true graph (a bdgraph.sim object) and outputs 
    the AUC, Pr- and Pr+. For iterative algorithms, these metrics 
    are provided at each iteration specified in the iter_vec_thin


5. The file "ggm_new5.cpp" 
    This file contains supporting code for the RJWWA algorithm 

