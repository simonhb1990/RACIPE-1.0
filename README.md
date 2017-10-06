# RACIPE-1.0
/***********************************************************************
Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
 
Paper to cite:
Huang, Bin, Mingyang Lu, Dongya Jia, Eshel Ben-Jacob, Herbert Levine, and Jose N. Onuchic. "Interrogating the topological robustness of gene regulatory circuits by randomization." PLoS computational biology 13, no. 3 (2017): e1005456.

This code has used PCG Random Number Generation script by Melissa O'Neill
************************************************************************/


1. Use "make" to compile and use "make clean" to clean all compiled files. The executebale file is named "RACIPE" as default.

2. Run the package:
  1> Run with .topo file:
    $ ./RACIPE *.topo [options]
    RACIPE generates two files according to the topology information (in .topo file):
    - Configure information file (.cfg file) storing all setting for simulation (can be used in the second way to run the package);
    - Parameter ranges file (.prs file) storing ranges of the parameters for randomiztion.

    If the "-flag" option is set to 0 (default), RACIPE will calculate the results for all RACIPE models, which means "_parameter.dat", "_solution_.dat", and "_T_test.dat" will be generated.
    If the "-flag" option is set to 1, only .cfg and .prs will be generated.

    Example:
    $make
    $./RACIPE TS.topo [options]

  2> Run with .cfg file
    $ ./RACIPE *.cfg [options]
    .prs file will be automatically regenerated according to the setting in .cfg file no matter old .prs file exists or not.
    
    To change the setting of the simulation, you can either change the .cfg file directly or use options. However, the options set in command line will always override the setting in .cfg file. New configuration file will be generated as "_tmp.cfg" in the same fold to store the simulation setting for the current run.

    Similarily,
    If the "-flag" option is set to 0 (default), RACIPE will calculate the results for all RACIPE models, which means "_parameter.dat", "_solution_.dat", and "_T_test.dat" will be generated.
    If the "-flag" option is set to 1, updated configuration file will be generated as "_tmp.cfg" in the same fold.

3. Options
Use "./RACIPE -h" to find all available options.

-h             : Show all available options.
-maxtime       : Maximum time f-maor the simulation (Default 23.5 h).
-flag          : Only produce .cfg file or not (Default 0, not only produce .cfg file).

-KDID          : Gene or link (See their ID in .cfg file) to be knocked down. 
-OEID          : Gene (See their ID in .cfg file) to be overexpressed. (follow by -OEFD)
-OEFD          : Fold change to overexpress a gene (-OEID must be set ahead, the value need to be bigger than 1). (Default 1) if the corresponding OEFD is not set up, it will be set to 1.
-DEID          : Gene (See their ID in .cfg file) to be downexpressed. (follow by -DEFD)
-DEFD          : Fold change to downexpress a gene (-DEID must be set ahead, the value need to be bigger than 1). (Default 1) if the corresponding DEFD is not set up, it will be set to 1.
               Attention: The .prs file will be same as the one without Knockdown/downexpression/overexpression treatments. Multiple genes or links can be treated at the same time by putting multiple -KDID, -OEID, -DEID in the command line or modifing the .cfg file. 

-dist          : Distribution used for randomization:
                 1 ---> Uniform Distribution (Default);
                 2 ---> Guassian Distribution;
                 3 ---> Exponential Distribution.
-SF            : Scale the distribution ranges of parameters except hill coefficients, should be smaller than 1 (Default 1).
-num_findT     : Number of simulations used to estimate threshold (Default 10000).
-num_paras     : Number of RACIPE models to generate (Default 100).
-num_ode       : Number of Random initial values to solve ODEs (Default 100).
-num_stability : Maximum number of stable states to save for one RACIPE model (Default 10).
-thrd          : Cutoff for convergence of steady states for numerically solving ODEs (Default 1.0).
-Toggle_f_p    : Save parameters of each RACIPE model or not (Default 1 (yes)).
-stepsize      : Stepsize for solving ODE (Default 0.1).
-maxiters      : Maximum of Iteration for solving ODE at each RIVs (Default 20).
-Toggle_T_test : Test threshold assumption or not (Default 1 (yes)).
-seed          : Set up random seed (Default 1). 
                 Attention: the actual seed used by the package is a function of the starting time and the seed set here. -seed is used for the case you run the package at the same time for the same circuit several times. 

-minP          : Minimum production  rate (Default 1.0).
-maxP          : Maximum production  rate (Default 100.0).
-minK          : Minimum degradation rate (Default 0.1).
-maxK          : Maximum degradation rate (Default 1.0).
-minN          : Minimum coefficient (Default 1.0).
-maxN          : Maximum coefficient (Default 6.0).
-minF          : Minimum fold change (Default 1.0).
-maxF          : Maximum fold change (Default 100.0).

4. Input file
  1> Topology information file (.topo file)
     Format of topology file, such as TS.topo:
     Source  Target  Type
     A       B       1
     ...     ...     ...

     Types of regualtions: 1 --> Activation; 2 --> Inhibition;

  Or

  2> Configure information file (.cfg file), generated by the package. See details below.

5. Output file
  1> Standard output on screen, which can be redirected to the other file.
     It contains the Topology information, the result of testing the half-functional rule for each gene, the information of system stability, and running time.

  2> Configure information file (.cfg file).
     It generated by the package, contains the settings for the simulation, the ID for each gene (two columns: ID, Gene name) and regulatory links (four columns: ID, Gene ID, Gene ID, Type of regulation (same as topology information file)) 

     Attention:
     1. When directly modifing the .cfg file for diffferent distribution, change the number after 'Distribution', and replace the name of the distribution after the number (optional).
     2. When directly modifing the .cfg file for knockdown genes or links, change the 'number_of_KDs' to be the total number of genes and links to be knocked down, and put the Genes' ID and and Links' ID after 'KD_ID' separeted by space or tab.
     3. When directly modifing the .cfg file for over/down expression of a gene, change the 'number_of_OEs' and/or 'number_of_DEs' to the genes to be treated, and put the gene IDs after 'OE_ID' and/or 'DE_ID', also need to set up the fold changes after 'OE_Fold_Change' and/or 'DE_Fold_Change'.

  3> Parameter ranges file (.prs file). 
     It generated by the package, contains ranges of the parameters for randomiztion. It can be diretly modified when you run with .cfg file.

  4> Configure information file (_tmp.cfg), generated when the package runs with .cfg file, contains the updated settings for the simulation.
  
  5> _parameter.dat storing the parameters of each RACIPE model.
     Format of _parameter.dat:
     Model_index Number_of_stable_states Parameters_of_model

     number_of_stable_states            : Number of stable steady state of the RACIPE model.
     Parameters_of_model                : The meaning of each column is in the same order as the parameters in .prs file.

  6> _solution_.dat storing the gene expression for each stable steady state in log2 scale.
     The models with different number of stable states are stored in different files, e.g. monostable models are stored in _solution_1.dat, ang bistable models are stored in _solution_2.dat.
     Format of _solution_.dat:
     Model_index Number_of_stable_states Solutions

     Solutions               : The meaning of each column is consistent with the order of gene IDs in .cfg file.

  7> _T_test.dat storing the test of threshold assumption for each RACIPE model.
     Format of _T_test.dat:
     Model_index Over_threshold_A Below_threshold_A ... 

     Over_threshold_A  : Number of stable states for current RACIPE model whose expression of gene A is larger than its threshold parameter of A.
     Below_threshold_A : Number of stable states for current RACIPE model whose expression of gene A is smaller than its threshold parameter of A.

     For any model, the probability for gene A's expression to be larger than its threshold equals to the sum of Over_threshold_A across all the models devided by the sum of both Over_threshold_A and Below_threshold_A across all the models.
