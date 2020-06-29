/***********************************************************************
 Random Circuit Perturbation (RACIPE) menthod
 
 Copyright 2016 BIN HUANG <bh14@rice.edu>
 
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

# include <stdlib.h>
# include <stdint.h>
# include <stdio.h>
# include <math.h>
# include <time.h>
# include <string.h>
# include <inttypes.h>

# include "RACIPELIB.h"
# include "pcg_basic.h"
# include "rkf45.h"

/*********Shared Functions*********/
// Generate random value in (minV, maxV) following uniform distribution
double randu(double minV, double maxV)
{
  double u;
  
  do {
    u = ((double)pcg32_random()/UINT32_MAX);
  } while (u==0);

  return (minV + (maxV - minV)*u); 
}

// Generate random value following guassian distribtuion N(m, stadvalue)
double randg(double m, double stdvalue)
{
  double u1 = 0.0, u2 = 0.0, rsq = 0.0, fac = 0.0;
  static double z1 = 0.0, z2 = 0.0;   
  static int call = 0;

  if (call == 1){
    call = !call;
    return z2*stdvalue + m;
  }

  do {
    u1 = 2.0*randu(0, 1) - 1;
    u2 = 2.0*randu(0, 1) - 1;
    rsq = u1*u1 + u2*u2;
  } while (rsq >= 1.0 || rsq == 0.0);
  
  fac = sqrt(-2.0*log(rsq)/rsq);
  
  z1 = u1*fac;
  z2 = u2*fac;

  call = !call; 

  return z1*stdvalue + m;  
}

// Generate random value following non-negative guassian distribtuion N(m, stadvalue)
double randpg(double m, double stdvalue)
{
  double u = 0.0;
  do {
    u = randg(m, stdvalue);
  } while (u <= 0);

  return u;
}

// Generate random value following exponential distribution with mean equals to m
double randexp(double m)
{
  double z = 0.0;
  double exp_value = 0.0;

  z = randu(0, 1);

  exp_value = -m * log(z); // Compute exponential random variable using inversion method

  return(exp_value);
}

// Generate random value for fold change parameter which need to be greater than or equal to 1
double randfd(double m, double stdvalue, int dist)
{
  double u = 0.0;

  switch (dist) {
    case 2:  // Guassian distribution
      do {
        u = randg(m, stdvalue);
      } while (u < 1);
      break;
    case 3:  // Exponential distribution
      do {
        u = randexp(m);
      } while (u < 1);
      break;
    }

    return u;
}

// Generate discrete integer between minN and maxN
double randd(double minN, double maxN, int dist)
{
  double u = 0.0;
  double z = 0.0;
  int cnt  = 0;
  int i    = 0;

  do {
    if (dist == 1) {
      u = randu(minN-1, maxN+1);
    }
    else if (dist == 2){
      u = randpg((minN+maxN)/2.0, (maxN-minN)/2.0);
    }
    else if (dist == 3){
      u = randexp((minN+maxN)/2.0);
    }

    for (i = 0; i < (maxN-minN+1); i++){
      if (u >= minN-0.5+i && u < minN+0.5+i){
        z = minN + i;
        cnt = 1;
        break;
      }
    }
  } while (u == maxN+0.5 || cnt == 0);

  return z; 
}

// Shifted Hill Function
double Hillshift (double x, double x0, double nx, double lamda)
{
  double out;
  out = lamda + (1.0 - lamda) * (1.0/(1.0 + pow((x/x0), nx)));
  return out;
}

// Find median of an array
double median(double *x, int n)
{
  double  tmp;
  int     i = 0;
  int     j = 0;

  // Sorting
  for (i = 0; i < n-1; i++) {
    for (j = i+1; j < n; j++){
      if (x[j] < x[i]) {
        tmp  = x[i];
        x[i] = x[j];
        x[j] = tmp;
      }
    }
  }
  
  if (n%2==0) {
    return ((x[n/2] + x[n/2 - 1]) / 2.0);
  }
  else {
    return x[n/2];
  }
}

// Distance between two solutions
double sumdelta (double *y, double *ytmp, int NEQN)
{
    int      i = 0;
    double out = 0.0;
    
    for (i = 0; i < NEQN; i++){
        out = out + (y[i] - ytmp[i])*(y[i] - ytmp[i]);
    }

    return sqrt(out);
}

/*********Preprocess the topology file (.topo) or configure file (.cfg) for RACIPE method*********/
void check_inputfile (int argc, char **argv, struct topo *topoinfo, struct opts *simu_opts, struct rlt *tmprlt)
{
  char *token;
  char *fileextension;
  char *inputfile;

  // Check the input arguments
  if (argc == 1) {
    printf("### Missing: Missing input file!\n");
    exit(3); 
  }
  else if (argc == 2 && strcmp(argv[1], "-h") == 0){
    printf("Available options:\n");
    printf("-h             : Show all available options.\n");
    printf("-maxtime       : Maximum time for the simulation (Default 23.5 h).\n");
    printf("-solver        : The integrator method (1 --> Euler or 2 --> RK45) to solve the ODEs equations (Default 1).\n");
    printf("-flag          : run RACIPE to produce the .cfg file only or the whole simulation (Default 0, perform the whole simulation).\n");
    printf("-KDID          : Gene or link (See their ID in the .cfg file) to be knocked down.\n");
    printf("-OEID          : Gene (See their ID in the .cfg file) to be overexpressed. (follow by -OEFD).\n");
    printf("-OEFD          : Fold change to overexpress a gene (-OEID must be first set in the option, the value need to be bigger than 1). (Default 1) if the corresponding OEFD is not set, it will be set to 1.\n");
    printf("-DEID          : Gene (See their ID in the .cfg file) to be downexpressed. (follow by -DEFD).\n");
    printf("-DEFD          : Fold change to downexpress a gene (-DEID must be first set in the option, the value need to be bigger than 1). (Default 1) if the corresponding DEFD is not set, it will be set to 1.\n");
    printf("-dist          : Distribution used for randomization\n");
    printf("                 1 ---> Uniform Distribution (Default)\n");
    printf("                 2 ---> Guassian Distribution\n");
    printf("                 3 ---> Exponential Distribution\n");
    printf("-SF            : Scale the distribution ranges of all the parameters except for the hill coefficients, should be smaller than 1 (Default 1).\n");
    printf("-num_findT     : The number of simulations used to estimate the thresholds (Default 10000).\n");
    printf("-num_paras     : The number of RACIPE models to generate (Default 100).\n");
    printf("-num_ode       : The number of Random initial values to solve ODEs (Default 100).\n");
    printf("-num_stability : The maximum number of stable states to save for one RACIPE model (Default 10).\n"); 
    printf("-thrd          : Cutoff for convergence of steady states for numerically solving ODEs (Default 1.0).\n");
    printf("-Toggle_f_p    : Save parameters of each RACIPE model or not (Default 1 (yes)).\n");
    printf("-stepsize      : Stepsize for solving ODEs (Default 0.1).\n");
    printf("-maxiters      : The maximum number of iterations for solving ODEs at each random initial condition times 1000 (Default 20).\n");
    printf("-Toggle_T_test : Test the threshold assumption or not (Default 1 (yes)).\n");
    printf("-SBML_model    : Output a model in the SBML format. The parameter will be the ID of the model (start from 1) to save (Default 0 (no SBML output)).\n");
    printf("-seed          : random seed (Default 1).\n");
    printf("-minP          : Minimum production rate (Default 1.0)\n");
    printf("-maxP          : Maximum production rate (Default 100.0)\n");
    printf("-minK          : Minimum degradation rate (Default 0.1)\n");
    printf("-maxK          : Maximum degradation rate (Default 1.0)\n");
    printf("-minN          : Minimum Hill coefficient (Default 1.0)\n");
    printf("-maxN          : Maximum Hill coefficient (Default 6.0)\n");
    printf("-minF          : Minimum fold change (Default 1.0)\n");
    printf("-maxF          : Maximum fold change (Default 100.0)\n");
    exit(6);
  }

  inputfile = strdup(argv[1]);
  token                   = strtok(inputfile, ".");
  topoinfo->modelname     = strdup(token);
  token                   = strtok(NULL, ".");
  fileextension           = strdup(token);  // Whether .topo or .cfg file is the input file

  if (strcmp(fileextension, "topo") == 0) {
    simu_opts->exts = 0; // flag for input file as topo file
    initial_simuopts (argc, argv, topoinfo, simu_opts);
    Model_generate   (argv[1], topoinfo, simu_opts, tmprlt);
  }
  else if ((strcmp(fileextension, "cfg") == 0)) {
    simu_opts->exts = 1; // flag for input file as cfg file
    read_cfg         (topoinfo, simu_opts, tmprlt);
    initial_simuopts (argc, argv, topoinfo, simu_opts);
    Model_generate   (argv[1], topoinfo, simu_opts, tmprlt);
  }
  else {
    printf("### Wrong: Input file is not recognized. It should be .topo or .cfg file!\n");
    exit(4);
  }
}

// Initialize the setting options for the simulations.
void initial_simuopts (int argc, char **argv, struct topo *topoinfo, struct opts *simu_opts)
{

  int i = 0;
  int j = 0;
  int KDmarker = 0;
  int OEmarker = 0;
  int DEmarker = 0;
  char *end;
  FILE *testexistence;
  char KDIDname[100] = "";
  char parasname[100] = "";

  if (simu_opts->exts == 1) {
    
    strcpy(parasname, topoinfo->modelname);
    if (simu_opts->numKD != 0){
      for (i = 0; i < simu_opts->numKD; i++){
        sprintf(KDIDname, "%d", simu_opts->KDID[i]);
        strcat(parasname, "_");
        strcat(parasname, KDIDname);
      }
    }
    strcat(parasname, ".prs");
    testexistence = fopen(parasname, "r");
    if (testexistence != NULL) {
      fclose(testexistence);
      remove(parasname);
    }
  }

  // Default setting
  if (simu_opts->exts == 0){
    simu_opts->maxtime       = 23.5;
    simu_opts->solver         = 1;
    simu_opts->flag          = 0;
    simu_opts->numKD         = 0;
    simu_opts->KDID          = (int *)calloc(20, sizeof(int));
    simu_opts->numOE         = 0;
    simu_opts->OEID          = (int *)calloc(20, sizeof(int));
    simu_opts->OEFD          = (double *)calloc(20, sizeof(double));
    for (i = 0; i < 20; i++){
      simu_opts->OEFD[i] = 1.0;
    }
    simu_opts->numDE         = 0;
    simu_opts->DEID          = (int *)calloc(20, sizeof(int));
    simu_opts->DEFD          = (double *)calloc(20, sizeof(double));
    for (i = 0; i < 20; i++){
      simu_opts->DEFD[i] = 1.0;
    }
    simu_opts->distname      = strdup("Uniform");
    simu_opts->dist          = 1;
    simu_opts->SF            = 1.0;
    simu_opts->num_findT     = 10000;
    simu_opts->num_paras     = 100;
    simu_opts->num_ode       = 100; 
    simu_opts->num_stability = 10;
    simu_opts->thrd          = 1.0;
    simu_opts->Toggle_f_p    = 1;
    simu_opts->stepsize      = 0.1;
    simu_opts->maxiters      = 20;
    simu_opts->Toggle_T_test = 1;
    simu_opts->SBML_model    = 0;
    simu_opts->myseed        = 1;
    simu_opts->minP          = 1.0;
    simu_opts->maxP          = 100.0;
    simu_opts->minK          = 0.1;
    simu_opts->maxK          = 1.0;
    simu_opts->minN          = 1.0;
    simu_opts->maxN          = 6.0;
    simu_opts->minF          = 1.0;
    simu_opts->maxF          = 100.0;

    if (argc == 2){
      printf("### Warning: Uniform distribution is used.\n");
    }
  }
 
  // Set up the customized options
  if (argc >= 3){
    if ((argc - 2)%2 != 0){
        printf("### Wrong: No enough input arguments!\n");
        exit(8);
    }

    for (i = 2; i < argc; i = i + 2){

      if (strcmp(argv[i], "-dist") == 0){

        switch (atoi(argv[i+1])) {
          case 1 :
            printf("Uniform distribution is used.\n");
            simu_opts->distname = strdup("Uniform");
            simu_opts->dist = 1;
            break;
          case 2 :
            printf("Guassion distribution is used.\n");
            simu_opts->distname = strdup("Guassian");
            simu_opts->dist = 2;
            break;
          case 3 :
            printf("Exponential distribution is used.\n");
            simu_opts->distname = strdup("Exponential");
            simu_opts->dist = 3;
            break;
          default:
            printf("The distribution you selected is not recognized! Please select the follwing distribution:\n");
            printf("1 ---> Uniform Distribution\n");
            printf("2 ---> Guassian Distribution\n");
            printf("3 ---> Exponential Distribution\n");
            exit(3);
        }
      }
      else if (strcmp(argv[i], "-SF") == 0){
        simu_opts->SF = atof(argv[i+1]);
        printf("### Warning: Hill coefficients will not be scaled!\n");
      }
      else if (strcmp(argv[i], "-solver") == 0){
        simu_opts->solver = atoi(argv[i+1]);
        if (simu_opts->solver == 1) {
          printf("1st Euler method is used to solve ODEs.\n");
        }
        else if (simu_opts->solver == 2) {
          printf("RK45 method is used to solve ODEs.\n");
        }
        else {
          printf("### Wrong: no such options for solver!\n");
          exit(0);
        }
        
      }
      else if (strcmp(argv[i], "-num_findT") == 0){
        simu_opts->num_findT = atoi(argv[i+1]);
      }
      else if (strcmp(argv[i], "-num_paras") == 0){
        simu_opts->num_paras = atoi(argv[i+1]);
      }
      else if (strcmp(argv[i], "-num_ode") == 0){
        simu_opts->num_ode = atoi(argv[i+1]);
      }
      else if (strcmp(argv[i], "-num_stability") == 0){
        simu_opts->num_stability = atoi(argv[i+1]);
      }
      else if (strcmp(argv[i], "-thrd") == 0){
        simu_opts->thrd = atof(argv[i+1]);
      }
      else if (strcmp(argv[i], "-Toggle_f_p") == 0){
        simu_opts->Toggle_f_p = atoi(argv[i+1]);
      }
      else if (strcmp(argv[i], "-stepsize") == 0){
        simu_opts->stepsize = atof(argv[i+1]);
      }
      else if (strcmp(argv[i], "-maxiters") == 0){
        simu_opts->maxiters = atoi(argv[i+1]);
      }
      else if (strcmp(argv[i], "-Toggle_T_test") == 0){
        simu_opts->Toggle_T_test = atoi(argv[i+1]);
      }
      else if (strcmp(argv[i], "-SBML_model") == 0){
        simu_opts->SBML_model = atoi(argv[i+1]);
      }
      else if (strcmp(argv[i], "-KDID") == 0){
        if (simu_opts->exts == 1 && KDmarker == 0){
          KDmarker = 1;
          simu_opts->numKD = 0;
          simu_opts->KDID  = (int *)calloc(20, sizeof(int));
        }
        simu_opts->KDID[simu_opts->numKD] = atoi(argv[i+1]);
        simu_opts->numKD = simu_opts->numKD + 1;
      }
      else if (strcmp(argv[i], "-OEID") == 0){
        if (simu_opts->exts == 1 && OEmarker == 0){
          OEmarker = 1;
          simu_opts->numOE = 0;
          simu_opts->OEID  = (int *)calloc(20, sizeof(int));
          simu_opts->OEFD  = (double *)calloc(20, sizeof(double));
          for (j = 0; j < 20; j++){
            simu_opts->OEFD[j] = 1.0;
          }
        }
        simu_opts->OEID[simu_opts->numOE] = atoi(argv[i+1]);
        simu_opts->numOE = simu_opts->numOE + 1;
      }
      else if (strcmp(argv[i], "-OEFD") == 0){
        simu_opts->OEFD[simu_opts->numOE-1] = atof(argv[i+1]);
        // for (j = 0; j < 20; j++){
        //   printf("%f\t", simu_opts->OEFD[j]);
        // }
        // printf("\n");
      }
      else if (strcmp(argv[i], "-DEID") == 0){
        if (simu_opts->exts == 1 && DEmarker == 0){
          DEmarker = 1;
          simu_opts->numDE = 0;
          simu_opts->DEID  = (int *)calloc(20, sizeof(int));
          simu_opts->DEFD  = (double *)calloc(20, sizeof(double));
          for (j = 0; j < 20; j++){
            simu_opts->DEFD[j] = 1.0;
          }
        }
        simu_opts->DEID[simu_opts->numDE] = atoi(argv[i+1]);
        simu_opts->numDE = simu_opts->numDE + 1;
      }
      else if (strcmp(argv[i], "-DEFD") == 0){
        simu_opts->DEFD[simu_opts->numDE-1] = atof(argv[i+1]);
      }
      else if (strcmp(argv[i], "-flag") == 0){
        simu_opts->flag = atoi(argv[i+1]);
      }
      else if (strcmp(argv[i], "-maxtime") == 0){
        simu_opts->maxtime = atof(argv[i+1]);
      }
      else if (strcmp(argv[i], "-seed") == 0){
        simu_opts->myseed = strtoull(argv[i+1], &end, 10);
      }
      else if (strcmp(argv[i], "-minP") == 0){
        simu_opts->minP = atof(argv[i+1]);
      }
      else if (strcmp(argv[i], "-maxP") == 0){
        simu_opts->maxP = atof(argv[i+1]);
      }
      else if (strcmp(argv[i], "-minK") == 0){
        simu_opts->minK = atof(argv[i+1]);
      }
      else if (strcmp(argv[i], "-maxK") == 0){
        simu_opts->maxK = atof(argv[i+1]);
      }
      else if (strcmp(argv[i], "-minN") == 0){
        simu_opts->minN = atof(argv[i+1]);
      }
      else if (strcmp(argv[i], "-maxN") == 0){
        simu_opts->maxN = atof(argv[i+1]);
      }
      else if (strcmp(argv[i], "-minF") == 0){
        simu_opts->minF = atof(argv[i+1]);
      }
      else if (strcmp(argv[i], "-maxF") == 0){
        simu_opts->maxF = atof(argv[i+1]);
      }
      else{
        printf("### Wrong: Can not recognize the input arguments.\n");
        printf("Please use one of the follwing:\n");
        printf("-h             : Show all available options.\n");
        printf("-maxtime       : Maximum time for the simulation (Default 23.5 h).\n");
        printf("-solver        : The integrator method (1 --> Euler or 2 --> RK45) to solve the ODEs equations (Default 1).\n");
        printf("-flag          : run RACIPE to produce the .cfg file only or the whole simulation (Default 0, perform the whole simulation).\n");
        printf("-KDID          : Gene or link (See their ID in the .cfg file) to be knocked down.\n");
        printf("-OEID          : Gene (See their ID in the .cfg file) to be overexpressed. (follow by -OEFD).\n");
        printf("-OEFD          : Fold change to overexpress a gene (-OEID must be first set in the option, the value need to be bigger than 1). (Default 1) if the corresponding OEFD is not set, it will be set to 1.\n");
        printf("-DEID          : Gene (See their ID in the .cfg file) to be downexpressed. (follow by -DEFD).\n");
        printf("-DEFD          : Fold change to downexpress a gene (-DEID must be first set in the option, the value need to be bigger than 1). (Default 1) if the corresponding DEFD is not set, it will be set to 1.\n");
        printf("-dist          : Distribution used for randomization\n");
        printf("                 1 ---> Uniform Distribution (Default)\n");
        printf("                 2 ---> Guassian Distribution\n");
        printf("                 3 ---> Exponential Distribution\n");
        printf("-SF            : Scale the distribution ranges of all the parameters except for the hill coefficients, should be smaller than 1 (Default 1).\n");
        printf("-num_findT     : The number of simulations used to estimate the thresholds (Default 10000).\n");
        printf("-num_paras     : The number of RACIPE models to generate (Default 100).\n");
        printf("-num_ode       : The number of Random initial values to solve ODEs (Default 100).\n");
        printf("-num_stability : The maximum number of stable states to save for one RACIPE model (Default 10).\n"); 
        printf("-thrd          : Cutoff for convergence of steady states for numerically solving ODEs (Default 1.0).\n");
        printf("-Toggle_f_p    : Save parameters of each RACIPE model or not (Default 1 (yes)).\n");
        printf("-stepsize      : Stepsize for solving ODEs (Default 0.1).\n");
        printf("-maxiters      : The maximum number of iterations for solving ODEs at each random initial condition times 1000 (Default 20).\n");
        printf("-Toggle_T_test : Test the threshold assumption or not (Default 1 (yes)).\n");
        printf("-SBML_model    : Output a model in the SBML format. The parameter will be the ID of the model (start from 1) to save (Default 0 (no SBML output)).\n");
        printf("-seed          : random seed (Default 1).\n");
        printf("-minP          : Minimum production rate (Default 1.0)\n");
        printf("-maxP          : Maximum production rate (Default 100.0)\n");
        printf("-minK          : Minimum degradation rate (Default 0.1)\n");
        printf("-maxK          : Maximum degradation rate (Default 1.0)\n");
        printf("-minN          : Minimum Hill coefficient (Default 1.0)\n");
        printf("-maxN          : Maximum Hill coefficient (Default 6.0)\n");
        printf("-minF          : Minimum fold change (Default 1.0)\n");
        printf("-maxF          : Maximum fold change (Default 100.0)\n");
        exit(6);
      }
    }
  }
}

// Generate the models, parameter ranges for randomization and the configure file.
void Model_generate (char *inputfile, struct topo *topoinfo, struct opts *simu_opts, struct rlt *tmprlt)
{
  int i = 0;
  char KDIDname[100] = "";
  char OEIDname[100] = "";
  char DEIDname[100] = "";

  // Section 1 parameters .cfg
  FILE *fin;
  FILE *fout;
  char configname[100] = "";
  
  // Section 2 parameters .prs
  FILE   *fparas = NULL; 
  char   parasname[100] = "";

  // Section 1 --- Printout the number of genes involved and transform the topo file
  if (simu_opts->exts == 0){
    fin  = fopen(inputfile, "r");

    if (fin == NULL){
      printf("### Wrong: No topology files are found!\n");
      exit(2);
    }

    strcpy(configname, topoinfo->modelname);
    if (simu_opts->numKD != 0){
      strcat(configname, "_KD");
      for (i = 0; i < simu_opts->numKD; i++){
        sprintf(KDIDname, "%d", simu_opts->KDID[i]);
        strcat(configname, "_");
        strcat(configname, KDIDname);
      }
    }
    if (simu_opts->numOE != 0){
      strcat(configname, "_OE");
      for (i = 0; i < simu_opts->numOE; i++){
        sprintf(OEIDname, "%d", simu_opts->OEID[i]);
        strcat(configname, "_");
        strcat(configname, OEIDname);
      }
    }
    if (simu_opts->numDE != 0){
      strcat(configname, "_DE");
      for (i = 0; i < simu_opts->numDE; i++){
        sprintf(DEIDname, "%d", simu_opts->DEID[i]);
        strcat(configname, "_");
        strcat(configname, DEIDname);
      }
    }

    strcat(configname, ".cfg");
    fout = fopen(configname, "w");
  }
  else{
    strcpy(configname, topoinfo->modelname);
    strcat(configname, "_tmp.cfg");
    fout = fopen(configname, "w");
    printf("### Warning: Current configure information is stored in %s file\n", configname);
  }

  // Parameters for the simulations, stored in .cfg file or _tmp.cfg file
  fprintf(fout, "Distribution\t%d\t%s\n",                  simu_opts->dist, simu_opts->distname);
  fprintf(fout, "ScaleFactor\t%f\n",                       simu_opts->SF);
  fprintf(fout, "NumberOfSimulationToFindThreshold\t%d\n", simu_opts->num_findT);
  fprintf(fout, "NumberOfRACIPEModels\t%d\n",              simu_opts->num_paras);
  fprintf(fout, "NumberOfRIVs\t%d\n",                      simu_opts->num_ode);
  fprintf(fout, "NumberOfStatesToStore\t%d\n",             simu_opts->num_stability);
  fprintf(fout, "ThresholdForConvergence\t%f\n",           simu_opts->thrd);
  fprintf(fout, "ToggleOfSaveParameters\t%d\n",            simu_opts->Toggle_f_p);
  fprintf(fout, "Stepsize\t%f\n",                          simu_opts->stepsize);
  fprintf(fout, "MaximumOfIterations\t%d\n",               simu_opts->maxiters);
  fprintf(fout, "TestThreshold\t%d\n",                     simu_opts->Toggle_T_test);
  fprintf(fout, "number_of_KDs\t%d\n",                     simu_opts->numKD);
  fprintf(fout, "number_of_OEs\t%d\n",                     simu_opts->numOE);
  fprintf(fout, "number_of_DEs\t%d\n",                     simu_opts->numDE);

  if (simu_opts->numKD == 0){
    fprintf(fout, "KD_ID\t%d\n",                           0);
  }
  else{
    fprintf(fout, "KD_ID");
    for (i = 0; i < simu_opts->numKD; i++){
      fprintf(fout, "\t%d",                                simu_opts->KDID[i]);
    }
    fprintf(fout, "\n");
  }

  if (simu_opts->numOE == 0){
    fprintf(fout, "OE_ID\t%d\n",                           0);
    fprintf(fout, "OE_Fold_Change\t%f\n",                  0.0);
  }
  else{
    fprintf(fout, "OE_ID");
    for (i = 0; i < simu_opts->numOE; i++){
      fprintf(fout, "\t%d",                                simu_opts->OEID[i]);
    }
    fprintf(fout, "\n");

    fprintf(fout, "OE_Fold_Change");
    for (i = 0; i < simu_opts->numOE; i++){
      fprintf(fout, "\t%f",                                simu_opts->OEFD[i]);
    }
    fprintf(fout, "\n");
  }

  if (simu_opts->numDE == 0){
    fprintf(fout, "DE_ID\t%d\n",                           0);
    fprintf(fout, "DE_Fold_Change\t%f\n",                  0.0);
  }
  else{
    fprintf(fout, "DE_ID");
    for (i = 0; i < simu_opts->numDE; i++){
      fprintf(fout, "\t%d",                                simu_opts->DEID[i]);
    }
    fprintf(fout, "\n");
    fprintf(fout, "DE_Fold_Change");
    for (i = 0; i < simu_opts->numDE; i++){
      fprintf(fout, "\t%f",                                simu_opts->DEFD[i]);
    }
    fprintf(fout, "\n");
  }
  
  fprintf(fout, "MaximumRunningTime\t%f\n",                simu_opts->maxtime);
  fprintf(fout, "Seed\t%lld\n",                              simu_opts->myseed);
  fprintf(fout, "minP\t%f\n",                              simu_opts->minP);
  fprintf(fout, "maxP\t%f\n",                              simu_opts->maxP);
  fprintf(fout, "minK\t%f\n",                              simu_opts->minK);
  fprintf(fout, "maxK\t%f\n",                              simu_opts->maxK);
  fprintf(fout, "minN\t%f\n",                              simu_opts->minN);
  fprintf(fout, "maxN\t%f\n",                              simu_opts->maxN);
  fprintf(fout, "minF\t%f\n",                              simu_opts->minF);
  fprintf(fout, "maxF\t%f\n",                              simu_opts->maxF);

  if (simu_opts->exts == 0){
    check_topo(fin, fout, topoinfo); // read in the topo file
    fclose(fin);

    // initial struct rlt *tmprlt
    tmprlt->Nstb      = 0;
    tmprlt->numover   = (int *)    calloc(topoinfo->numG,                           sizeof(int));
    tmprlt->numdown   = (int *)    calloc(topoinfo->numG,                           sizeof(int));
    tmprlt->cnt_store = (int *)    calloc(simu_opts->num_stability,                 sizeof(int));
    tmprlt->y_store   = (double *) calloc(simu_opts->num_ode*topoinfo->numG,        sizeof(double));
    tmprlt->soln      = (double *) calloc(simu_opts->num_stability*topoinfo->numG,  sizeof(double));
    tmprlt->paras     = (double *) calloc(3*topoinfo->numR+2*topoinfo->numG,        sizeof(double));
  }

  fclose(fout); 
  fin    = NULL;
  fout   = NULL;
  
  // S2 --- Generate the randomization range
  strcpy(parasname, topoinfo->modelname);
  if (simu_opts->exts == 0){
    if (simu_opts->numKD != 0){
      strcat(parasname, "_KD");
      for (i = 0; i < simu_opts->numKD; i++){
        sprintf(KDIDname, "%d", simu_opts->KDID[i]);
        strcat(parasname, "_");
        strcat(parasname, KDIDname);
      }
    }
    if (simu_opts->numOE != 0){
      strcat(parasname, "_OE");
      for (i = 0; i < simu_opts->numOE; i++){
        sprintf(OEIDname, "%d", simu_opts->OEID[i]);
        strcat(parasname, "_");
        strcat(parasname, OEIDname);
      }
    }
    if (simu_opts->numDE != 0){
      strcat(parasname, "_DE");
      for (i = 0; i < simu_opts->numDE; i++){
        sprintf(DEIDname, "%d", simu_opts->DEID[i]);
        strcat(parasname, "_");
        strcat(parasname, DEIDname);
      }
    }
  }
  strcat(parasname, ".prs");

  if (simu_opts->exts == 0){
    fparas = fopen(parasname, "w+"); 
    generate_random_range(fparas, topoinfo, simu_opts);
  }
  else{
    fparas = fopen(parasname, "r");
    if (fparas == NULL){
      printf("### Warning: .prs file will be regenerated according to current .cfg file\n");
      fparas = fopen(parasname, "w+"); 
      generate_random_range(fparas, topoinfo, simu_opts);
    }
    else{
      printf("### Wrong: old .prs file is detected and needs to be deleted!\n");
      exit(11);
    }
  }
  fclose(fparas);
  fparas = NULL;

  if (simu_opts->flag == 1){
    printf(".cfg file is produced!\n");
    exit(10);
  }
}

// Read in the topology information
void check_topo(FILE *f_in, FILE *f_out, struct topo *topoinfo)
{
  int  tmpnumG  = 0;
  int  tmpnumR  = 0; 
  int  RT       = 0;  //regulation type
  int  count    = 0;
  int  i        = 0;
  int  j        = 0;
  int  cntP     = 0;
  int  matchnum = 0;
  char G1[100]  = "";
  char G2[100]  = "";

  topoinfo->Gname = (char**)malloc(sizeof(char *));

  if (f_in == NULL){
    printf("Topology file is missing!\n");
    exit(1);
  }

  rewind(f_in);
  
  fscanf(f_in, "%*[^\n]\n", NULL); //skip first line of topo file
  
  while (fscanf(f_in, "%s\t%s\t%d\n", G1, G2, &RT) == 3) {
    tmpnumR++;
                
    // Check the types of regualtions: 1 --> Activation; 2 --> Inhibition;                
    if (RT != 1 && RT != 2){
      printf("ERROR: Number %d regulation is not recognized\n", tmpnumR);
      exit(1);
    } 
    
    if (count == 0) { // second line of topo file
      
      if (strcmp(G1, G2) == 0) {
        count = count + 1;
        topoinfo->Gname = (char**)realloc(topoinfo->Gname, (count)*sizeof(char *));
        topoinfo->Gname[count-1] = (char*)malloc(sizeof(G1));
        strcpy(topoinfo->Gname[count-1], G1);
        tmpnumG = tmpnumG + 1;
      }
      else {
        count = count + 2;
        topoinfo->Gname = (char**)realloc(topoinfo->Gname, (count)*sizeof(char *));
        topoinfo->Gname[count-2] = (char*)malloc(sizeof(G1));
        topoinfo->Gname[count-1] = (char*)malloc(sizeof(G2));
        strcpy(topoinfo->Gname[count-2], G1);
        strcpy(topoinfo->Gname[count-1], G2);
        tmpnumG = tmpnumG + 2;
      }
    }
    else {
      matchnum = 0;
      // compare G1 first
      for (i = 0; i < count; i++) {
        if (strcmp(G1, topoinfo->Gname[i]) == 0) {
          matchnum = matchnum + 1;
        }
      }
      
      if (matchnum == 0) {
        count = count + 1;
        topoinfo->Gname = (char**)realloc(topoinfo->Gname, (count)*sizeof(char *));
        topoinfo->Gname[count-1] = (char*)malloc(sizeof(G1));
        strcpy(topoinfo->Gname[count-1], G1);
        tmpnumG = tmpnumG + 1;
      }
      
      if (strcmp(G1, G2) != 0) {
        matchnum = 0;
        for (i = 0; i < count; i++) {
          if (strcmp(G2, topoinfo->Gname[i]) == 0) {
            matchnum = matchnum + 1;
          }
        }
      
        if (matchnum == 0) {
          count = count + 1;
          topoinfo->Gname = (char**)realloc(topoinfo->Gname, (count)*sizeof(char *));
          topoinfo->Gname[count-1] = (char*)malloc(sizeof(G2));
          strcpy(topoinfo->Gname[count-1], G2);
          tmpnumG = tmpnumG + 1;
        }
      }
    }

  }
  
  // Screen printout
  printf("\n-------------------------------------------\n");
  printf("The topo file contains the following genes:\n");
  printf("Gene_ID -- Gene_Name\n");
  for (i = 0; i < count; i++) {
    printf("%d -- %s\n", i+1, topoinfo->Gname[i]);
  } 
  printf("-------------------------------------------\n");
  
  topoinfo->numR = tmpnumR;
  topoinfo->numG = tmpnumG;
  
  printf("Total number of genes = %d\n",       tmpnumG);
  printf("Total number of regulations = %d\n", tmpnumR);

  printf("-------------------------------------------\n");
  
  // transformation of the topo file  
  topoinfo->SourceG   = (int *)calloc(tmpnumR, sizeof(int));
  topoinfo->TargetG   = (int *)calloc(tmpnumR, sizeof(int));
  topoinfo->TypeR     = (int *)calloc(tmpnumR, sizeof(int));
  topoinfo->ParasPos  = (int *)calloc(tmpnumR, sizeof(int));
  count   = 0;
  
  rewind(f_in);
  fscanf(f_in, "%*[^\n]\n", NULL); //skip first line of topo file
  while (fscanf(f_in, "%s\t%s\t%d\n", G1, G2, &RT) == 3) {
    
    count++;
    
    for (i = 0; i < tmpnumG; i++) {
      if (strcmp(G1, topoinfo->Gname[i]) == 0) {
        topoinfo->SourceG[count-1] = i;
      }
      
      if (strcmp(G2, topoinfo->Gname[i]) == 0) {
        topoinfo->TargetG[count-1] = i;
      }
    }
    
    topoinfo->TypeR[count-1] = RT;
  }
  
  fprintf(f_out, "NumberOfRegulations\t%d\n",       topoinfo->numR);
  fprintf(f_out, "NumberOfGenes\t%d\n",             topoinfo->numG);
    
  for (i = 0; i < topoinfo->numG; i++) {
    fprintf(f_out, "%d\t%s\n", i+1, topoinfo->Gname[i]);
  } 

  for (i = 0; i < topoinfo->numR; i++){
    fprintf(f_out, "%d\t%d\t%d\t%d\n", topoinfo->numG+i+1, topoinfo->SourceG[i]+1, topoinfo->TargetG[i]+1, topoinfo->TypeR[i]);
  }

  for(i = 0; i < topoinfo->numG; i++){
    for (j = 0; j < topoinfo->numR; j++){
      if (topoinfo->TargetG[j] == i){ 
        topoinfo->ParasPos[j] = 3*cntP + 2*topoinfo->numG;  // Position of threshold parameters for each regulation.
        cntP = cntP + 1;
      }
    }
  }

  // for (i = 0; i < topoinfo->numR; i++){
  //   printf("%d\n", topoinfo->ParasPos[i]);
  // }
}

// Generate the parameter ranges for randomization and .prs file
void generate_random_range(FILE *f_paras, struct topo *topoinfo, struct opts *simu_opts)
{
  char   tmpparasname[100] = "";
  double tmpmin = 0.0;
  double tmpmax = 0.0;
  int    typeR  = 0;

  int    num    = simu_opts->num_findT;
  int    dist   = simu_opts->dist;
  double SF     = simu_opts->SF;

  /**** Default ranges for each class of parameters ****/
  double minP_d = simu_opts->minP;
  double maxP_d = simu_opts->maxP;
  double minK_d = simu_opts->minK;
  double maxK_d = simu_opts->maxK;
  double minN_d = simu_opts->minN;
  double maxN_d = simu_opts->maxN;
  double minF_d = simu_opts->minF;
  double maxF_d = simu_opts->maxF;
  /**** Default ranges for each class of parameters ****/

  double minP = 0.0;
  double maxP = 0.0;
  double minK = 0.0;
  double maxK = 0.0;
  double minN = 0.0;
  double maxN = 0.0;
  double minF = 0.0;
  double maxF = 0.0;

  double meanP = 0.0;
  double stdP  = 0.0;
  double meanK = 0.0;
  double stdK  = 0.0;
  double meanF = 0.0;
  double stdF  = 0.0;

  int i = 0;
  int j = 0;

  double *minT;
  double *maxT;
  minT = (double *)calloc(topoinfo->numG, sizeof(double));
  maxT = (double *)calloc(topoinfo->numG, sizeof(double));

  double *meanT;
  double *stdT;
  meanT = (double *)calloc(topoinfo->numG, sizeof(double));
  stdT  = (double *)calloc(topoinfo->numG, sizeof(double));

  double *amplifyfold; // The fold changes need to applied in order to make the minT bigger than 0.01
  double tmpfold;
  amplifyfold = (double *)calloc(topoinfo->numG, sizeof(double));

  switch (dist) {
    case 1 :
      minP = ((minP_d + maxP_d) - SF * (maxP_d - minP_d))/2.0;
      maxP = ((minP_d + maxP_d) + SF * (maxP_d - minP_d))/2.0;
      minK = ((minK_d + maxK_d) - SF * (maxK_d - minK_d))/2.0;
      maxK = ((minK_d + maxK_d) + SF * (maxK_d - minK_d))/2.0;
      minN = minN_d;
      maxN = maxN_d;
      minF = ((minF_d + maxF_d) - SF * (maxF_d - minF_d))/2.0;
      maxF = ((minF_d + maxF_d) + SF * (maxF_d - minF_d))/2.0;

      fprintf(f_paras, "Parameter\tMinimum_value\tMaximum_Value\tRegulation_type\n");

      // Format of parameter files: name minV maxV type_of_regulation
      // estimate the threshold
      printf("Amplification of the parameter ranges (production rates and thresholds)\n");
      for(i = 0; i < topoinfo->numG; i++){
        // estimate threshold and printout
        estimate_threshold(num, i, minP, maxP, minK, maxK, minN, maxN, minF, maxF, minT, maxT, topoinfo, dist, SF);
        // printf("%s = %f\n", topoinfo->Gname[i], (minT+maxT)/2.0);
        if (minT[i] < 0.01){
          tmpfold = 10.0;
          while ((minT[i]*tmpfold) < 0.01){
            tmpfold = tmpfold*10.0;
          }
          amplifyfold[i] = tmpfold;
        }
        else{
          amplifyfold[i] = 1;
        }
        printf("%s\t%f\n", topoinfo->Gname[i], amplifyfold[i]);
      }

      // production rate  
      for (i = 0; i < topoinfo->numG; i++){
        fprintf(f_paras, "Prod_of_%s\t%f\t%f\t%d\n", topoinfo->Gname[i], minP*amplifyfold[i], maxP*amplifyfold[i], 0);
      }
      
      // degradation rate
      for (i = 0; i < topoinfo->numG; i++){
        fprintf(f_paras, "Deg_of_%s\t%f\t%f\t%d\n", topoinfo->Gname[i], minK, maxK, 0);
      }

      for(i = 0; i < topoinfo->numG; i++){
        for (j = 0; j < topoinfo->numR; j++){
          if (topoinfo->TargetG[j] == i){ 
            // Threshold
            fprintf(f_paras,   "Trd_of_%sTo%s\t%f\t%f\t%d\n", topoinfo->Gname[topoinfo->SourceG[j]], topoinfo->Gname[topoinfo->TargetG[j]], minT[topoinfo->SourceG[j]]*amplifyfold[topoinfo->SourceG[j]], maxT[topoinfo->SourceG[j]]*amplifyfold[topoinfo->SourceG[j]], 0);
            // Number of binding sites
            fprintf(f_paras,   "Num_of_%sTo%s\t%f\t%f\t%d\n", topoinfo->Gname[topoinfo->SourceG[j]], topoinfo->Gname[topoinfo->TargetG[j]], minN, maxN, 0);
            // Fold change of a regulation
            if (topoinfo->TypeR[j] == 1) {  //Activation
              fprintf(f_paras, "Act_of_%sTo%s\t%f\t%f\t%d\n", topoinfo->Gname[topoinfo->SourceG[j]], topoinfo->Gname[topoinfo->TargetG[j]], minF, maxF, topoinfo->TypeR[j]);
            }
            else if (topoinfo->TypeR[j] == 2) { //Inhibition
              fprintf(f_paras, "Inh_of_%sTo%s\t%f\t%f\t%d\n", topoinfo->Gname[topoinfo->SourceG[j]], topoinfo->Gname[topoinfo->TargetG[j]], minF, maxF, topoinfo->TypeR[j]);
            }
          }
        }
      }

      free(minT);
      free(maxT);
      free(meanT);
      free(stdT);
      break;
    case 2 :
      // Format of parameter files: name meanV stdV
      meanP = (minP_d+maxP_d)/2.0;
      stdP  = (maxP_d-minP_d)*SF/2.0;
      meanK = (minK_d+maxK_d)/2.0;
      stdK  = (maxK_d-minK_d)*SF/2.0;
      meanF = (minF_d+maxF_d)/2.0;
      stdF  = (maxF_d-minF_d)*SF/2.0;

      fprintf(f_paras, "Parameter\tMean\tStandard_deviation\tRegulation_type\n");

      printf("Amplification of the parameter ranges (production rates and thresholds)\n");
      for(i = 0; i < topoinfo->numG; i++){
        // estimate threshold and printout
        estimate_threshold(num, i, meanP, stdP, meanK, stdK, minN_d, maxN_d, meanF, stdF, meanT, stdT, topoinfo, dist, SF);
        // printf("%s = %f\n", topoinfo->Gname[i], (minT+maxT)/2.0);
        if (minT[i] < 0.01){
          tmpfold = 10.0;
          while ((minT[i]*tmpfold) < 0.01){
            tmpfold = tmpfold*10.0;
          }
          amplifyfold[i] = tmpfold;
        }
        else{
          amplifyfold[i] = 1;
        }
        printf("%s\t%f\n", topoinfo->Gname[i], amplifyfold[i]);
      }

      // production rate  
      for (i = 0; i < topoinfo->numG; i++){
        fprintf(f_paras, "Prod_of_%s\t%f\t%f\t%d\n", topoinfo->Gname[i], meanP*amplifyfold[i], stdP*amplifyfold[i], 0);
      }
      
      // degradation rate
      for (i = 0; i < topoinfo->numG; i++){
        fprintf(f_paras, "Deg_of_%s\t%f\t%f\t%d\n", topoinfo->Gname[i], meanK, stdK, 0);
      }

      for(i = 0; i < topoinfo->numG; i++){
        for (j = 0; j < topoinfo->numR; j++){
          if (topoinfo->TargetG[j] == i){ 
            // Threshold
            fprintf(f_paras,   "Trd_of_%sTo%s\t%f\t%f\t%d\n", topoinfo->Gname[topoinfo->SourceG[j]], topoinfo->Gname[topoinfo->TargetG[j]], meanT[topoinfo->SourceG[j]]*amplifyfold[topoinfo->SourceG[j]], stdT[topoinfo->SourceG[j]]*amplifyfold[topoinfo->SourceG[j]], 0);
            // Number of binding sites
            fprintf(f_paras,   "Num_of_%sTo%s\t%f\t%f\t%d\n", topoinfo->Gname[topoinfo->SourceG[j]], topoinfo->Gname[topoinfo->TargetG[j]], minN_d, maxN_d, 0);
            // Fold change of a regulation
            if (topoinfo->TypeR[j] == 1) {  //Activation
              fprintf(f_paras, "Act_of_%sTo%s\t%f\t%f\t%d\n", topoinfo->Gname[topoinfo->SourceG[j]], topoinfo->Gname[topoinfo->TargetG[j]], meanF, stdF, topoinfo->TypeR[j]);
            }
            else if (topoinfo->TypeR[j] == 2) { //Inhibition
              fprintf(f_paras, "Inh_of_%sTo%s\t%f\t%f\t%d\n", topoinfo->Gname[topoinfo->SourceG[j]], topoinfo->Gname[topoinfo->TargetG[j]], meanF, stdF, topoinfo->TypeR[j]);
            }
          }
        }
      }

      free(minT);
      free(maxT);
      free(meanT);
      free(stdT);
      break;
    case 3 :
      // Format of parameter files: name meanV stdV
      meanP = (minP_d+maxP_d)*SF/2.0;
      meanK = (minK_d+maxK_d)*SF/2.0;
      meanF = (minF_d+maxF_d)*SF/2.0;

      fprintf(f_paras, "Parameter\tMean\tNo_sense\tRegulation_type\n");
      
      printf("Amplification of the parameter ranges (production rates and thresholds)\n");
      for(i = 0; i < topoinfo->numG; i++){
        // estimate threshold and printout
        estimate_threshold(num, i, meanP, stdP, meanK, stdK, minN_d, maxN_d, meanF, stdF, meanT, stdT, topoinfo, dist, SF);
        // printf("%s = %f\n", topoinfo->Gname[i], (minT+maxT)/2.0);
        if (minT[i] < 0.01){
          tmpfold = 10.0;
          while ((minT[i]*tmpfold) < 0.01){
            tmpfold = tmpfold*10.0;
          }
          amplifyfold[i] = tmpfold;
        }
        else{
          amplifyfold[i] = 1;
        }
        printf("%s\t%f\n", topoinfo->Gname[i], amplifyfold[i]);
      }

      // production rate  
      for (i = 0; i < topoinfo->numG; i++){
        fprintf(f_paras, "Prod_of_%s\t%f\t%f\t%d\n", topoinfo->Gname[i], meanP*amplifyfold[i], stdP, 0);
      }
      
      // degradation rate
      for (i = 0; i < topoinfo->numG; i++){
        fprintf(f_paras, "Deg_of_%s\t%f\t%f\t%d\n", topoinfo->Gname[i], meanK, stdK, 0);
      }

      for(i = 0; i < topoinfo->numG; i++){
        for (j = 0; j < topoinfo->numR; j++){
          if (topoinfo->TargetG[j] == i){ 
            // Threshold
            fprintf(f_paras,   "Trd_of_%sTo%s\t%f\t%f\t%d\n", topoinfo->Gname[topoinfo->SourceG[j]], topoinfo->Gname[topoinfo->TargetG[j]], meanT[topoinfo->SourceG[j]]*amplifyfold[topoinfo->SourceG[j]], stdT[topoinfo->SourceG[j]], 0);
            // Number of binding sites
            fprintf(f_paras,   "Num_of_%sTo%s\t%f\t%f\t%d\n", topoinfo->Gname[topoinfo->SourceG[j]], topoinfo->Gname[topoinfo->TargetG[j]], minN_d, maxN_d, 0);
            // Fold change of a regulation
            if (topoinfo->TypeR[j] == 1) {  //Activation
              fprintf(f_paras, "Act_of_%sTo%s\t%f\t%f\t%d\n", topoinfo->Gname[topoinfo->SourceG[j]], topoinfo->Gname[topoinfo->TargetG[j]], meanF, stdF, topoinfo->TypeR[j]);
            }
            else if (topoinfo->TypeR[j] == 2) { //Inhibition
              fprintf(f_paras, "Inh_of_%sTo%s\t%f\t%f\t%d\n", topoinfo->Gname[topoinfo->SourceG[j]], topoinfo->Gname[topoinfo->TargetG[j]], meanF, stdF, topoinfo->TypeR[j]);
            }
          }
        }
      }

      free(minT);
      free(maxT);
      free(meanT);
      free(stdT);
      break;
  }

  topoinfo->prsrandrange    = (double **)calloc(3,                                        sizeof(double *));
  topoinfo->prsrandrange[0] = (double *) calloc(3*topoinfo->numR+2*topoinfo->numG,        sizeof(double));
  topoinfo->prsrandrange[1] = (double *) calloc(3*topoinfo->numR+2*topoinfo->numG,        sizeof(double));
  topoinfo->prsrandrange[2] = (double *) calloc(3*topoinfo->numR+2*topoinfo->numG,        sizeof(double));

  rewind(f_paras);
  fscanf(f_paras, "%*[^\n]\n", NULL);
  i = 0;

  while (fscanf(f_paras, "%s\t%lf\t%lf\t%d", tmpparasname, &tmpmin, &tmpmax, &typeR) == 4){
    topoinfo->prsrandrange[0][i] = tmpmin;
    topoinfo->prsrandrange[1][i] = tmpmax;
    topoinfo->prsrandrange[2][i] = (double)typeR;
    // printf("%f\t%f\t%f\n", topoinfo->prsrandrange[0][i], topoinfo->prsrandrange[1][i], topoinfo->prsrandrange[2][i]);
    i++;
  }
}

// Estimate the threshold range for randomization, called in generate_random_range
void estimate_threshold(int num, int ID, double minP, double maxP, double minK, double maxK, double minN, double maxN, double minF, double maxF, double *minT, double *maxT, struct topo *topoinfo, int dist, double SF)
{
  int    i    = 0;
  int    j    = 0;
  int    numA = 0;
  int    numI = 0;

  double g      = 0.0;
  double k      = 0.0;
  double T      = 0.0;
  double n      = 0.0;
  double lambda = 0.0;
  double MA     = 0.0;
  double MB     = 0.0;
  double f1     = 0.0;
  double f2     = 0.0;

  double *A; // A is a standalone gene
  double *B;

  A  = (double *)calloc(num, sizeof(double));
  B  = (double *)calloc(num, sizeof(double));

  for (i = 0; i < topoinfo->numR; i++){
    if (topoinfo->TargetG[i] == ID){
      if (topoinfo->TypeR[i] == 1){ //Activation
        numA = numA + 1;        
      }
      else if (topoinfo->TypeR[i] == 2){ //Inhibition
        numI = numI + 1;
      }
    }
  }

  f1 = (2.0 - SF*1.96)/2.0;
  f2 = (2.0 + SF*1.96)/2.0;

  switch (dist){
    case 1:    // Uniform distribution
      for (i = 0; i < num; i++){
        g    = randu(minP, maxP);
        k    = randu(minK, maxK);
        A[i] = g/k;
      }   

      MA = median(A, num); 

      for (i = 0; i < num; i++){
        g    = randu(minP, maxP);
        k    = randu(minK, maxK);
        B[i] = g/k;
        
        if (numA != 0){
          for (j = 0; j < numA; j++){
            g      = randu(minP, maxP);
            k      = randu(minK, maxK);
            n      = randd(minN, maxN, dist);
            T      = randu(MA*f1, MA*f2);
            lambda = randu(minF, maxF);

            B[i] = B[i]*Hillshift(g/k, T, n, lambda)/lambda;
          }
        }
        
        if (numI != 0){
          for (j = 0; j < numI; j++){
            g      = randu(minP, maxP);
            k      = randu(minK, maxK);
            n      = randd(minN, maxN, dist);
            T      = randu(MA*f1, MA*f2);
            lambda = 1.0/randu(minF, maxF);

            B[i] = B[i]*Hillshift(g/k, T, n, lambda);
          }
        }
      }

      MB = median(B, num);
      minT[ID] = MB*f1;
      maxT[ID] = MB*f2;
      free(A);
      free(B);
      break;
    case 2:    // non-negative Guassian distribution
      for (i = 0; i < num; i++){
          g    = randpg(minP, maxP);
          k    = randpg(minK, maxK);
          A[i] = g/k;
        }   

        MA = median(A, num); 

        for (i = 0; i < num; i++){
          g    = randpg(minP, maxP);
          k    = randpg(minK, maxK);
          B[i] = g/k;
          
          if (numA != 0){
            for (j = 0; j < numA; j++){
              g      = randpg(minP, maxP);
              k      = randpg(minK, maxK);
              n      = randd(minN, maxN, dist);
              T      = randpg(MA, (MA*f2-MA*f1)/2.0);
              lambda = randfd(minF, maxF, dist);

              B[i] = B[i]*Hillshift(g/k, T, n, lambda)/lambda;
            }
          }
          
          if (numI != 0){
            for (j = 0; j < numI; j++){
              g      = randpg(minP, maxP);
              k      = randpg(minK, maxK);
              n      = randd(minN, maxN, dist);
              T      = randpg(MA, (MA*f2-MA*f1)/2.0);
              lambda = 1.0/randfd(minF, maxF, dist);

              B[i] = B[i]*Hillshift(g/k, T, n, lambda);
            }
          }
        }

        MB = median(B, num);
        minT[ID] = MB;  //mean
        maxT[ID] = (MB*f2-MB*f1)/2.0; //standard deviation
        free(A);
        free(B);
        break;
    case 3:    // Exponential distribution
      for (i = 0; i < num; i++){
        g    = randexp(minP);
        k    = randexp(minK);
        A[i] = g/k;
      }   

      MA = median(A, num); 

      for (i = 0; i < num; i++){
        g    = randexp(minP);
        k    = randexp(minK);
        B[i] = g/k;
        
        if (numA != 0){
          for (j = 0; j < numA; j++){
            g      = randexp(minP);
            k      = randexp(minK);
            n      = randd(minN, maxN, dist);
            T      = randexp(MA);
            lambda = randfd(minF, 0, dist);

            B[i] = B[i]*Hillshift(g/k, T, n, lambda)/lambda;
          }
        }
        
        if (numI != 0){
          for (j = 0; j < numI; j++){
            g      = randexp(minP);
            k      = randexp(minK);
            n      = randd(minN, maxN, dist);
            T      = randexp(MA);
            lambda = 1.0/randfd(minF, 0, dist);

            B[i] = B[i]*Hillshift(g/k, T, n, lambda);
          }
        }
      }

      MB = median(B, num);
      minT[ID] = MB;
      maxT[ID] = 0.0;
      free(A);
      free(B);
      break;
    }
}

// Read in the configure information
void read_cfg(struct topo *topoinfo, struct opts *simu_opts, struct rlt *tmprlt)
{
  int tmpID = 0;
  int cntP  = 0;
  int i     = 0;
  int j     = 0;

  int tmpSourceG = 0;
  int tmpTargetG = 0;

  char configname[100] = "";

  strcpy(configname, topoinfo->modelname);
  strcat(configname, ".cfg");

  FILE *fcfg;
  char tmpparasname[1000]  = "";
  char tmpparasname2[1000]  = "";

  fcfg = fopen(configname, "r");

  if (fcfg == NULL){
    printf("No configure file provided!\n");
    exit(2);
  }

  rewind(fcfg);

  // simulation settings
  fscanf(fcfg, "%s\t%d\t%s\n",       tmpparasname,  &simu_opts->dist, tmpparasname2);
  fscanf(fcfg, "%s\t%lf\n",          tmpparasname,  &simu_opts->SF);
  fscanf(fcfg, "%s\t%d\n",           tmpparasname,  &simu_opts->num_findT);
  fscanf(fcfg, "%s\t%d\n",           tmpparasname,  &simu_opts->num_paras);
  fscanf(fcfg, "%s\t%d\n",           tmpparasname,  &simu_opts->num_ode);
  fscanf(fcfg, "%s\t%d\n",           tmpparasname,  &simu_opts->num_stability);
  fscanf(fcfg, "%s\t%lf\n",          tmpparasname,  &simu_opts->thrd);
  fscanf(fcfg, "%s\t%d\n",           tmpparasname,  &simu_opts->Toggle_f_p);
  fscanf(fcfg, "%s\t%lf\n",          tmpparasname,  &simu_opts->stepsize);
  fscanf(fcfg, "%s\t%d\n",           tmpparasname,  &simu_opts->maxiters);
  fscanf(fcfg, "%s\t%d\n",           tmpparasname,  &simu_opts->Toggle_T_test);
  fscanf(fcfg, "%s\t%d\n",           tmpparasname,  &simu_opts->numKD);
  fscanf(fcfg, "%s\t%d\n",           tmpparasname,  &simu_opts->numOE);
  fscanf(fcfg, "%s\t%d\n",           tmpparasname,  &simu_opts->numDE);

  if (simu_opts->numKD == 0){
    simu_opts->KDID = (int *)calloc(1, sizeof(int));
    fscanf(fcfg, "%s\t%d\n",         tmpparasname,  &simu_opts->KDID[0]);
  }
  else{
    simu_opts->KDID = (int *)calloc(simu_opts->numKD, sizeof(int));
    fscanf(fcfg, "%s", tmpparasname);
    for (i = 0; i < simu_opts->numKD; i++){
      fscanf(fcfg, "%d",                         &simu_opts->KDID[i]);
      // printf("%d\n",  simu_opts->KDID[i]);
    }
  }

  if (simu_opts->numOE == 0){
    simu_opts->OEID = (int *)   calloc(1, sizeof(int));
    simu_opts->OEFD = (double *)calloc(1, sizeof(double));
    fscanf(fcfg, "%s\t%d\n",          tmpparasname,  &simu_opts->OEID[0]);
    fscanf(fcfg, "%s\t%lf\n",         tmpparasname,  &simu_opts->OEFD[0]);
  }
  else{
    simu_opts->OEID = (int *)   calloc(simu_opts->numOE, sizeof(int));
    simu_opts->OEFD = (double *)calloc(simu_opts->numOE, sizeof(double));
    fscanf(fcfg, "%s", tmpparasname);
    for (i = 0; i < simu_opts->numOE; i++){
      fscanf(fcfg, "%d",                         &simu_opts->OEID[i]);
      // printf("%d\n",  simu_opts->KDID[i]);
    }
    fscanf(fcfg, "%s", tmpparasname);
    for (i = 0; i < simu_opts->numOE; i++){
      fscanf(fcfg, "%lf",                        &simu_opts->OEFD[i]);
      // printf("%d\n",  simu_opts->KDID[i]);
    }
  }

  if (simu_opts->numDE == 0){
    simu_opts->DEID = (int *)   calloc(1, sizeof(int));
    simu_opts->DEFD = (double *)calloc(1, sizeof(double));
    fscanf(fcfg, "%s\t%d\n",          tmpparasname,  &simu_opts->DEID[0]);
    fscanf(fcfg, "%s\t%lf\n",         tmpparasname,  &simu_opts->DEFD[0]);
  }
  else{
    simu_opts->DEID = (int *)   calloc(simu_opts->numDE, sizeof(int));
    simu_opts->DEFD = (double *)calloc(simu_opts->numDE, sizeof(double));
    fscanf(fcfg, "%s", tmpparasname);
    for (i = 0; i < simu_opts->numDE; i++){
      fscanf(fcfg, "%d",                         &simu_opts->DEID[i]);
      // printf("%d\n",  simu_opts->KDID[i]);
    }
    fscanf(fcfg, "%s", tmpparasname);
    for (i = 0; i < simu_opts->numDE; i++){
      fscanf(fcfg, "%lf",                        &simu_opts->DEFD[i]);
      // printf("%d\n",  simu_opts->KDID[i]);
    }
  }

  fscanf(fcfg, "%s\t%lf\n",          tmpparasname,  &simu_opts->maxtime);
  fscanf(fcfg, "%s\t%lld\n",         tmpparasname,  &simu_opts->myseed);
  fscanf(fcfg, "%s\t%lf\n",          tmpparasname,  &simu_opts->minP);
  fscanf(fcfg, "%s\t%lf\n",          tmpparasname,  &simu_opts->maxP);
  fscanf(fcfg, "%s\t%lf\n",          tmpparasname,  &simu_opts->minK);
  fscanf(fcfg, "%s\t%lf\n",          tmpparasname,  &simu_opts->maxK);
  fscanf(fcfg, "%s\t%lf\n",          tmpparasname,  &simu_opts->minN);
  fscanf(fcfg, "%s\t%lf\n",          tmpparasname,  &simu_opts->maxN);
  fscanf(fcfg, "%s\t%lf\n",          tmpparasname,  &simu_opts->minF);
  fscanf(fcfg, "%s\t%lf\n",          tmpparasname,  &simu_opts->maxF);

  simu_opts->distname = strdup(tmpparasname2);
  simu_opts->flag     = 0;

  // topology information
  fscanf(fcfg, "%s\t%d\n",           tmpparasname,  &topoinfo->numR);
  fscanf(fcfg, "%s\t%d\n",           tmpparasname,  &topoinfo->numG);

  topoinfo->SourceG         = (int *)    calloc(topoinfo->numR,                           sizeof(int));
  topoinfo->TargetG         = (int *)    calloc(topoinfo->numR,                           sizeof(int));
  topoinfo->TypeR           = (int *)    calloc(topoinfo->numR,                           sizeof(int));
  topoinfo->ParasPos        = (int *)    calloc(topoinfo->numR,                           sizeof(int));
  topoinfo->Gname           = (char**)   calloc(topoinfo->numG,                           sizeof(char *));
  for (i = 0; i < topoinfo->numG; i++){
    topoinfo->Gname[i]      = (char*)    calloc(100,                                      sizeof(char));
  }

  tmprlt->Nstb              = 0;
  tmprlt->numover           = (int *)    calloc(topoinfo->numG,                           sizeof(int));
  tmprlt->numdown           = (int *)    calloc(topoinfo->numG,                           sizeof(int));
  tmprlt->cnt_store         = (int *)    calloc(simu_opts->num_stability,                 sizeof(int));
  tmprlt->y_store           = (double *) calloc(simu_opts->num_ode*topoinfo->numG,        sizeof(double));
  tmprlt->soln              = (double *) calloc(simu_opts->num_stability*topoinfo->numG,  sizeof(double));
  tmprlt->paras             = (double *) calloc(3*topoinfo->numR+2*topoinfo->numG,        sizeof(double));

  for (i = 0; i < topoinfo->numG; i++) {
    // fscanf(fcfg, "%*[^\n]\n", NULL);
    fscanf(fcfg, "%d\t%s\n", &tmpID,  topoinfo->Gname[i]);
    // printf("cfg_%d\t%s\n", tmpID,  topoinfo->Gname[i]);
  }

  for (i = 0; i < topoinfo->numR; i++) {
    fscanf(fcfg, "%d\t%d\t%d\t%d\n",  &tmpID, &tmpSourceG, &tmpTargetG, topoinfo->TypeR + i);
    topoinfo->SourceG[i] = tmpSourceG - 1;
    topoinfo->TargetG[i] = tmpTargetG - 1;
    // printf("%d\t%d\t%d\t%d\n",        tmpID, topoinfo->SourceG[i], topoinfo->TargetG[i], topoinfo->TypeR[i]);
  }

  for(i = 0; i < topoinfo->numG; i++){
    for (j = 0; j < topoinfo->numR; j++){
      if (topoinfo->TargetG[j] == i){ 
        topoinfo->ParasPos[j] = 3*cntP + 2*topoinfo->numG;  // Position of threshold parameters for each regulation.
        cntP = cntP + 1;
      }
    }
  }
}

/*********RACIPE Functions*********/
void run_RACIPE(struct opts *simu_opts, struct topo *topoinfo, struct rlt *tmprlt)
{

  clock_t begin, end;
  double  time_spent = 0.0;
  int     i          = 0;
  int     j          = 0;

  begin = clock();

  // Seed the random value generator
  pcg32_srandom((uint64_t)(((uint64_t)time(NULL) ^ (intptr_t)&printf)*simu_opts->myseed), 54u);

  for (i = 1; i <= simu_opts->num_paras; i++){

    set_parameters(simu_opts, topoinfo, tmprlt);
    
    // printf("RACIPELIB 1514");

    for (j = 0; j < simu_opts->num_ode; j++){
      if (simu_opts->solver == 1) {
        solve_ODE_euler(j, simu_opts, topoinfo, tmprlt);
      }
      else {
        solve_ODE_rk45(j, simu_opts, topoinfo, tmprlt);
      }
    }

    count_state (simu_opts, topoinfo, tmprlt);

    save_model_paras(simu_opts, topoinfo, tmprlt, i);
    save_model_solns(simu_opts, topoinfo, tmprlt, i);
    T_test          (simu_opts, topoinfo, tmprlt, i);

    if (simu_opts->SBML_model == i){
      export_SBML_model(simu_opts, topoinfo, tmprlt, i);
    }

    end = clock();
    time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    if (time_spent/3600 >= simu_opts->maxtime){
      printf("### Warning: Time-out!\n");
      break;
    }
  }

  // Screen printout
  if (simu_opts->Toggle_T_test == 1){
    printf("\n-------------------T_test------------------\n");
    printf("Gene_ID -- Probs_over_T\n");
    for (i = 0; i < topoinfo->numG; i++){
       printf("%d -- %f\n", i+1, (double)tmprlt->numover[i]/(double)(tmprlt->numover[i]+tmprlt->numdown[i]));
    }
  }

  printf("\n-----------------Stability-----------------\n");
  printf("#states -- Count\n");
  for (i = 0; i < simu_opts->num_stability; i++){
     printf("%d -- %d\n", i+1, tmprlt->cnt_store[i]);
  }

  end = clock();
  time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
  printf("---> Actual running time : %f seconds ( %f hours)\n", time_spent, time_spent/3600.0);
  printf("The maximum running time is %f.\n", simu_opts->maxtime);

  release_memory  (simu_opts, topoinfo, tmprlt);
}

void save_model_paras(struct opts *simu_opts, struct topo *topoinfo, struct rlt *tmprlt, int num)
{
  static FILE *f_p = NULL;                  
  char fpname [100] = "";

  int i = 0;

  int cnt = tmprlt->Nstb;

  char KDIDname[100] = "";
  char OEIDname[100] = "";
  char DEIDname[100] = "";

  if (simu_opts->Toggle_f_p == 1) {
    if (f_p == NULL) {
      strcpy(fpname, topoinfo->modelname);
      if (simu_opts->exts == 0){
        if (simu_opts->numKD != 0){
          strcat(fpname, "_KD");
          for (i = 0; i < simu_opts->numKD; i++){
            sprintf(KDIDname, "%d", simu_opts->KDID[i]);
            strcat(fpname, "_");
            strcat(fpname, KDIDname);
          }
        }
        if (simu_opts->numOE != 0){
          strcat(fpname, "_OE");
          for (i = 0; i < simu_opts->numOE; i++){
            sprintf(OEIDname, "%d", simu_opts->OEID[i]);
            strcat(fpname, "_");
            strcat(fpname, OEIDname);
          }
        }
        if (simu_opts->numDE != 0){
          strcat(fpname, "_DE");
          for (i = 0; i < simu_opts->numDE; i++){
            sprintf(DEIDname, "%d", simu_opts->DEID[i]);
            strcat(fpname, "_");
            strcat(fpname, DEIDname);
          }
        }
      }
      strcat(fpname, "_parameters");
      strcat(fpname, ".dat");
      f_p   = fopen(fpname,"w");
    }

    fprintf(f_p, "%d\t%d", num, cnt);
    for (i = 0; i < 3*topoinfo->numR+2*topoinfo->numG; i++){
      fprintf(f_p, "\t%f", tmprlt->paras[i]);
    }
    fprintf(f_p, "\n");

    if (simu_opts->num_paras == num){
      fclose(f_p);
      f_p = NULL;
    }
  }
}

void save_model_solns(struct opts *simu_opts, struct topo *topoinfo, struct rlt *tmprlt, int num)
{

  static FILE **f_s;

  int cnt = tmprlt->Nstb;

  if (num == 1) {
    f_s = (FILE **) calloc (simu_opts->num_stability, sizeof(FILE *));
  }
  
  char fsname[simu_opts->num_stability][100];
  char tmpparasname[100] = "";

  int i = 0;
  int h = 0;
  int h2= 0;

  char KDIDname[100] = "";
  char OEIDname[100] = "";
  char DEIDname[100] = "";

  if (f_s[cnt-1] == NULL) {
    sprintf(tmpparasname, "%d", cnt);
    strcpy (fsname[cnt-1], topoinfo->modelname);
    if (simu_opts->exts == 0){
      if (simu_opts->numKD != 0){
        strcat(fsname[cnt-1], "_KD");
        for (i = 0; i < simu_opts->numKD; i++){
          sprintf(KDIDname, "%d", simu_opts->KDID[i]);
          strcat(fsname[cnt-1], "_");
          strcat(fsname[cnt-1], KDIDname);
        }
      }
      if (simu_opts->numOE != 0){
        strcat(fsname[cnt-1], "_OE");
        for (i = 0; i < simu_opts->numOE; i++){
          sprintf(OEIDname, "%d", simu_opts->OEID[i]);
          strcat(fsname[cnt-1], "_");
          strcat(fsname[cnt-1], OEIDname);
        }
      }
      if (simu_opts->numDE != 0){
        strcat(fsname[cnt-1], "_DE");
        for (i = 0; i < simu_opts->numDE; i++){
          sprintf(DEIDname, "%d", simu_opts->DEID[i]);
          strcat(fsname[cnt-1], "_");
          strcat(fsname[cnt-1], DEIDname);
        }
      }
    }
    strcat (fsname[cnt-1], "_solution");
    strcat (fsname[cnt-1], "_");
    strcat (fsname[cnt-1], tmpparasname);
    strcat (fsname[cnt-1], ".dat");
    f_s[cnt-1] = fopen(fsname[cnt-1],"w");
  }

  fprintf(f_s[cnt-1], "%d\t%d", num, cnt);
    for (h = 0; h < cnt; h++){
        h2 = 1;
        while (h2 <= topoinfo->numG) {
            fprintf(f_s[cnt-1], "\t%f", log2(tmprlt->soln[topoinfo->numG*h + h2 - 1]));
            h2++;
        }
    }
  fprintf(f_s[cnt-1], "\n");

  if (simu_opts->num_paras == num) {
    for (i = 0; i < simu_opts->num_stability; i++){
      if (f_s[i] != NULL){
        fclose(f_s[i]);
      }
    }
    free(f_s);
  }
}

void export_SBML_model (struct opts *simu_opts, struct topo *topoinfo, struct rlt *tmprlt, int modelID)
{
  FILE *f_sbml = NULL;                  
  char fpname [100] = "";

  int h = 0;
  int i = 0;
  int j = 0;

  int cnt = tmprlt->Nstb;

  char KDIDname[100] = "";
  char OEIDname[100] = "";
  char DEIDname[100] = "";
  char modelIDname[100] = "";
  char nstbname[100] = "";

  for (h = 1; h <= cnt; h++){

    // printf("stable state is %d\n", cnt);

    if (f_sbml == NULL) {
      strcpy(fpname, topoinfo->modelname);
      if (simu_opts->exts == 0){
        if (simu_opts->numKD != 0){
          strcat(fpname, "_KD");
          for (i = 0; i < simu_opts->numKD; i++){
            sprintf(KDIDname, "%d", simu_opts->KDID[i]);
            strcat(fpname, "_");
            strcat(fpname, KDIDname);
          }
        }
        if (simu_opts->numOE != 0){
          strcat(fpname, "_OE");
          for (i = 0; i < simu_opts->numOE; i++){
            sprintf(OEIDname, "%d", simu_opts->OEID[i]);
            strcat(fpname, "_");
            strcat(fpname, OEIDname);
          }
        }
        if (simu_opts->numDE != 0){
          strcat(fpname, "_DE");
          for (i = 0; i < simu_opts->numDE; i++){
            sprintf(DEIDname, "%d", simu_opts->DEID[i]);
            strcat(fpname, "_");
            strcat(fpname, DEIDname);
          }
        }
      }
      strcat(fpname, "_sbml");
      sprintf(modelIDname, "%d",  modelID);
      strcat(fpname, "_");
      strcat(fpname, modelIDname);

      sprintf(nstbname, "%d", h);
      strcat(fpname, "_");
      strcat(fpname, nstbname);
      strcat(fpname, ".xml");
      f_sbml   = fopen(fpname,"w");
    }

    // write the header
    fprintf(f_sbml, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
    fprintf(f_sbml, "<sbml level=\"2\" version=\"3\" xmlns=\"http://www.sbml.org/sbml/level2/version3\">\n");
    fprintf(f_sbml, "\t<model name=\"No_%d_model_%d\">\n", modelID, h);

    // write the defined function
    fprintf(f_sbml, "\t\t<listOfFunctionDefinitions>\n");
    fprintf(f_sbml, "\t\t\t<functionDefinition id=\"hillfunction\">\n");
    fprintf(f_sbml, "\t\t\t\t<math xmlns=\"http://www.w3.org/1998/Math/MathML\"\n");
    fprintf(f_sbml, "\t\t\t\t\txmlns:sbml=\"http://www.sbml.org/sbml/level3/version2/core\">\n");
    fprintf(f_sbml, "\t\t\t\t\t<lambda>\n");
    fprintf(f_sbml, "\t\t\t\t\t\t<bvar><ci> x </ci></bvar>\n");
    fprintf(f_sbml, "\t\t\t\t\t\t<bvar><ci> x0 </ci></bvar>\n");
    fprintf(f_sbml, "\t\t\t\t\t\t<bvar><ci> lamda </ci></bvar>\n");
    fprintf(f_sbml, "\t\t\t\t\t\t<bvar><ci> nx </ci></bvar>\n");
    fprintf(f_sbml, "\t\t\t\t\t\t<apply>\n");
    fprintf(f_sbml, "\t\t\t\t\t\t\t<plus/>\n");
    fprintf(f_sbml, "\t\t\t\t\t\t\t<ci> lamda </ci>\n");
    fprintf(f_sbml, "\t\t\t\t\t\t\t<apply>\n");
    fprintf(f_sbml, "\t\t\t\t\t\t\t\t<times/>\n");
    fprintf(f_sbml, "\t\t\t\t\t\t\t\t<apply>\n");
    fprintf(f_sbml, "\t\t\t\t\t\t\t\t\t<minus/>\n");
    fprintf(f_sbml, "\t\t\t\t\t\t\t\t\t<cn>1.0</cn>\n");
    fprintf(f_sbml, "\t\t\t\t\t\t\t\t\t<ci>lamda</ci>\n");
    fprintf(f_sbml, "\t\t\t\t\t\t\t\t</apply>\n");
    fprintf(f_sbml, "\t\t\t\t\t\t\t\t<apply>\n");
    fprintf(f_sbml, "\t\t\t\t\t\t\t\t\t<divide/>\n");
    fprintf(f_sbml, "\t\t\t\t\t\t\t\t\t<cn>1.0</cn>\n");
    fprintf(f_sbml, "\t\t\t\t\t\t\t\t\t<apply>\n");
    fprintf(f_sbml, "\t\t\t\t\t\t\t\t\t\t<plus/>\n");
    fprintf(f_sbml, "\t\t\t\t\t\t\t\t\t\t<cn>1.0</cn>\n");
    fprintf(f_sbml, "\t\t\t\t\t\t\t\t\t\t<apply>\n");
    fprintf(f_sbml, "\t\t\t\t\t\t\t\t\t\t\t<power/>\n");
    fprintf(f_sbml, "\t\t\t\t\t\t\t\t\t\t\t<apply>\n");
    fprintf(f_sbml, "\t\t\t\t\t\t\t\t\t\t\t\t<divide/>\n");
    fprintf(f_sbml, "\t\t\t\t\t\t\t\t\t\t\t\t<ci>x</ci>\n");
    fprintf(f_sbml, "\t\t\t\t\t\t\t\t\t\t\t\t<ci>x0</ci>\n");
    fprintf(f_sbml, "\t\t\t\t\t\t\t\t\t\t\t</apply>\n");
    fprintf(f_sbml, "\t\t\t\t\t\t\t\t\t\t\t<ci>nx</ci>\n");
    fprintf(f_sbml, "\t\t\t\t\t\t\t\t\t\t</apply>\n");
    fprintf(f_sbml, "\t\t\t\t\t\t\t\t\t</apply>\n");
    fprintf(f_sbml, "\t\t\t\t\t\t\t\t</apply>\n");
    fprintf(f_sbml, "\t\t\t\t\t\t\t</apply>\n");
    fprintf(f_sbml, "\t\t\t\t\t\t</apply>\n");
    fprintf(f_sbml, "\t\t\t\t\t</lambda>\n");
    fprintf(f_sbml, "\t\t\t\t</math>\n");
    fprintf(f_sbml, "\t\t\t</functionDefinition>\n");
    fprintf(f_sbml, "\t\t</listOfFunctionDefinitions>\n");

    // write the list of species
    fprintf(f_sbml, "\t\t<listOfSpecies>\n");
    for (i=0; i<topoinfo->numG; i++){
      fprintf(f_sbml, "\t\t\t<species id=\"x%d\"  initialAmount=\"%f\"    name=\"%s\"/>\n", i, tmprlt->soln[topoinfo->numG*(h-1) + i], topoinfo->Gname[i]);
    }
    fprintf(f_sbml, "\t\t</listOfSpecies>\n");

    // write the list of parameters
    fprintf(f_sbml, "\t\t<listOfParameters>\n");
    // production
    for (i=0; i<topoinfo->numG; i++){
      fprintf(f_sbml, "\t\t\t<parameter id=\"g%d\"  value=\"%f\"/>\n", i, tmprlt->paras[i]);
    }
    // degradation
    for (i=0; i<topoinfo->numG; i++){
      fprintf(f_sbml, "\t\t\t<parameter id=\"k%d\"  value=\"%f\"/>\n", i, tmprlt->paras[i+topoinfo->numG]);
    }
    // Threshold
    for (i = 0; i < topoinfo->numR; i++){
      fprintf(f_sbml, "\t\t\t<parameter id=\"T%d\"  value=\"%f\"/>\n", 3*i + 2*topoinfo->numG, tmprlt->paras[3*i + 2*topoinfo->numG]);
    }

    // Coefficient
    for (i = 0; i < topoinfo->numR; i++){
      fprintf(f_sbml, "\t\t\t<parameter id=\"n%d\"  value=\"%f\"/>\n", 3*i + 2*topoinfo->numG + 1, tmprlt->paras[3*i + 2*topoinfo->numG + 1]);
    }

    // lambda
    for (i = 0; i < topoinfo->numR; i++){
        if (topoinfo->prsrandrange[2][3*i + 2 + 2*topoinfo->numG] == 1) {      // Activation
          fprintf(f_sbml, "\t\t\t<parameter id=\"lambda%d\"  value=\"%f\"/>\n", 3*i + 2*topoinfo->numG + 2, tmprlt->paras[3*i + 2*topoinfo->numG + 2]);
        }
        else if (topoinfo->prsrandrange[2][3*i + 2 + 2*topoinfo->numG] == 2) { // Inhibition
          fprintf(f_sbml, "\t\t\t<parameter id=\"lambda%d\"  value=\"%f\"/>\n", 3*i + 2*topoinfo->numG + 2, tmprlt->paras[3*i + 2*topoinfo->numG + 2]);
        }
    }

    fprintf(f_sbml, "\t\t</listOfParameters>\n");

    // write the list of reastions
    fprintf(f_sbml, "\t\t<listOfReactions>\n");
    for(i = 0; i < topoinfo->numG; i++){
      fprintf(f_sbml, "\t\t\t<reaction id=\"eq_%d\">\n", i+1);
      fprintf(f_sbml, "\t\t\t\t<listOfProducts>\n");
      fprintf(f_sbml, "\t\t\t\t\t<speciesReference species=\"x%d\" />\n", i);
      fprintf(f_sbml, "\t\t\t\t</listOfProducts>\n");
                  

      fprintf(f_sbml, "\t\t\t\t<kineticLaw>\n");
      fprintf(f_sbml, "\t\t\t\t\t<math xmlns=\"http://www.w3.org/1998/Math/MathML\">\n");
      fprintf(f_sbml, "\t\t\t\t\t\t<apply>\n");
      fprintf(f_sbml, "\t\t\t\t\t\t\t<minus/>\n");
      fprintf(f_sbml, "\t\t\t\t\t\t\t<apply>\n");
      fprintf(f_sbml, "\t\t\t\t\t\t\t\t<times/>\n");

      // production
      fprintf(f_sbml, "\t\t\t\t\t\t\t\t<ci>g%d</ci>\n", i);
                    
      // regulation
      for (j = 0; j < topoinfo->numR; j++){
        if (topoinfo->TargetG[j] == i){
          if (topoinfo->TypeR[j] == 1){ //Activation       
            // fprintf(f_model, "*(Hillshift(ytmp(%d), p(%d), p(%d), p(%d)/p(%d)))", topoinfo->SourceG[j]+1, 2*topoinfo->numG+3*(count-1)+1, 2*topoinfo->numG+3*(count-1)+1+1, 2*topoinfo->numG+3*(count-1)+2+1, 2*topoinfo->numG+3*(count-1)+2+1);
            fprintf(f_sbml, "\t\t\t\t\t\t\t\t<apply>\n");
            fprintf(f_sbml, "\t\t\t\t\t\t\t\t\t<divide/>\n");
            fprintf(f_sbml, "\t\t\t\t\t\t\t\t\t<apply>\n");
            fprintf(f_sbml, "\t\t\t\t\t\t\t\t\t\t<ci>hillfunction</ci>\n");
            fprintf(f_sbml, "\t\t\t\t\t\t\t\t\t\t<ci>x%d</ci>\n", topoinfo->SourceG[j]);
            fprintf(f_sbml, "\t\t\t\t\t\t\t\t\t\t<ci>T%d</ci>\n", 2*topoinfo->numG+3*j); //threshold
            fprintf(f_sbml, "\t\t\t\t\t\t\t\t\t\t<ci>lambda%d</ci>\n", 2*topoinfo->numG+3*j+2); //fold change
            fprintf(f_sbml, "\t\t\t\t\t\t\t\t\t\t<ci>n%d</ci>\n", 2*topoinfo->numG+3*j+1); //coefficient
            fprintf(f_sbml, "\t\t\t\t\t\t\t\t\t</apply>\n");
            fprintf(f_sbml, "\t\t\t\t\t\t\t\t\t<ci>lambda%d</ci>\n", 2*topoinfo->numG+3*j+2); //fold change
            fprintf(f_sbml, "\t\t\t\t\t\t\t\t</apply>\n");
          }
          else if (topoinfo->TypeR[j] == 2) { //Inhibition
            // fprintf(f_model, "*Hillshift(ytmp(%d), p(%d), p(%d), p(%d))",        topoinfo->SourceG[j]+1, 2*topoinfo->numG+3*(count-1)+1, 2*topoinfo->numG+3*(count-1)+1+1, 2*topoinfo->numG+3*(count-1)+2+1);
            fprintf(f_sbml, "\t\t\t\t\t\t\t\t<apply>\n");
            fprintf(f_sbml, "\t\t\t\t\t\t\t\t\t<ci>hillfunction</ci>\n");
            fprintf(f_sbml, "\t\t\t\t\t\t\t\t\t<ci>x%d</ci>\n", topoinfo->SourceG[j]);
            fprintf(f_sbml, "\t\t\t\t\t\t\t\t\t<ci>T%d</ci>\n", 2*topoinfo->numG+3*j); //threshold
            fprintf(f_sbml, "\t\t\t\t\t\t\t\t\t<ci>lambda%d</ci>\n", 2*topoinfo->numG+3*j+2); //fold change
            fprintf(f_sbml, "\t\t\t\t\t\t\t\t\t<ci>n%d</ci>\n", 2*topoinfo->numG+3*j+1); //coefficient
            fprintf(f_sbml, "\t\t\t\t\t\t\t\t</apply>\n");
          }
        }
        
      }

      fprintf(f_sbml, "\t\t\t\t\t\t\t</apply>\n");

      // degradation
      fprintf(f_sbml, "\t\t\t\t\t\t\t<apply>\n");
      fprintf(f_sbml, "\t\t\t\t\t\t\t<times/>\n");
      fprintf(f_sbml, "\t\t\t\t\t\t\t\t<ci>k%d</ci>\n", i);
      fprintf(f_sbml, "\t\t\t\t\t\t\t\t<ci>x%d</ci>\n", i);
      fprintf(f_sbml, "\t\t\t\t\t\t\t</apply>\n");
      fprintf(f_sbml, "\t\t\t\t\t\t</apply>\n");


      fprintf(f_sbml, "\t\t\t\t\t</math>\n");
      fprintf(f_sbml, "\t\t\t\t</kineticLaw>\n");

      fprintf(f_sbml, "\t\t\t</reaction>\n");
    }

    fprintf(f_sbml, "\t\t</listOfReactions>\n");

    // close the blocks
    fprintf(f_sbml, "\t</model>\n");
    fprintf(f_sbml, "</sbml>\n");

    fclose(f_sbml);
    f_sbml = NULL;
  }

}

void set_parameters (struct opts *simu_opts, struct topo *topoinfo, struct rlt *tmprlt)
{ 

  int i = 0;

  switch (simu_opts->dist) {
    case 1:
      // Production rate
      for (i = 0; i < topoinfo->numG; i++){
          tmprlt->paras[i] = randu(topoinfo->prsrandrange[0][i], topoinfo->prsrandrange[1][i]);
      }

      // Degradation rate
      for (i = 0; i < topoinfo->numG; i++){
          tmprlt->paras[i + topoinfo->numG] = randu(topoinfo->prsrandrange[0][i + topoinfo->numG], topoinfo->prsrandrange[1][i + topoinfo->numG]);
      }

      // Threshold
      for (i = 0; i < topoinfo->numR; i++){
          tmprlt->paras[3*i + 2*topoinfo->numG] = randu(topoinfo->prsrandrange[0][3*i + 2*topoinfo->numG], topoinfo->prsrandrange[1][3*i + 2*topoinfo->numG]);
      }

      // Coefficient
      for (i = 0; i < topoinfo->numR; i++){
          tmprlt->paras[3*i + 1 + 2*topoinfo->numG] = randd(topoinfo->prsrandrange[0][3*i + 1 + 2*topoinfo->numG], topoinfo->prsrandrange[1][3*i + 1 + 2*topoinfo->numG], simu_opts->dist);
      }

      // lambda
      for (i = 0; i < topoinfo->numR; i++){
          if (topoinfo->prsrandrange[2][3*i + 2 + 2*topoinfo->numG] == 1) {      // Activation
            tmprlt->paras[3*i + 2 + 2*topoinfo->numG] =     randu(topoinfo->prsrandrange[0][3*i + 2 + 2*topoinfo->numG], topoinfo->prsrandrange[1][3*i + 2 + 2*topoinfo->numG]);
          }
          else if (topoinfo->prsrandrange[2][3*i + 2 + 2*topoinfo->numG] == 2) { // Inhibition
            tmprlt->paras[3*i + 2 + 2*topoinfo->numG] = 1.0/randu(topoinfo->prsrandrange[0][3*i + 2 + 2*topoinfo->numG], topoinfo->prsrandrange[1][3*i + 2 + 2*topoinfo->numG]);
          }
      }
      break;
    case 2:
      // Production rate
      for (i = 0; i < topoinfo->numG; i++){
          tmprlt->paras[i] = randpg(topoinfo->prsrandrange[0][i], topoinfo->prsrandrange[1][i]);
      }

      // Degradation rate
      for (i = 0; i < topoinfo->numG; i++){
          tmprlt->paras[i + topoinfo->numG] = randpg(topoinfo->prsrandrange[0][i + topoinfo->numG], topoinfo->prsrandrange[1][i + topoinfo->numG]);
      }

      // Threshold
      for (i = 0; i < topoinfo->numR; i++){
          tmprlt->paras[3*i + 2*topoinfo->numG] = randpg(topoinfo->prsrandrange[0][3*i + 2*topoinfo->numG], topoinfo->prsrandrange[1][3*i + 2*topoinfo->numG]);
      }

      // Coefficient
      for (i = 0; i < topoinfo->numR; i++){
          tmprlt->paras[3*i + 1 + 2*topoinfo->numG] = randd(topoinfo->prsrandrange[0][3*i + 1 + 2*topoinfo->numG], topoinfo->prsrandrange[1][3*i + 1 + 2*topoinfo->numG], simu_opts->dist);
      }

      // lambda
      for (i = 0; i < topoinfo->numR; i++){
          
          if (topoinfo->prsrandrange[2][3*i + 2 + 2*topoinfo->numG] == 1) {      // Activation
            tmprlt->paras[3*i + 2 + 2*topoinfo->numG] =     randfd(topoinfo->prsrandrange[0][3*i + 2 + 2*topoinfo->numG], topoinfo->prsrandrange[1][3*i + 2 + 2*topoinfo->numG], simu_opts->dist);
          }
          else if (topoinfo->prsrandrange[2][3*i + 2 + 2*topoinfo->numG] == 2) { // Inhibition
            tmprlt->paras[3*i + 2 + 2*topoinfo->numG] = 1.0/randfd(topoinfo->prsrandrange[0][3*i + 2 + 2*topoinfo->numG], topoinfo->prsrandrange[1][3*i + 2 + 2*topoinfo->numG], simu_opts->dist);
          }
      }
      break;
    case 3:
      // Production rate
      for (i = 0; i < topoinfo->numG; i++){
          tmprlt->paras[i] = randexp(topoinfo->prsrandrange[0][i]);
      }

      // Degradation rate
      for (i = 0; i < topoinfo->numG; i++){
          tmprlt->paras[i + topoinfo->numG] = randexp(topoinfo->prsrandrange[0][i + topoinfo->numG]);
      }

      // Threshold
      for (i = 0; i < topoinfo->numR; i++){
          tmprlt->paras[3*i + 2*topoinfo->numG] = randexp(topoinfo->prsrandrange[0][3*i + 2*topoinfo->numG]);
      }

      // Coefficient
      for (i = 0; i < topoinfo->numR; i++){
          tmprlt->paras[3*i + 1 + 2*topoinfo->numG] = randd(topoinfo->prsrandrange[0][3*i + 1 + 2*topoinfo->numG], topoinfo->prsrandrange[1][3*i + 1 + 2*topoinfo->numG], simu_opts->dist);
      }

      // lambda
      for (i = 0; i < topoinfo->numR; i++){
          
          if (topoinfo->prsrandrange[2][3*i + 2 + 2*topoinfo->numG] == 1) {      // Activation
            tmprlt->paras[3*i + 2 + 2*topoinfo->numG] =     randfd(topoinfo->prsrandrange[0][3*i + 2 + 2*topoinfo->numG], 0, simu_opts->dist);
          }
          else if (topoinfo->prsrandrange[2][3*i + 2 + 2*topoinfo->numG] == 2) { // Inhibition
            tmprlt->paras[3*i + 2 + 2*topoinfo->numG] = 1.0/randfd(topoinfo->prsrandrange[0][3*i + 2 + 2*topoinfo->numG], 0, simu_opts->dist);
          }
      }
      break;
  }
  
  if (simu_opts->numKD != 0) {
    for (i = 0; i < simu_opts->numKD; i++){
      if (simu_opts->KDID[i] <= topoinfo->numG) {
        tmprlt->paras[simu_opts->KDID[i] - 1] = 0.0;
      }
      else {
        tmprlt->paras[topoinfo->ParasPos[simu_opts->KDID[i] - topoinfo->numG - 1] + 2] = 1.0;
      }
    }
  }

  if (simu_opts->numOE != 0) {
    for (i = 0; i < simu_opts->numOE; i++){
      if (simu_opts->OEID[i] <= topoinfo->numG) {
        // printf("%d\t%f\t", simu_opts->OEID[i], tmprlt->paras[simu_opts->OEID[i] - 1]);
        tmprlt->paras[simu_opts->OEID[i] - 1] = tmprlt->paras[simu_opts->OEID[i] - 1]*simu_opts->OEFD[i];
        // printf("%f\n", tmprlt->paras[simu_opts->OEID[i] - 1]);
      }
      else {
        printf("### Wrong: no gene (#%d) is selected for overexpression.\n", simu_opts->OEID[i]);
      }
    }
  }

  if (simu_opts->numDE != 0) {
    for (i = 0; i < simu_opts->numDE; i++){
      if (simu_opts->DEID[i] <= topoinfo->numG) {
        // printf("%d\t%f\t", simu_opts->DEID[i], tmprlt->paras[simu_opts->DEID[i] - 1]);
        tmprlt->paras[simu_opts->DEID[i] - 1] = tmprlt->paras[simu_opts->DEID[i] - 1]/simu_opts->DEFD[i];
        // printf("%f\n", tmprlt->paras[simu_opts->DEID[i] - 1]);
      }
      else {
        printf("### Wrong: no gene (#%d) is selected for overexpression.\n", simu_opts->DEID[i]);
      }
    }
  }
}

void model_ODE(double t,  double *ytmp, double *yp, double *p, struct topo *topoinfo)
{
  int i     = 0;

  for(i = 0; i < topoinfo->numG; i++){
    yp[i] = p[i];
  }

  for (i = 0; i < topoinfo->numR; i++){
    if (topoinfo->TypeR[i] == 1){  //Activation
      yp[topoinfo->TargetG[i]] = yp[topoinfo->TargetG[i]]*(Hillshift(ytmp[topoinfo->SourceG[i]], p[topoinfo->ParasPos[i]], p[topoinfo->ParasPos[i]+1], p[topoinfo->ParasPos[i]+2])/p[topoinfo->ParasPos[i]+2]);
    }
    else if (topoinfo->TypeR[i] == 2){  //Inhibition
      yp[topoinfo->TargetG[i]] = yp[topoinfo->TargetG[i]]*Hillshift(ytmp[topoinfo->SourceG[i]], p[topoinfo->ParasPos[i]], p[topoinfo->ParasPos[i]+1], p[topoinfo->ParasPos[i]+2]);
    }
  }

  for(i = 0; i < topoinfo->numG; i++){
    yp[i] = yp[i] - p[i + topoinfo->numG]*ytmp[i];
  }
}

void RIVs(double *y, double *ytmp, double *p, struct topo *topoinfo)
{
  int    i         = 0;
  int    j         = 0;
  double minV      = 0;
  double maxV      = 0;

  // printf("test_RIVs\n");

  for(i = 0; i < topoinfo->numG; i++){

    if (p[i] != 0){
      minV = p[i];
      maxV = p[i];

      for (j = 0; j < topoinfo->numR; j++){
        if (topoinfo->TargetG[j] == i){
          if (topoinfo->TypeR[j] == 1){ //Activation  
            minV = minV*(1.0/p[topoinfo->ParasPos[j] + 2]); 
            // printf("%f\t", p[topoinfo->ParasPos[j] + 2]);
          }
          else if (topoinfo->TypeR[j] == 2) { //Inhibition
            minV = minV*(p[topoinfo->ParasPos[j] + 2]);
            // printf("%f\t", p[topoinfo->ParasPos[j] + 2]);
          }       
        }
      }
      // printf("\n");

      minV = minV/p[topoinfo->numG + i];

      maxV = maxV/p[topoinfo->numG + i];

      y[i] = exp2(randu(log2(minV), log2(maxV)));
    }
    else {
      y[i]    = 0.0;
      ytmp[i] = 0.0;
    }
  }

  // printf("--------------------------------\n");
  // printf("\n");
}

// solve the ODE by 1st Euler method
void solve_ODE_euler (int j, struct opts *simu_opts, struct topo *topoinfo, struct rlt *tmprlt)
{

  int    n_step     = 1000;
  int    i_step     = 1;
  int    i          = 0;
  double testdelta  = 0.0;
  double t          = 0.0;
  double t_start    = 0.0;
  double t_stop     = 0.0;
  double *y;
  double *yp;
  double *ytmp;

  y     = (double *)calloc(topoinfo->numG,      sizeof(double));
  yp    = (double *)calloc(topoinfo->numG,      sizeof(double));
  ytmp  = (double *)calloc(topoinfo->numG,      sizeof(double));

  for (i = 0; i < topoinfo->numG; i++){
    ytmp[i] = 2000.0;
  }
 
  int cnt_loop = 0;
  
  RIVs(y, ytmp, tmprlt->paras, topoinfo);
    
  testdelta = sumdelta(y, ytmp, topoinfo->numG);
  
  while (testdelta != 0 && cnt_loop < simu_opts->maxiters) {
    t_start = t_stop;
    t_stop  = t_stop + 100;
    
    cnt_loop = cnt_loop + 1;
  
    for ( i_step = 1; i_step <= n_step; i_step++ )
    { 
      for (i = 0; i < topoinfo->numG; i++){
          ytmp[i] = y[i];
      }
        
      model_ODE ( t, ytmp, yp, tmprlt->paras, topoinfo );
      
      for (i = 0; i < topoinfo->numG; i++){
              y[i] = ytmp[i] + yp[i]*simu_opts->stepsize;
          }
          
          t = t + simu_opts->stepsize;
      }
      
      testdelta = sumdelta(y, ytmp, topoinfo->numG);
    }
  
  for (i = 0; i < topoinfo->numG; i++){
    tmprlt->y_store[topoinfo->numG*j + i] = y[i];
  }
  
  free(y);
  free(yp);
  free(ytmp);
}

// solve the ODE by RK-45
void solve_ODE_rk45 (int j, struct opts *simu_opts, struct topo *topoinfo, struct rlt *tmprlt)
{
  double abserr;
  double relerr;
  int    flag;
  int    n_step     = 100;
  int    i_step     = 1;
  int    i          = 0;
  double testdelta  = 0.0;
  double t          = 0.0;
  double t_out      = 0.0;
  double t_start    = 0.0;
  double t_stop     = 0.0;
  double *y;
  double *yp;
  double *ytmp;

  abserr = r4_epsilon ( ); //sqrt ( r4_epsilon ( ) );
  relerr = r4_epsilon ( ); //sqrt ( r4_epsilon ( ) );

  y     = (double *)calloc(topoinfo->numG,      sizeof(double));
  yp    = (double *)calloc(topoinfo->numG,      sizeof(double));
  ytmp  = (double *)calloc(topoinfo->numG,      sizeof(double));

  for (i = 0; i < topoinfo->numG; i++){
    ytmp[i] = 2000.0;
  }
 
  int cnt_loop = 0;
  
  RIVs(y, ytmp, tmprlt->paras, topoinfo);
    
  testdelta = sumdelta(y, ytmp, topoinfo->numG);
  
  while (testdelta != 0 && cnt_loop < simu_opts->maxiters) {
    t_start = t_stop;
    t_stop  = t_stop + 100;
    
    cnt_loop = cnt_loop + 1;

    model_ODE ( t, ytmp, yp, tmprlt->paras, topoinfo);
    flag = 1;
  
    // printf("RACIPELIB 2008");

    for ( i_step = 1; i_step <= n_step; i_step++ )
    { 
        for (i = 0; i < topoinfo->numG; i++){
            ytmp[i] = y[i];
        }
    
        t = ( ( double ) ( n_step - i_step + 1 ) * t_start  
            + ( double ) (          i_step - 1 ) * t_stop )
            / ( double ) ( n_step              );

        t_out = ( ( double ) ( n_step - i_step ) * t_start  
                + ( double ) (          i_step ) * t_stop )
                / ( double ) ( n_step          );

        flag = r4_rkf45 (model_ODE, topoinfo->numG, y, yp, tmprlt->paras, topoinfo, &t, t_out, &relerr, abserr, flag );
        // printf ( "%4d  %12f\n", flag, t); 

        // if (i_step == n_step - 1) {
        //     printf ( "%4d  %12f  %12f  %12f\n", flag, t, y[0], y[1]); 
        // }
    }

    // printf ( "%4d  %12f  %12f  %12f\n", flag, t, y[0], y[1]);
      
    testdelta = sumdelta(y, ytmp, topoinfo->numG);
  }
  
  for (i = 0; i < topoinfo->numG; i++){
    tmprlt->y_store[topoinfo->numG*j + i] = y[i];
  }
  
  free(y);
  free(yp);
  free(ytmp);
}


void count_state (struct opts *simu_opts, struct topo *topoinfo, struct rlt *tmprlt)
// detect the stable states from all the solutions from different RIVs
{
    int i         = 0;
    int j         = 0;
    int h         = 0;
    int count     = 0;
    int cnt       = 1;
    double delta  = 0.0;
    double sumpow = 0.0;
    
    for (h = 1; h <= topoinfo->numG; h++){
        tmprlt->soln[h - 1] = tmprlt->y_store[h - 1]; 
    }
    
    for (i = 2; i <= simu_opts->num_ode; i++){
        count = 0;
        
        for (j = 1; j <= cnt; j++){
            h = 1;
            sumpow = 0.0;
            while (h <= topoinfo->numG) {
                sumpow = sumpow + pow((tmprlt->y_store[topoinfo->numG*(i-1) + h - 1] - tmprlt->soln[topoinfo->numG*(j-1) + h - 1]), 2);
                h++;
            }
        
            delta = sqrt(sumpow);
            
            if (delta > simu_opts->thrd){
                count = count + 1;
            } 
        }    
        
        if (count == cnt){
            cnt = cnt + 1;
            
            if (cnt <= simu_opts->num_stability) {
                for (h = 1; h <= topoinfo->numG; h++){
                    tmprlt->soln[(cnt-1)*topoinfo->numG + h - 1] = tmprlt->y_store[topoinfo->numG*(i-1) + h - 1]; 
                }
            }
            else{
                cnt = simu_opts->num_stability;
                break;
            }
        }
    }
    
    tmprlt->cnt_store[cnt-1] = tmprlt->cnt_store[cnt-1] + 1;
    tmprlt->Nstb = cnt;
}

void T_test(struct opts *simu_opts, struct topo *topoinfo, struct rlt *tmprlt, int num)
{
  static FILE *f_test = NULL;
  char ftname [100] = "";

  int i = 0;
  int j = 0;
  int h = 0;
  double tmp = 0.0;
  int *localnumover;
  int *localnumdown;

  int cnt = tmprlt->Nstb;

  char KDIDname[100] = "";
  char OEIDname[100] = "";
  char DEIDname[100] = "";

  if (simu_opts->Toggle_T_test == 1) {
    if (f_test == NULL) {
      strcpy(ftname, topoinfo->modelname);
      if (simu_opts->exts == 0){
        if (simu_opts->numKD != 0){
          strcat(ftname, "_KD");
          for (i = 0; i < simu_opts->numKD; i++){
            sprintf(KDIDname, "%d", simu_opts->KDID[i]);
            strcat(ftname, "_");
            strcat(ftname, KDIDname);
          }
        }
        if (simu_opts->numOE != 0){
          strcat(ftname, "_OE");
          for (i = 0; i < simu_opts->numOE; i++){
            sprintf(OEIDname, "%d", simu_opts->OEID[i]);
            strcat(ftname, "_");
            strcat(ftname, OEIDname);
          }
        }
        if (simu_opts->numDE != 0){
          strcat(ftname, "_DE");
          for (i = 0; i < simu_opts->numDE; i++){
            sprintf(DEIDname, "%d", simu_opts->DEID[i]);
            strcat(ftname, "_");
            strcat(ftname, DEIDname);
          }
        }
      }
      strcat(ftname, "_T_test");
      strcat(ftname, ".dat");
      f_test   = fopen(ftname, "w");
    }

    localnumover = (int *) calloc(topoinfo->numG, sizeof(int));
    localnumdown = (int *) calloc(topoinfo->numG, sizeof(int));

    for (i = 0; i < cnt; i++){
      for (j = 0; j < topoinfo->numG; j++){
        for (h = 0; h < topoinfo->numR; h++){
          if (topoinfo->SourceG[h] == j){
            tmp = tmprlt->soln[topoinfo->numG*i + j]/tmprlt->paras[topoinfo->ParasPos[h]];
            if (tmp >= 1){
              localnumover[j] = localnumover[j] + 1;
            }
            else {
              localnumdown[j] = localnumdown[j] + 1;
            }
          }
        }
      }
    }

    fprintf(f_test, "%d", num);
    for (j = 0; j < topoinfo->numG; j++){
      fprintf(f_test, "\t%d\t%d", localnumover[j], localnumdown[j]);
      tmprlt->numover[j] = tmprlt->numover[j] + localnumover[j];
      tmprlt->numdown[j] = tmprlt->numdown[j] + localnumdown[j];
    }
    fprintf(f_test, "\n");

    free(localnumover);
    free(localnumdown);
  
    if (simu_opts->num_paras == num){
      fclose(f_test);
      f_test = NULL;
    }
  
  }
}

void release_memory(struct opts *simu_opts, struct topo *topoinfo, struct rlt *tmprlt)
{
  int i = 0;

  free(topoinfo->SourceG);
  free(topoinfo->TargetG);
  free(topoinfo->TypeR);
  free(topoinfo->ParasPos);
  free(topoinfo->modelname);

  if (simu_opts->exts == 0){
    for (i = 0; i < topoinfo->numG; i++) {
      free(topoinfo->Gname[i]);
    }
    free(topoinfo->Gname);
  }

  free(topoinfo->prsrandrange[0]);
  free(topoinfo->prsrandrange[1]);
  free(topoinfo->prsrandrange[2]);
  free(topoinfo->prsrandrange);
  
  free(tmprlt->numover);
  free(tmprlt->numdown);
  free(tmprlt->cnt_store);
  free(tmprlt->y_store);
  free(tmprlt->soln);
  free(tmprlt->paras);

  free(simu_opts->KDID);
  free(simu_opts->OEID);
  free(simu_opts->OEFD);
  free(simu_opts->DEID);
  free(simu_opts->DEFD);

}