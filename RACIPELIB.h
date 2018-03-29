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

#ifndef RACIPELIB_h
#define RACIPELIB_h

#include <stdio.h>
#include <inttypes.h>

struct topo {
	int    numG;            // Number of genes in the circuit
	int    numR;            // Number of regulations in the circuit
	int    *SourceG;        // Source genes (as Gene ID, start from 0)
	int    *TargetG;        // Target genes (as Gene ID, start from 0)
	int    *TypeR;          // Type of regulations, 1 for activation and 2 for inhibition
    int    *ParasPos;       // Position of threshold parameters for each regulation
	char   *modelname;      // temporary store for the model name
	char   **Gname;         // Gene names
	double **prsrandrange;  // Min (first row) and max (second row) value of the range of parameters for randomization. The last row is the type of regulation (0 for all the parameters for production and degradation, but 1 or 2 for parameters relevant with regulatons)
};

struct rlt {
	int    Nstb;			// number of stable states. 
	int    *numover;        // number of gene expression of a gene to be bigger than its threshold;
    int    *numdown;        // number of gene expression of a gene to be bigger than its threshold;
    int    *cnt_store;      // the number of each stability
    double *y_store;        // y_store to store the solutions for each RIVs
  	double *soln;           // soln to store the solutions for each RACIPE model
    double *paras;          // paras to store the parameters for each RACIPE model
};

struct opts {
	double   maxtime;         // 										  - Maximum running time
	int 	 solver;		  // 										  - solver for ODEs								
	int      flag; 			  //										  - Only produce .cfg file or not
	int      exts; 			  //										  - Extensiton of inputfile: 0 for .topo and 1 for .cfg
	int      numKD; 		  // 										  - Number of genes or links' KD
	int*     KDID;			  // 										  - The genes or links to knockdown	
	int      numOE; 		  // 										  - Number of genes overexpression
	int*     OEID;			  // 										  - The genes or links to overexpression
	double*  OEFD; 			  //										  - The fold change to overexpress a gene
	int      numDE; 		  // 										  - Number of genes downexpression
	int*     DEID;			  // 										  - The genes or links to downexpression
	double*  DEFD; 			  //										  - The fold change to downexpress a gene
	char     *distname;
	int      dist;      	  // "Distribution"            			      - Distribution used for randomization
	double   SF;    		  //                           			      - Scale factor
	int      num_findT;       // "NumberOfSimulationToFindThreshold"      - Number of simulaitons used to estimate threshold;
	int      num_paras;       // "NumberOfRACIPEModels"    			      - Number of RACIPE models to generate
	int      num_ode;      	  // "NumberOfRIVs"            			      - Number of Random initial values to solve ODEs
	int      num_stability;   // "NumberOfStatesToStore"   			      - Number of stable states to count
	double   thrd;    		  // "ThresholdForConvergence" 			      - Cutoff for convergence of states
	int      Toggle_f_p;      // "ToggleOfSaveParameters"  			      - Toggle to save parameter or not
	double   stepsize;    	  // "Stepsize"                			      - Stepsize for solving ODE;
	int      maxiters;        // "MaximumOfIterations"     			      - Maximum of Iteration for solving ODE at each RIVs
	int      Toggle_T_test;   // "TestThreshold"           			      - Toggle to test threshold assumption or not
	int 	 SBML_model;	  //            			      			  - Export which model in SBML format
	uint64_t myseed;          //  										  - Random seed set by user
	double   minP;			  //  										  - Minimum production rate
  	double   maxP;			  //  										  - Maximum production rate
  	double   minK;			  //  										  - Minimum degradation rate
  	double   maxK;			  //  										  - Maximum degradation rate
  	double   minN;			  //  										  - Minimum coefficient
  	double   maxN;			  //  										  - Maximum coefficient
  	double   minF;			  //  										  - Minimum fold change
  	double   maxF;			  //  										  - Maximum fold change
};

extern double randu(double minV, double maxV);
extern double randg(double m, double stdvalue);
extern double randpg(double m, double stdvalue);
extern double randexp(double m);
extern double randfd(double m, double stdvalue, int dist);
extern double randd(double minN, double maxN, int dist);
extern double Hillshift (double x, double x0, double nx, double lamda);
extern double median(double *x, int n);
extern double sumdelta(double *y, double *ytmp, int NEQN);

extern void   check_inputfile (int argc, char **argv, struct topo *topoinfo, struct opts *simu_opts, struct rlt *tmprlt);
extern void   initial_simuopts (int argc, char **argv, struct topo *topoinfo, struct opts *simu_opts);
extern void   Model_generate (char *inputfile, struct topo *topoinfo, struct opts *simu_opts, struct rlt *tmprlt);
extern void   check_topo(FILE *f_in, FILE *f_out, struct topo *topoinfo);
extern void   generate_random_range(FILE *f_paras, struct topo *topoinfo, struct opts *simu_opts);
extern void   estimate_threshold(int num, int ID, double minP, double maxP, double minK, double maxK, double minN, double maxN, double minF, double maxF, double *minT, double *maxT, struct topo *topoinfo, int dist, double SF);
extern void   read_cfg(struct topo *topoinfo, struct opts *simu_opts, struct rlt *tmprlt);

extern void   run_RACIPE(struct opts *simu_opts, struct topo *topoinfo, struct rlt *tmprlt);
extern void   save_model_paras(struct opts *simu_opts, struct topo *topoinfo, struct rlt *tmprlt, int num);
extern void   save_model_solns(struct opts *simu_opts, struct topo *topoinfo, struct rlt *tmprlt, int num);
extern void	  export_SBML_model (struct opts *simu_opts, struct topo *topoinfo, struct rlt *tmprlt, int modelID);
extern void   set_parameters (struct opts *simu_opts, struct topo *topoinfo, struct rlt *tmprlt);
extern void   model_ODE(double t,  double *ytmp, double *yp, double *p, struct topo *topoinfo);
extern void   RIVs(double *y, double *ytmp, double *p, struct topo *topoinfo);
extern void   solve_ODE_euler (int j, struct opts *simu_opts, struct topo *topoinfo, struct rlt *tmprlt);
extern void   solve_ODE_rk45 (int j, struct opts *simu_opts, struct topo *topoinfo, struct rlt *tmprlt);
extern void   count_state (struct opts *simu_opts, struct topo *topoinfo, struct rlt *tmprlt);
extern void   T_test(struct opts *simu_opts, struct topo *topoinfo, struct rlt *tmprlt, int num);
extern void   release_memory(struct opts *simu_opts, struct topo *topoinfo, struct rlt *tmprlt);

#endif
