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
# include <stdio.h>
# include <stdint.h>
# include <math.h>
# include <time.h>
# include <string.h>

# include "RACIPELIB.h"
# include "pcg_basic.h"

int main (int argc, char **argv)
{  
  struct  topo   topoinfo;    // store topology relevant information, details in RACIPELIB.h.
  struct  opts   simu_opts;   // store the options of the simulation, details in RACIPELIB.h.
  struct  rlt    tmprlt;      // tmporarily store the results, details in RACIPELIB.h.

  // Preprocess the topology file (.topo) or configure file (.cfg) for RACIPE method
  check_inputfile (argc, argv, &topoinfo, &simu_opts, &tmprlt); 
  // Run RACIPE to produce results
  run_RACIPE      (&simu_opts, &topoinfo, &tmprlt);

  return 0;
}


