# Evometrate
Code from: Hin &amp; De Roos Evolution of Size-Dependent Intraspecific Competition Predicts Body Size Scaling of Metabolic Rate. Functional Ecology

Execution of the code requires the installation of PSPManalysis package
See: https://bitbucket.org/amderoos/pspmanalysis/src/master/ for more info

 - Hin&DeRoos_EvoMetrate_analysis.R contains the commands to run the model and only this file should be executed
 - Hin&DeRoos_EvoMetrate_functions.R contains functions and parameter settings and is sourced from (..)_analysis.R
 - Hin&DeRoos_EvoMetrate_plotting.R creates Figures 1 - 4 from the above publication is sourced from (..)_analysis.R
 - Indet_growth_2exp.h is a C-header file that contains the model implmentation and should not be changed.
 - EBT_EQ_P.minmax.txt and EBT_EQ_Q.minmax contain output of numerical simulations and are required to plot Figure1

Code to plot figure 5 of the above publication can be found in Dryad Digital Repository:
https://doi.org/doi:10.5061/dryad.2fc6677

Code is released under GNU Public License, Version 3. A copy of the license is attached.

Please cite the above paper when reusing any part of the code in a publication
 
Last modified: VH - 23 November 2018
