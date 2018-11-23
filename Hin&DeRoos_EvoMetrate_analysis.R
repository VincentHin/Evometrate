###########################################################################
#
#    R script to run analysis from the paper:
# 
#    Evolution of Size-Dependent Intraspecific Competition Predicts Body Size Scaling of Metabolic Rate
#    Hin, V & De Roos, A.M.
#    Functional Ecology
# 
#    Code is released under GNU Public License, Version 3
#    Please cite the above paper when reusing any part of the code in a publication
# 
#    Execution of the code requires the installation of PSPManalysis package
#    See: https://bitbucket.org/amderoos/pspmanalysis/src/master/ for more info
#
#    Code to plot figure 5 can be found in Dryad Digital Repository:
#    https://doi.org/doi:10.5061/dryad.2fc6677
#
#    Last modified: VH - 23 November 2018
#
###########################################################################
# Install PSPManalysis:
devtools::install_bitbucket("amderoos/PSPManalysis", subdir = "R/")

# Load packages and functions
library('plyr'); library('PSPManalysis'); library('Hmisc')
funcs <- new.env()
source("Hin&DeRoos_EvoMetrate_functions.R", local = funcs)
attach(funcs); rm(funcs)

Q <- seq(0.3, 1.5, 0.01)
P <- seq(0.3, 1.5, 0.01)
RmaintQ_XB <- 0.3 * 0.1^1 / (0.5*0.1^Q - 0.1*0.1^1)
RmaintQ_XM <- 0.3 * 10^1 / (0.5 * 10^Q - 0.1*10^1)
RmaintP_XB <- 0.3 * 0.1^P / (0.5 * 0.1^1 - 0.1*0.1^P)
RmaintP_XM <- 0.3 * 10^P / (0.5 * 10^1 - 0.1*10^P)

## Calculate the equilibrium point for default parameter value
eq.pnt <- equipoint(bifpar = "RMAX")

## Calculate equilibrium curve as a function (Figure 1)
## of Q (max. ingestion exponent)
eq.Q <- equicurve(bifpar = 'Q', inipoint = eq.pnt$curvepoints, par1 = c(0,1.15,0.1))
eq.Q$curvepoints <- cal_rates(eq.Q$curvepoints)
## and P (maintenance rate exponent)
eq.P <- equicurve(bifpar = 'P', inipoint = eq.pnt$curvepoints, par1 = c(0.85,2.0,0.1))
eq.P$curvepoints <- cal_rates(eq.P$curvepoints)

## Calculate the Evolutionary Isoclines in the Phase Plane of Q and P (Figure 2)
ess.QP <- ess.scan.2D(bifpar1 = 'Q', bifpar2 = 'P', 
                      ESSpnt1 = eq.Q$bifpoints, 
                      ESSpnt2 = eq.P$bifpoints,
                      par1 = c(0,2,0.1), 
                      par2 = c(0,2,0.1))

## The 2-dimensional CSS of Q and P is detected:
CSS.QP <- ess.QP$bifpoints1[[1]][c(2:4,1)]

## --- ESS (Q & P) ~ SB ---- ##
pars <- pars.default
ess.QP_SB <- twodim.ess.cont(bifpar = 'XB', parbnds = c(5E-3,0.5,0.1), 
                             ESSpar1 = 'Q', ESSpar2 = 'P', 
                             ESSpnt = ess.QP$bifpoints1[[1]], 
                             ESSbnds1 = c(0.4,1.3), ESSbnds2 = c(0.4,1.3))
ess.QP_SB$curvepoints <- cal_rates(ess.QP_SB$curvepoints)
ess.QP_SB$curvedesc

## ---- ESS (Q & P) ~ SM ---- ##
pars <- pars.default
ess.QP_SM <- twodim.ess.cont(bifpar = 'XM', parbnds = c(1,100,0.1), 
                             ESSpar1 = 'Q', ESSpar2 = 'P', 
                             ESSpnt = ess.QP$bifpoints1[[1]], 
                             ESSbnds1 = c(0.4,1.3), ESSbnds2 = c(0.4,1.3))
ess.QP_SM$curvepoints <- cal_rates(ess.QP_SM$curvepoints)
ess.QP_SM$curvedesc
##

## --- ESS (Q & P) ~ MUJ ---- ##
ess.QP_MUJ <- twodim.ess.cont(bifpar = 'MUJ', 
                              parbnds = c(0,1,0.1), 
                              ESSpar1 = 'Q', ESSpar2 = 'P', 
                              ESSpnt = ess.QP$bifpoints1[[1]], 
                              ESSbnds1 = c(0,2), ESSbnds2 = c(0,2))
ess.QP_MUJ$curvepoints <- cal_rates(ess.QP_MUJ$curvepoints)
ess.QP_MUJ$curvedesc

## ---- EQ ~ MUJ | CSS (Q & P) ---- ##
pars <- setpars$pars.default
pars[setpars$pars.names == 'Q'] <- as.numeric(ess.QP$bifpoints1[[1]][4])
pars[setpars$pars.names == 'P'] <- as.numeric(ess.QP$bifpoints1[[1]][1])
eq.MUJ <- equicurve(bifpar = "MUJ", inipoint = eq.pnt$curvepoints, p = pars)
eq.MUJ$curvepoints <- cal_rates(eq.MUJ$curvepoints)
eq.MUJ$curvedesc

## --- ESS (Q & P) ~ MUJ ---- ##
ess.QP_MUJ <- twodim.ess.cont(bifpar = 'MUJ',
                              parbnds = c(0,0.4,0.1), 
                              ESSpar1 = 'Q', ESSpar2 = 'P', 
                              ESSpnt = ess.QP$bifpoints1[[1]], 
                              ESSbnds1 = c(0,2), ESSbnds2 = c(0,2))
ess.QP_MUJ$curvepoints <- cal_rates(ess.QP_MUJ$curvepoints)
ess.QP_MUJ$curvedesc

## ---- ESS (Q & P) | SB = 0.05 ---- ##
pars <- setpars$pars.default
pars[setpars$pars.names == 'XB'] <- 0.05

eq.pnt_SB0.05 <- equipoint(bifpar = "RMAX", p = pars)
eq.Q_SB0.05 <- equicurve(bifpar = 'Q', 
                         inipoint = eq.pnt_SB0.05$curvepoints, 
                         p = pars,
                         par1 = c(0.4,1.15,0.1))
eq.P_SB0.05 <- equicurve(bifpar = 'P', 
                         inipoint = eq.pnt_SB0.05$curvepoints,
                         p = pars,
                         par1 = c(0.86,2,0.1))

ess.QP_SB0.05 <- ess.scan.2D(bifpar1 = 'Q', bifpar2 = 'P', 
                             ESSpnt1 = eq.Q_SB0.05$bifpoints[1,], 
                             ESSpnt2 = eq.P_SB0.05$bifpoints,
                             p = pars,
                             par1 = c(0.4,1.5,0.1), 
                             par2 = c(0.4,1.5,0.1))
##

## ---- EQ ~ MUA | ESS (Q & P) & SB = 0.05 ---- ##
pars <- setpars$pars.default
pars[setpars$pars.names == 'XB'] <- 0.05
pars[setpars$pars.names == 'P'] <- as.numeric(ess.QP_SB0.05$bifpoints1[[1]][1])
pars[setpars$pars.names == 'Q'] <- as.numeric(ess.QP_SB0.05$bifpoints1[[1]][4])

eq.MUA_SB0.05 <- equicurve(bifpar = "MUA", 
                           inipoint = eq.pnt$curvepoints, 
                           p = pars, 
                           par1 = c(0,0.04,0.1))
eq.MUA_SB0.05$curvepoints <- cal_rates(eq.MUA_SB0.05$curvepoints)
eq.MUA_SB0.05$curvedesc

## --- ESS (Q & P) ~ MUA | SB = 0.05 ---- ##
ess.QP_MUA_SB0.05 <- twodim.ess.cont(bifpar = 'MUA', parbnds = c(0,1,0.1),
                                     ESSpar1 = 'Q', 
                                     ESSpar2 = 'P', 
                                     ESSpnt = ess.QP_SB0.05$bifpoints1[[1]],
                                     ESSbnds1 = c(0.4,1.3), 
                                     ESSbnds2 = c(0.4,1.3), 
                                     p = pars)
ess.QP_MUA_SB0.05$curvepoints <- cal_rates(ess.QP_MUA_SB0.05$curvepoints)
ess.QP_MUA_SB0.05$curvedesc

# Force clean all output files in working directory:
PSPMclean('F')

## Make all plots by sourcing (...)_plotting.R
source('Hin&DeRoos_EvoMetrate_plotting.R')