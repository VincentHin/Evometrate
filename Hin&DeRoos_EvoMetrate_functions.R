###########################################################################
# This script contains functions and settings that are used by:
#   Hin&DeRoos_EvoMetrate_analysis.R
# and should not be run directly, as it is sourced from the above file
###########################################################################

par.defaults <- par()
par.defaults$cin <- par.defaults$cra <- par.defaults$csi <- par.defaults$cxy <- par.defaults$din <- par.defaults$page <- NULL

# Load simulation results from ebttool
EBT_EQ_Q_minmax <- read.table('EBT_EQ_Q_minmax.txt', header = TRUE)
EBT_EQ_P_minmax <- read.table('EBT_EQ_P_minmax.txt', header = TRUE)

cycle_Q <- 1.03
resource_Q <- 0.9011366
cycle_P <- 0.96
resource_P <- 4.42194

setpars <- list(model.name = 'Indet_growth_2exp',
                env.names = 'resource',
                impact.names = c('ingest', 'juv_bio', 'ad_bio', 'tot_bio', 'birth_size', 'rep_rate', 'juv_prod', 'birth_surv', 'birth_age'),
                ebt.names = c('time', 'resource', 'juv_bio', 'ad_bio',  'total_bio', 'juv_number', 'ad_number', 'total_number', 'ingestion', 'birthrate', 
                              'rep_rate', 'mat_rate', 'age_mat', 'juv_surv',  'max_size', 'max_age', 'bifparone', 'period1', 'period2'),
                bifpar.names = NULL,
                istate.names = 'size',
                pars.default = c(30, 0.1, 1, 1.0, 0.1, 1.0, 0.015, 0, 0, 0.1, 1.0, 1.0, 10.0, 0.5, 3.0),
                pars.names = c('RMAX', 'DELTA', 'M', 'Q', 'T', 'P', 'MU', 'MUJ', 'MUA', 'XB', 'XJ', 'Xref', 'XM', 'SIGMA', 'H'),
                parameters = list(RMAX = 30, DELTA = 0.1, M = 1, Q = 1.0, Tcon = 0.1, P = 1.0, MU = 0.015, MUJ = 0, 
                                  MUA = 0, XB = 0.1, XJ = 1, Xref = 1, XM = 10, SIGMA = 0.5, H = 3))

##
myaxis <- function(y.at = NULL, x.at = NULL, side = 2, lw = 0.5, tck.lw = 1, xlabels = TRUE, ylabels = TRUE){
  box(lwd = lw)
  
  if(is.logical(xlabels) && xlabels == TRUE) xlabels <- x.at * xlabels
  axis(side=1, at = x.at, tcl = -0.2, labels = FALSE, lwd = lw, lwd.ticks = tck.lw)
  axis(side=1, at = xlabels, labels = xlabels, line = -0.3, lwd = 0)
  
  if(is.logical(ylabels) && ylabels == TRUE) ylabels <- y.at * ylabels
  axis(side=side, at = y.at, tcl = -0.2, labels = FALSE, lwd = lw, lwd.ticks = tck.lw)
  axis(side=side, at = ylabels, labels = ylabels, line = -0.3, lwd = 0, las = 2)
}
##

# Function to transform output
cal_rates <- function(x){
  mat_rate = x$juv_prod + x$rep_rate; 
  juv_surv = exp(-(x$birth_surv / x$birthrate));
  max_size = (x$birth_size / x$birthrate) + 1;
  age_mat = (x$birth_age / x$birthrate) + 1;
  cbind(x, mat_rate, juv_surv, max_size, age_mat)
}


# : Names can consist of: bifpar.names, env.names, birthrate, impact.names, R0, R0_x, eig_J, eig_H, eig_J+J, ZT_C01_Z, RHS_norm
Give.Name <- function(type = 'EQ', 
                      evodims = 1, 
                      bifpar.names = setpars$bifpar.names, 
                      env.names = setpars$env.names, 
                      impact.names = setpars$impact.names, 
                      istate.names = setpars$istate.names){
  if (is.null(impact.names)) return('Please specify a vector of impact names, either directly or within a list called setpars')
  
  # Set names of bifurcation parameters.
  bifpar.nms <- c('bifparone', 'bifpartwo', 'bifparthree', 'bifparfour', 'bifparfive')
  
  if (is.null(bifpar.names)){ 
    print('Using default bifpar.names')
  } else {
    bifpar.nms <- c('bifparone', 'bifpartwo', 'bifparthree', 'bifparfour', 'bifparfive')
    bifpar.nms[1:length(bifpar.names)] <- bifpar.names
  }
  
  if (is.null(env.names)){ 
    print('Using default env.names')
    env.nms <- 'resource'
  } else {
    env.nms <- env.names
  }
  
  if (type == "EQ"){
    name = switch((evodims + 1), 
                  c(bifpar.nms[1], env.nms, 'birthrate', impact.names, 'R0', 'RHS_norm'), 
                  c(bifpar.nms[1], env.nms, 'birthrate', impact.names, 'R0', 'R0_x', 'RHS_norm'))
  } else if (type == "ESS"){
    name = switch(evodims, 
                  c(bifpar.nms[1], env.nms, 'birthrate', bifpar.nms[2], impact.names, 'R0', 'R0_x', 'R0_xx', 'R0_yy', 'RHS_norm'),
                  c(bifpar.nms[1], env.nms, 'birthrate', bifpar.nms[2:3], impact.names,'R0', 'R0_x', 'eig_J', 'eig_H', 'eig_J+J_0.5', 'ZT_C01_Z', 'RHS_norm'),
                  c(bifpar.nms[1], env.nms, 'birthrate', bifpar.nms[2:4], impact.names,'R0', 'R0_x', 'eig_J', 'eig_H', 'eig_J+J_0.5', 'ZT_C01_Z', 'RHS_norm'),
                  c(bifpar.nms[1], env.nms, 'birthrate', bifpar.nms[2:5], impact.names,'R0', 'R0_x', 'eig_J', 'eig_H', 'eig_J+J_0.5', 'ZT_C01_Z', 'RHS_norm'))
  } else if (type == "EVO"){
    name = switch(evodims,  
                  c('evotime', env.nms, 'birthrate', bifpar.nms[1], impact.names, 'R0', 'RHS_norm'),
                  c('evotime', env.nms, 'birthrate', bifpar.nms[1:2], impact.names, 'R0', 'RHS_norm'),
                  c('evotime', env.nms, 'birthrate', bifpar.nms[1:3], impact.names, 'R0', 'RHS_norm'),
                  c('evotime', env.nms, 'birthrate', bifpar.nms[1:4], impact.names, 'R0', 'RHS_norm'),
                  c('evotime', env.nms, 'birthrate', bifpar.nms[1:5], impact.names, 'R0', 'RHS_norm'))
  } else if (type == "BP" || type == "LP"){
    name = c(bifpar.nms[1], env.nms, 'birthrate', bifpar.nms[2], impact.names, 'R0', 'RHS_norm')
  } else if (type == "PIP"){
    name = c(bifpar.nms[1], env.nms, 'birthrate', bifpar.nms[2], impact.names, 'R0', 'R0_0', 'R0_1', 'RHS_norm')
  } else if (type == "IND"){
    name = c('survival', istate.names, paste0('cum_',c(impact.names, 'R0')))
  } else {
    name = "INVALID NAME"
    print("Type is not known, choose either EQ, ESS, EVO, BP, LP or PIP")
  }
  cat('Names set to: ', name, '\n\n')
  
  return(name)
}
##

## --------- FUNCTION: EQUIPOINT() --------- ##
## : -> Find an equilibrium point from a boundary equilibrium for the parameters in 'pars' -- ##
equipoint <- function(modeln = setpars$model.name, 
                      trivpar = 'RMAX',
                      bifpar = 'RMAX', 
                      p = setpars$pars.default, 
                      p.names = setpars$pars.names, 
                      par1 = c(0.1, p[which(p.names == bifpar)], 0.1),
                      impact.names = setpars$impact.names, ...){
  if (is.null(bifpar)) return("Please specify a bifurcation parameter for continuation of trivial equilibrium")
  if (length(par1) != 3) return("Par1 vector is not of length three: minimum parameter value, maximum parameter value, stepsize")
  if (length(p) != length(p.names)) return("Number of parameters is unequal to number of parameter names")
  
  # Set minpar, maxpar and stepsize
  minpar <- par1[1]
  maxpar <- par1[2]
  stepsize <- par1[3]
  
  if(bifpar == trivpar) p[p.names == bifpar] <- minpar
  
  # Calculate trivial equilibrium from zero initial condition
  ini <- c(minpar, p[p.names == trivpar])
  parini <- c(which(p.names == bifpar) - 1, 0.0, maxpar)
  triv <- PSPMequi(modeln, 'EQ', as.numeric(ini), stepsize, parini, p, options = c("popZE", "0"))
  if (is.null(triv$biftypes)) return('No branching point found')
  bound.point <- as.numeric(triv$bifpoints[triv$biftypes == "BP #0"])
  
  # Restore default value of bifurcation parameter
  p[p.names == bifpar] <- maxpar
  
  # Calculate non-trivial equilibrium from branching point
  ini <- c(bound.point[1:2],0)
  parini <- c(which(p.names == bifpar) - 1, 0.0, maxpar)
  nontriv <- PSPMequi(modeln, 'EQ', ini, stepsize, parini, p, ...)
  
  index <- which.min(abs(nontriv$curvepoints[,1] - p[p.names == bifpar]))
  nontriv$curvepoints <- nontriv$curvepoints[index,]
  if(!is.null(impact.names)) names(nontriv$curvepoints) <- Give.Name('EQ', evodims = 0, impact.names = impact.names)
  
  return(nontriv)
}
##

## --------- FUNCTION: EQUICURVE() --------- ##
## : -> Function to calculate bifurcation from parameters in pars for bifurcation parameter bifpar -- ##
equicurve <- function(modeln = setpars$model.name, 
                      bifpar = NULL, 
                      inipoint = NULL, 
                      p = setpars$pars.default, 
                      p.names = setpars$pars.names, 
                      par1 = c(0,100,0.1), ...){
  if (is.null(bifpar)) return("Please specify a bifurcation parameter")
  if(length(which(p.names == bifpar)) == 0) return('Bifpar not found')
  if (is.null(inipoint)) return("Please specify an initial point")
  if (length(par1) != 3) return("Par1 vector is not of length three: minimum parameter value, maximum parameter value, stepsize")
  if (length(p) != length(p.names)) return("Number of parameters is unequal to number of parameter names")
  
  if(is.list(inipoint)){
    cat('Initial point specified as list, looking for $curvepoints...\n\n')
    if(is.numeric(inipoint$curvepoints)){
      cat('Initial point taken from curvepoints element of list \n\n')
      inipoint <- inipoint$curvepoints
    }
  }
  
  # Get some names
  nms <- Give.Name(type = "EQ")
  
  # Set minpar, maxpar and stepsize
  minpar = par1[1]
  maxpar = par1[2]
  stepsize = par1[3]
  
  # Store default values of the bifurcation parameters
  old.bifpar = p[p.names == bifpar]
  
  # Set initial point & parameter boundaries
  env.dim <- length(setpars$env.names)
  sel <- c(2:(2+env.dim))
  start.eq <- inipoint[sel]
  ini <- c(p[p.names == bifpar], start.eq)
  parini <- c(which(p.names == bifpar) - 1, minpar, maxpar)
  
  # Run forwards and backwards
  pbifplus <- PSPMequi(modeln, 'EQ', as.numeric(ini), stepsize, parini, p, options = c("popEVO", "0"), ...)
  pbifmin  <- PSPMequi(modeln, 'EQ', as.numeric(ini), -1*stepsize, parini, p, options = c("popEVO", "0"))
  
  lst <- list('curvepoints' = NULL, 'bifpoints' = NULL, 'curvedesc' = NULL)
  
  # Store all points, sort and select
  pbif <- rbind(pbifmin$curvepoints[nrow(pbifmin$curvepoints):1,], pbifplus$curvepoints)
  pbif[which(pbif == "NaN")] <- 0.0
  pbif <- as.data.frame(pbif)
  names(pbif) <- nms
  
  lst$curvepoints <- pbif
  lst$curvedesc <- pbifplus$curvedesc
  
  bifpoints <- rbind(pbifmin$bifpoints, pbifplus$bifpoints)
  bifpoints[which(bifpoints == "NAN" | bifpoints == "NaN")] <- 0.0
  bifpoints <- as.data.frame(bifpoints)
  indx <- sapply(bifpoints, is.factor)
  bifpoints[indx] <- lapply(bifpoints[indx], function(x) as.numeric(as.character(x)))
  biftype <- as.data.frame(c(pbifmin$biftypes, pbifplus$biftypes))
  bifpoints <- cbind(bifpoints, biftype)
  lst$bifpoints <- as.data.frame(bifpoints)
  if (nrow(bifpoints) > 0) names(lst$bifpoints) <- c(nms, "biftype")
  
  return(lst)
}
##

## --------- FUNCTION: ESS.SCAN.2D() --------- ##
## : -> Function to calculate evolutionary isocline in phase-plane -- ##
ess.scan.2D <- function(modeln = setpars$model.name,
                        bifpar1 = NULL, 
                        bifpar2 = NULL, 
                        ESSpnt1 = NULL, 
                        ESSpnt2 = NULL, 
                        p = setpars$pars.default, 
                        p.names = setpars$pars.names, 
                        par1 = c(0.0,100,0.1),
                        par2 = c(0,1,0.1)){
  if (is.null(bifpar1)) return("Please specify bifurcation parameter 1")
  if(length(which(p.names == bifpar1)) == 0) return('Bifpar1 not found')
  if (is.null(bifpar2)) return("Please specify bifurcation parameter 2")
  if(length(which(p.names == bifpar2)) == 0) return('Bifpar2 not found')
  if (length(par1) != 3) return("Par1 vector is not of length three: minimum parameter value, maximum parameter value, stepsize of the run")
  if (length(par2) != 3) return("Par2 vector is not of length three: minimum parameter value, maximum parameter value, stepsize of the run")
  if (is.null(ESSpnt1)) return("Please specify ESS point for parameter 1")
  if (is.null(ESSpnt2)) print("Only calculating first evolutionary isocline. ESSpnt2 taken from p-array")
  if (length(p) != length(p.names)) return("Number of parameters is unequal to number of parameter names")
  env.dim <- length(setpars$env.names)
  sel <- c(2:(2+env.dim),1)
  
  # Set minpar, maxpar and stepsize
  minpar1 <- par1[1]
  maxpar1 <- par1[2]
  stepsize1 <- par1[3]
  minpar2 <- par2[1]
  maxpar2 <- par2[2]
  stepsize2 <- par2[3]
  
  nms <- Give.Name(type="ESS")
  
  # Store all output in a list
  lst <- list('curvepoints1' = NULL, 'bifpoints1' = NULL, 'curvedesc1' = NULL,
              'curvepoints2' = NULL, 'bifpoints2' = NULL, 'curvedesc2' = NULL)
  
  # Set initial point
  if (!is.data.frame(ESSpnt1)){
    ESSpnt1 <- ldply(ESSpnt1, data.frame)
    ESSpnt1$.id <- sapply(strsplit(ESSpnt1$.id, "="), "[[", 2)
    idx.bftp <- grep(pattern = "biftype", names(ESSpnt1), ignore.case = TRUE)
    ini1 <- ESSpnt1[which(ESSpnt1[,idx.bftp] != "BP #0" & ESSpnt1[,idx.bftp] != "LP"), c(1, sel + 1)]
    parini1 <- c(which(p.names == bifpar2) - 1, minpar2, maxpar2, 0, which(p.names == bifpar1) - 1, minpar1, maxpar1)
  } else {
    idx.bftp <- grep(pattern = "biftype", names(ESSpnt1), ignore.case = TRUE)
    ini1 <- ESSpnt1[which(ESSpnt1[,idx.bftp] != "BP #0" & ESSpnt1[,idx.bftp] != "LP"), sel]
    ini1 <- cbind(bifparone = p[which(p.names == bifpar2)], ini1)
    parini1 <- c(which(p.names == bifpar2) - 1, minpar2, maxpar2, 0, which(p.names == bifpar1) - 1, minpar1, maxpar1)
  }
  
  for (i in 1:nrow(ini1))
  {
    ESSpar1 <- PSPMequi(modeln, 'ESS', as.numeric(ini1[i,]), stepsize1, parini1, p, options = c("popEVO", "0"))
    ESSpar1_back <- PSPMequi(modeln, 'ESS', as.numeric(ini1[i,]), -1*stepsize1, parini1, p, options = c("popEVO", "0"))
    
    # Store all points, sort and select
    pnts1 <- rbind(ESSpar1_back$curvepoints[nrow(ESSpar1_back$curvepoints):1,], ESSpar1$curvepoints)
    pnts1[which(pnts1 == "NAN" | pnts1 == "NaN")] <- 0.0
    pnts1 <- as.data.frame(pnts1)
    indx <- sapply(pnts1, is.factor)
    pnts1[indx] <- lapply(pnts1[indx], function(x) as.numeric(as.character(x)))
    names(pnts1) <- nms
    lst$curvepoints1[[length(lst$curvepoints1) + 1]] <- pnts1
    
    bifpoints <- rbind(ESSpar1_back$bifpoints, ESSpar1$bifpoints)
    bifpoints[which(bifpoints == "NAN" | bifpoints == "NaN")] <- 0.0
    bifpoints <- as.data.frame(bifpoints)
    indx <- sapply(bifpoints, is.factor)
    bifpoints[indx] <- lapply(bifpoints[indx], function(x) as.numeric(as.character(x)))
    biftype <- as.data.frame(c(ESSpar1_back$biftypes, ESSpar1$biftypes))
    bifpoints <- cbind(bifpoints, biftype)
    if(nrow(bifpoints) > 0) names(bifpoints) <- c(nms, "biftype")
    lst$bifpoints1[[length(lst$bifpoints1) + 1]] <- bifpoints
    
    lst$curvedesc1[[length(lst$curvedesc1) + 1]] <- ESSpar1$curvedesc
  }
  
  # Initialize CSS bifpar2
  if (!is.null(ESSpnt2)){
    if (!is.data.frame(ESSpnt2)){
      ESSpnt2 = ldply(ESSpnt2, data.frame)
      ESSpnt2$.id = as.numeric(sapply(strsplit(ESSpnt2$.id, "="), "[[", 2))
      idx.bftp <- grep(pattern = "biftype", names(ESSpnt2), ignore.case = TRUE)
      ini2 = ESSpnt2[which(ESSpnt2[,idx.bftp] != "BP #0" & ESSpnt2[,idx.bftp] != "LP"), c(1, sel + 1)]
      parini2 = c(which(p.names == bifpar1) - 1, minpar1, maxpar1, 0, which(p.names == bifpar2) - 1, minpar2, maxpar2)
    } else {
      idx.bftp <- grep(pattern = "biftype", names(ESSpnt2), ignore.case = TRUE)
      ini2 = ESSpnt2[which(ESSpnt2[,idx.bftp] != "BP #0" & ESSpnt2[,idx.bftp] != "LP"), sel]
      ini2 = cbind(p[which(p.names == bifpar1)], ini2)
      parini2 = c(which(p.names == bifpar1) - 1, minpar1, maxpar1, 0, which(p.names == bifpar2) - 1, minpar2, maxpar2)
    }
    
    for (i in 1:nrow(ini2)){
      ESSpar2 <- PSPMequi(modeln, 'ESS', as.numeric(ini2[i,]), stepsize2, parini2, p, options = c("popEVO", "0"))
      ESSpar2_back <- PSPMequi(modeln, 'ESS', as.numeric(ini2[i,]), -1*stepsize2, parini2, p, options = c("popEVO", "0"))
      
      # Store all points, sort and select
      pnts2 <- rbind(ESSpar2_back$curvepoints[nrow(ESSpar2_back$curvepoints):1,], ESSpar2$curvepoints)
      pnts2[which(pnts2 == "NAN" | pnts2 == "NaN")] <- 0.0
      pnts2 <- as.data.frame(pnts2)
      indx <- sapply(pnts2, is.factor)
      pnts2[indx] <- lapply(pnts2[indx], function(x) as.numeric(as.character(x)))
      names(pnts2) <- nms
      lst$curvepoints2[[length(lst$curvepoints2) + 1]] <- pnts2
      
      bifpoints <- rbind(ESSpar2_back$bifpoints, ESSpar2$bifpoints)
      bifpoints[which(bifpoints == "NAN" | bifpoints == "NaN")] <- 0.0
      bifpoints <- as.data.frame(bifpoints)
      indx <- sapply(bifpoints, is.factor)
      bifpoints[indx] <- lapply(bifpoints[indx], function(x) as.numeric(as.character(x)))
      biftype <- as.data.frame(c(ESSpar2_back$biftypes, ESSpar2$biftypes))
      bifpoints <- cbind(bifpoints, biftype)
      if(nrow(bifpoints) > 0) names(bifpoints) <- c(nms, "biftype")
      lst$bifpoints2[[length(lst$bifpoints2) + 1]] <- bifpoints
      
      lst$curvedesc2[[length(lst$curvedesc2) + 1]] <- ESSpar2$curvedesc
    }
  }
  return(lst)
}
##

## --------- FUNCTION: BP.SCAN.2D() --------- ##
# Function to continue branching point
bp.scan.2D <- function(modeln = setpars$model.name,
                       bifpar1 = NULL,
                       bifpar2 = NULL,
                       env = NULL,
                       BPpnt = c(-1,-1),
                       p = setpars$pars.default,
                       p.names = setpars$pars.names,
                       par1 = c(0.0,50,0.1),
                       par2 = c(0,50,0.1)){
  if (is.null(env)) return("Please specify trivial equilibrium density of the environment in 'env'")
  if (length(env) != length(setpars$env.names)) return('env is not of correct dimension')
  if (is.null(bifpar1)) return("Please specify bifurcation parameter 1")
  if (is.null(bifpar2)) return("Please specify bifurcation parameter 2")
  if(length(which(p.names == bifpar1)) == 0) return('Bifpar1 not found')
  if(length(which(p.names == bifpar2)) == 0) return('Bifpar2 not found')
  if (length(par1) != 3) return("Par1 vector is not of length three: minimum parameter value, maximum parameter value, stepsize of the run")
  if (length(par2) != 3) return("Par2 vector is not of length three: minimum parameter value, maximum parameter value, stepsize of the run")
  if (sum(BPpnt) == -2) return("Please specify BP point")
  if (length(BPpnt) != 2) return("BPpnt is not correctly specified: use c(bifpar1, bifpar1)")
  
  nms <- Give.Name(type = "BP")
  
  # Set minpar, maxpar and stepsize
  minpar1 <- par1[1]
  maxpar1 <- par1[2]
  stepsize1 <- par1[3]
  minpar2 <- par2[1]
  maxpar2 <- par2[2]
  
  #initialize list
  lst <- list(curvepoints = NULL, curvedesc = NULL)
  
  # initialize BP
  ini <- c(BPpnt[1], env, BPpnt[2])
  parini <- as.numeric(c(which(p.names == bifpar1) - 1, minpar1, maxpar1, which(p.names == bifpar2) - 1, minpar2, maxpar2))
  
  # Run forwards and backwards ESS par1
  BPpar <- PSPMequi(modeln, 'BP', as.numeric(ini), stepsize1, parini, p, options = c("popBP", "0"))
  BPpar_back <- PSPMequi(modeln, 'BP', as.numeric(ini), -1*stepsize1, parini, p, options = c("popBP", "0"))
  
  pnts <- rbind(BPpar_back$curvepoints[nrow(BPpar_back$curvepoints):1,], BPpar$curvepoints)
  pnts[which(pnts == "NAN" | pnts == "NaN")] <- 0.0
  pnts <- as.data.frame(pnts)
  indx <- sapply(pnts, is.factor)
  pnts[indx] <- lapply(pnts[indx], function(x) as.numeric(as.character(x)))
  names(pnts) <- nms  
  
  lst$curvepoints <- pnts
  lst$curvedesc <- BPpar$curvedesc
  
  return(lst)
}
##

## --------- FUNCTION: TWODIM.ESS.CONT() --------- ##
# Function to calculate value of two-dimensional ESS point with respect to third parameter
twodim.ess.cont <- function(modeln = setpars$model.name,
                            bifpar = NULL, 
                            parbnds = NULL, 
                            ESSpar1 = NULL, 
                            ESSpar2 = NULL, 
                            ESSpnt = NULL, 
                            p = setpars$pars.default, 
                            p.names = setpars$pars.names, 
                            ESSbnds1 = c(0,2), 
                            ESSbnds2 = c(0,2)){
  if (is.null(bifpar)) return("Please specify bifurcation parameter")
  if (is.null(ESSpar1)) return("Please specify ESS parameter 1")
  if (is.null(ESSpar2)) return("Please specify ESS parameter 2")
  if(length(which(p.names == bifpar)) == 0) return('Bifpar not found')
  if(length(which(p.names == ESSpar1)) == 0) return('ESSpar1 not found')
  if(length(which(p.names == ESSpar2)) == 0) return('ESSpar2 not found')
  if (is.null(ESSpnt)) return("Please specify ESS point for two parameters")
  if (length(parbnds) != 3) return("Parbnds vector is not of length three: minimum parameter value, maximum parameter value, stepsize")
  if (length(ESSbnds1) != 2) return("ESSbnds1 vector is not of length two: minimum parameter value, maximum parameter value")
  if (length(ESSbnds2) != 2) return("ESSbnds2 vector is not of length two: minimum parameter value, maximum parameter value")
  if (length(p) != length(p.names)) return("Number of parameters is unequal to number of parameter names")
  
  env.dim <- length(setpars$env.names)
  sel <- c(2:c(env.dim+3),1)
  
  # Set minpar, maxpar and stepsize
  minpar <- parbnds[1]
  maxpar <- parbnds[2]
  stepsize <- parbnds[3]
  minESSpar1 <- ESSbnds1[1]
  maxESSpar1 <- ESSbnds1[2]
  minESSpar2 <- ESSbnds2[1]
  maxESSpar2 <- ESSbnds2[2]
  
  nms <- Give.Name(type = "ESS", evodims = 2)
  
  # Store all output in a list
  lst <- list('curvepoints' = NULL, 'bifpoints' = NULL, 'curvedesc' = NULL)
  
  # Set initial point
  idx.bftp <- grep(pattern = "biftype", names(ESSpnt), ignore.case = TRUE)
  ini <- ESSpnt[which(ESSpnt[,idx.bftp] != "BP #0" & ESSpnt[,idx.bftp] != "LP"), sel]
  ini <- cbind(bifparone <- p[which(p.names == bifpar)], ini)
  parini <- c(which(p.names == bifpar) - 1, minpar, maxpar, 0, which(p.names == ESSpar1) - 1, minESSpar1, maxESSpar1, 0, which(p.names == ESSpar2) - 1, minESSpar2, maxESSpar2)
  
  ESSpar <- PSPMequi(modeln, 'ESS', as.numeric(ini), stepsize, parini, p, options = c('popEVO', '0'))
  ESSpar_back <- PSPMequi(modeln, 'ESS', as.numeric(ini), -1*stepsize, parini, p, options = c('popEVO', '0'))
  
  # Store all points, sort and select
  pnts <- rbind(ESSpar_back$curvepoints[nrow(ESSpar_back$curvepoints):1,], ESSpar$curvepoints)
  pnts[which(pnts == "NAN" | pnts == "NaN")] <- 0.0
  pnts <- as.data.frame(pnts)
  indx <- sapply(pnts, is.factor)
  pnts[indx] <- lapply(pnts[indx], function(x) as.numeric(as.character(x)))
  names(pnts) <- nms
  lst$curvepoints <- pnts
  
  bifpoints <- rbind(ESSpar_back$bifpoints, ESSpar$bifpoints)
  bifpoints[which(bifpoints == "NAN" | bifpoints == "NaN")] <- 0.0
  bifpoints <- as.data.frame(bifpoints)
  indx <- sapply(bifpoints, is.factor)
  bifpoints[indx] <- lapply(bifpoints[indx], function(x) as.numeric(as.character(x)))
  biftype <- as.data.frame(c(ESSpar_back$biftypes, ESSpar$biftypes))
  bifpoints <- cbind(bifpoints, biftype)
  if(nrow(bifpoints) > 0) names(bifpoints) <- c(nms, "biftype")
  lst$bifpoints <- bifpoints
  
  lst$curvedesc <- ESSpar$curvedesc
  
  return(lst)
}
##
