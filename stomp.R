#############################################################
### FUNCTIONS TO GENERATE STOMP PHYTOPLANKTON COMMUNITIES ###
#############################################################

if (!require('nleqslv')) install.packages('nleqslv'); library('nleqslv')

# community model
stomp <- function(X, vars) {
  zm <- 7.7 # depth
  phi <- vars$phi # photosynthetic efficiency
  mort <- vars$mort # mortality
  abs.spec <- vars$abs.spec # absorption spectra
  I.abs <- vars$I.abs # absorption of light by water
  I.in <- vars$I.in # incoming light
  I.zm <- I.in * exp(-(abs.spec %*% X + I.abs) * zm) # light throughout depth
  gamma.zm <- t(abs.spec) %*% I.zm # absorption of light by species
  return(phi / zm * gamma.zm - mort)
}

# calculate species' carrying capacity
monocalc <- function(n.spp, vars, extinctions = TRUE) {
  sapply(1:n.spp, function(i) { 
    vars.isol <- list('phi' = vars$phi[i],
                      'mort' = vars$mort,
                      'abs.spec' = array(vars$abs.spec[,i],
                                         dim = c(nrow(vars$abs.spec), 1)),
                      'I.abs' = vars$I.abs,
                      'I.in' = vars$I.in)
    if (extinctions == TRUE) {
      equi.isol <- tryCatch(
        uniroot(stomp, interval = c(0, 1e11), vars = vars.isol)$root, 
        error = function(i) {0})
    } else {
      equi.isol <- tryCatch(
        uniroot(stomp, interval = c(1e9, 1e11), vars = vars.isol)$root, 
        error = function(i) {NA})
    }
    return(equi.isol)
  }) 
}

# calculate equilibrium in mixture
polycalc <- function(n.spp, equi.isol, vars, extinctions = TRUE) {
  equi.isol[is.na(equi.isol)] <- 0
  equi <- nleqslv(x = equi.isol / 10, fn = stomp, vars = vars,
                  control = list(maxit = 1e5, allowSingular = TRUE))$x
  
  # allow for species extinctions?
  if (extinctions == TRUE) {
    while (any(equi < 0)) {
      extant <- (1:n.spp)[equi > 0]
      vars.ext <- list('phi' = vars$phi[extant],
                       'mort' = vars$mort,
                       'abs.spec' = array(vars$abs.spec[,extant],
                                          dim = c(nrow(vars$abs.spec),
                                                  length(extant))),
                       'I.abs' = vars$I.abs,
                       'I.in' = vars$I.in)
      equi <- rep(0, n.spp)
      if (length(extant) == 0) { # if all extinct, return zeros
      } else if (length(extant) == 1) { # if only 1 spp, use uniroot
        equi[extant] <- monocalc(1, vars.ext)
      } else {
        equi[extant] <- nleqslv(
          x = (equi.isol / 10)[extant], fn = stomp, vars = vars.ext, 
          control = list(maxit = 1e5, allowSingular = TRUE))$x
      }
    }
  }
  return(equi)
}

# generate community
generate <- function(n.spp) {
  if (n.spp > 4) {stop('Stomp communities can only have 4 species')}
  
  I.in <- read.csv('sunlight.csv')[,2] # incident light
  I.in <- I.in / sum(I.in) * 4e5 
  I.abs <- read.csv('abs_bg.csv')[,2] # background light absorption
  
  # photosynthetic info
  pigs <- as.matrix(read.csv('pigments.csv')[,-1])
  pig.spp <- read.csv('pigment_algae_table.csv')[,-c(1,2)]
  
  repeat { # generating growth parameters
    spp.id <- sample(1:ncol(pig.spp), n.spp, replace = FALSE)
    phi <- runif(n.spp, 1, 3) * 1e6
    mort <- 0.014 / 3600 
    abs.spec <- pigs %*% (as.matrix(pig.spp[,spp.id]) * 
                            runif(9 * n.spp, 0.5, 1.5))
    for (i in 1:n.spp) {abs.spec[,i] <- abs.spec[,i] / sum(abs.spec[,i])}
    abs.spec <- abs.spec * 2e-7
    vars <- list('phi' = phi, 'mort' = mort, 'abs.spec' = abs.spec, 
                 'I.abs' = I.abs, 'I.in' = I.in)
    equi.isol <- monocalc(n.spp, vars, FALSE) # carrying capacities
    
    if (!anyNA(equi.isol)) { 
      equi <- polycalc(n.spp, equi.isol, vars, FALSE) # equilibria in community
      if (all(equi >= 0.01 * equi.isol & equi <= equi.isol)) {break}
    }
  }
  return(list('phi' = phi, 'abs.spec' = abs.spec, 'mort' = mort, 
              'spp.id' = spp.id, 'I.abs' = I.abs, 'I.in' = I.in,
              'equi' = equi, 'equi.isol' = equi.isol))
}
