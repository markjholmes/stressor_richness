### FUNCTIONS TO GENERATE STOMP PHYTOPLANKTON COMMUNITIES ###

if (!require('nleqslv')) install.packages('nleqslv'); library('nleqslv')
if (!require('fda.usc')) install.packages('fda.usc'); library('fda.usc')

# community model  =====
stomp <- function(X, vars) {
  with(vars, {
    abs <- c(abs.spec %*% X + I.abs) * zm # total light absorption
    gamma <- int.simpson(t(abs.spec / abs * I.in * (1 - exp(-abs)))) # spp abs
    dXdt <- phi * gamma - mort # growth rate
    return(dXdt) 
  })
}

# function to calculate species' carrying capacity  =====
monocalc <- function(n.spp, vars, extinctions = TRUE) {
  sapply(1:n.spp, function(i) { 
    vars.isol <- list('phi' = vars$phi[i],
                      'mort' = vars$mort,
                      'abs.spec' = matrix(vars$abs.spec[,i], ncol = 1),
                      'I.abs' = vars$I.abs,
                      'I.in' = vars$I.in,
                      'zm' = vars$zm)
    # if extinctions permitted, allow lower equilibrium values
    if (extinctions == TRUE) {
      equi.isol <- tryCatch(
        uniroot(stomp, interval = c(0, 1e10), vars = vars.isol)$root, 
        error = function(i) {0})
    } else {
      equi.isol <- tryCatch(
        uniroot(stomp, interval = c(1e6, 1e10), vars = vars.isol)$root, 
        error = function(i) {NA})
    }
    return(equi.isol)
  }) 
}

# function to calculate equilibrium in mixture  =====
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
                       'I.in' = vars$I.in,
                       'zm' = vars$zm)
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

# function to generate community  =====
generate <- function(n.spp) {
  if (n.spp > 4) {stop('Stomp communities can only have 4 species')}
  
  I.in <- read.csv('sunlight.csv')[,2] # incident light
  I.in <- I.in / sum(I.in) * 40 # rescale
  I.abs <- read.csv('abs_bg.csv')[,2] # background light absorption
  zm <- 50
  
  # photosynthetic info
  pigs <- as.matrix(read.csv('pigments.csv')[,-1])
  pig.spp <- read.csv('pigment_algae_table.csv')[,-c(1,2)]
  
  repeat { # generating growth parameters
    # select species
    spp.id <- sample(1:ncol(pig.spp), n.spp, replace = FALSE) 
    
    # growth parameters
    phi <- runif(n.spp, 1, 3) * 1e6 # generate photosynthetic efficiency
    mort <- 0.003 # mortality/loss
    # add randomisation to absorption
    abs.spec <- pigs %*% as.matrix(pig.spp[,spp.id] * runif(9 * n.spp, 1, 2))
    # rescale absorption
    for (i in 1:n.spp) { 
      abs.spec[,i] <- abs.spec[,i] / sum(abs.spec[,i])
    }
    # smaller randomisation
    abs.spec <- abs.spec * rep(rnorm(n.spp, 2e-7, 2e-8), each = 301)
    vars <- list('phi' = phi, 'mort' = mort, 'abs.spec' = abs.spec, 
                 'I.abs' = I.abs, 'I.in' = I.in, 'zm' = zm)
    equi.isol <- monocalc(n.spp, vars, FALSE) # carrying capacities
    
    if (!anyNA(equi.isol)) { 
      # community equilibria 
      equi <- polycalc(n.spp, equi.isol, vars, FALSE) 
      
      if (all(equi >= 0.01 * equi.isol & equi <= equi.isol)) {break}
    }
  }
  return(list('phi' = phi, 'abs.spec' = abs.spec, 'mort' = mort, 
              'spp.id' = spp.id, 'I.abs' = I.abs, 'I.in' = I.in, 'zm' = zm,
              'equi' = equi, 'equi.isol' = equi.isol))
}
