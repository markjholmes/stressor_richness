########################################################
### FUNCTIONS TO GENERATE LOTKA-VOLTERRA COMMUNITIES ###
########################################################

if (!require('nleqslv')) install.packages('nleqslv'); library('nleqslv')

# community model
lotkavolterra <- function(X, vars) {
  mu <- vars$mu # exponential growth rate
  alphas <- vars$alpha # limitation
  return(mu - (alphas %*% X))
}

#  calculate species' carrying capacity
monocalc <- function(n.spp, vars, extinctions = TRUE) {
  sapply(1:n.spp, function(i) { # monoculture equilibria
    vars.isol <- list('mu' = vars$mu[i],
                      'alphas' = vars$alphas[i,i])
    if (extinctions == TRUE) {
      equi.isol <- tryCatch(
        uniroot(lotkavolterra, interval = c(0, 1e3), vars = vars.isol)$root, 
        error = function(i) {0})
    } else {
      equi.isol <- tryCatch(
        uniroot(lotkavolterra, interval = c(10, 1e3), vars = vars.isol)$root, 
        error = function(i) {NA})
    }
    return(equi.isol)
  })
}

# calculate equilibrium in mixture
polycalc <- function(n.spp, equi.isol, vars, extinctions = TRUE) {
  equi.isol[is.na(equi.isol)] <- 0
  equi <- nleqslv(x = equi.isol / 10, fn = lotkavolterra, vars = vars,
                  control = list(maxit = 1e5, allowSingular = TRUE))$x
  
  # allow for species extinctions?
  if (extinctions == TRUE) {
    while (any(equi < 0)) {
      extant <- (1:n.spp)[equi > 0]
      vars.ext <- list('mu' = vars$mu[extant],
                       'alphas' = as.matrix(vars$alphas[extant,extant]))
      equi <- rep(0, n.spp)
      if (length(extant) == 0) { # if all extinct, return zeros
      } else if (length(extant) == 1) { # if only 1 spp, use uniroot
        equi[extant] <- monocalc(1, vars.ext)
      } else {
        equi[extant] <- nleqslv(
          x = (equi.isol / 10)[extant], fn = lotkavolterra, vars = vars.ext, 
          control = list(maxit = 1e5, allowSingular = TRUE))$x
      }
    }
  }
  return(equi)
}

# generate community
generate <- function(n.spp) {
  if (n.spp > 16) {warning('Lots of species, may be slow')}
  
  repeat { # generating growth parameters
    mu <- abs(rnorm(n.spp, 1, 0.1)) 
    alphas <- matrix(rnorm(n.spp ^ 2, 1, 0.1), 
                     nrow = n.spp, ncol = n.spp) / (400)
    diag(alphas) <- rnorm(n.spp, 2, 0.2) / (300)
    vars <- list('mu' = mu, 'alphas' = alphas)
    equi.isol <- monocalc(n.spp, vars, FALSE) # carrying capacities
    if (!anyNA(equi.isol)) { # equilibria in community
      equi <- polycalc(n.spp, equi.isol, vars, FALSE)
      if (all(equi >= 0.01 * equi.isol & equi <= equi.isol)) {break}
    }
  }
  return(list('mu' = mu, 'alphas' = alphas, 
              'equi' = equi, 'equi.isol' = equi.isol))
}