### FUNCTIONS TO GENERATE MACARTHUR COMMUNITIES ###

if (!require('nleqslv')) install.packages('nleqslv'); library('nleqslv')

# community model =====
macarthur <- function(X, vars) {
  with(vars, {
    # resources
    R <- c(K - (X %*% t(u)))
    # growth
    return(b * (c(t(u * rep(w, each = n.res)) %*% R) - m))
  })
}


# function to calculate species' carrying capacity =====
monocalc <- function(n.spp, vars, extinctions = TRUE) {
  sapply(1:n.spp, function(i) { # monoculture equilibria
    vars.isol <- list('w' = vars$w[i],
                      'b' = vars$b[i],
                      'u' = as.matrix(vars$u[,i], ncol = 1),
                      'm' = vars$m[i],
                      'r' = vars$r,
                      'K' = vars$K)
    if (extinctions == TRUE) {
      equi.isol <- tryCatch(
        uniroot(macarthur, interval = c(0, 1e3), vars = vars.isol)$root, 
        error = function(i) {0})
    } else {
      equi.isol <- tryCatch(
        uniroot(macarthur, interval = c(10, 1e3), vars = vars.isol)$root, 
        error = function(i) {NA})
    }
    return(equi.isol)
  })
}

# function to calculate equilibrium in mixture =====
polycalc <- function(n.spp, equi.isol, vars, extinctions = TRUE) {
  equi.isol[is.na(equi.isol)] <- 0
  equi <- nleqslv(x = equi.isol / 10, fn = macarthur, vars = vars,
                  control = list(maxit = 1e5, allowSingular = TRUE))$x
  
  # allow for species extinctions?
  if (extinctions == TRUE) {
    while (any(equi < 0)) {
      extant <- (1:n.spp)[equi > 0]
      vars.ext <- list('w' = vars$w[extant],
                       'b' = vars$b[extant],
                       'u' = as.matrix(vars$u[,extant], ncol = 1),
                       'm' = vars$m[extant],
                       'r' = vars$r,
                       'K' = vars$K)
      equi <- rep(0, n.spp)
      if (length(extant) == 0) { # if all extinct, return zeros
      } else if (length(extant) == 1) { # if only 1 spp, use uniroot
        equi[extant] <- monocalc(1, vars.ext)
      } else {
        equi[extant] <- nleqslv(
          x = (equi.isol / 10)[extant], fn = macarthur, vars = vars.ext, 
          control = list(maxit = 1e5, allowSingular = TRUE))$x
      }
    }
  }
  return(equi)
}

# function to generate community =====
generate <- function(n.spp, n.res = 16) {
  if (n.spp > 16) {warning('Lots of species, may be slow')}
  repeat { # generating growth parameters
    b <- rep(1, n.spp) 
    m <- runif(n.spp, 1, 3) / 1000 # maintenance / mortality
    w <- runif(n.spp, 1, 3)  # weighting of resources 
    
    # niche-based consumption to ensure coexistence
    pref <- sample(1:n.res, n.spp, replace = FALSE) # define spp niches
    
    u <- sapply(1:n.spp, function(i) {
      u.sd <- runif(1, 0.7, 1.1) * n.res / n.spp
      u.pref <- pref[i]
      dist <- dnorm(1:n.res, u.pref, sd = u.sd)
      dist.2 <- dnorm(1:n.res, u.pref + n.res, sd = u.sd)
      dist.3 <- dnorm(1:n.res, u.pref - n.res, sd = u.sd)
      sapply(1:n.res, function(j) {max(dist[j],dist.2[j],dist.3[j])})
    }, simplify = 'array')
    
    rescale <- array(rep(0.01 / colSums(u), each = n.res), dim = dim(u))
    u <- u * rescale # all species consume equally
    
    r <- K <- rep(0.25, n.res) # resource pars 
    vars <- list('w' = w, 'b' = b, 'u' = u, 'm' = m, 'r' = r, 'K' = K)
    
    # monoculture equilibria
    equi.isol <- monocalc(n.spp, vars, FALSE) # carrying capacities
    
    if (!anyNA(equi.isol)) { 
      equi <- polycalc(n.spp, equi.isol, vars, FALSE) # equilibria in community
      if (all(equi >= 0.01 * equi.isol & equi <= equi.isol)) {break}
    }
  }
  return(list('w' = w, 'b' = b, 'u' = u, 'm' = m, 'r' = r, 'K' = K, 
              'equi' = equi, 'equi.isol' = equi.isol))
}
