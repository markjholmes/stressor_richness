#==============================================================================
# set up lotka volterra sim 
#==============================================================================

# import required functions 
library('nleqslv')

# growth function
lotkavolterra <- function(X, vars) {
  with(vars, {
    return(mu - (alphas %*% X))
  })
}

# get working sets of pars
repeat {
  #mu <- abs(rnorm(n.spp, 1, 0.1)) 
  mu <- runif(n.spp, 0.8, 1.2)
  alphas <- matrix(rnorm(n.spp ^ 2, 1, 0.1), 
                   nrow = n.spp, ncol = n.spp) / (400)
  diag(alphas) <- rnorm(n.spp, 2, 0.2) / (300)
  
  vars <- list('mu' = mu, 'alphas' = alphas)
  
  equi.isol <- sapply(1:n.spp, function(i) {
    vars.isol <- list('mu' = mu[i],
                 'alphas' = alphas[i,i])
    tryCatch(uniroot(lotkavolterra, interval = c(10, 1e3), 
                     vars = vars.isol)$root, error = function(i) {NA})
  })
  
  if (!anyNA(equi.isol)) {
    equi <- nleqslv(x = equi.isol / 10, fn = lotkavolterra, vars = vars,
                    control = list(maxit = 1e5, allowSingular = TRUE))$x
    
    if (all(equi >= 0.01 * equi.isol & equi <= equi.isol)) {break}
  }
}
