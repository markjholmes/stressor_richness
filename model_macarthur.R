# =============================================================================
# set up MacArthur sim
# =============================================================================

# import required functions 
library('nleqslv')

macarthur <- function(X, vars) {
  with(vars, {
    # resources
    R <- c(K - (X %*% t(u)))
    # growth
    return(b * (c(t(u * rep(w, each = n.res)) %*% R) - m))
  })
}

# get working sets of pars
repeat {
  # consumer parameters
  b <- rep(1, n.spp) #rnorm(n.spp, 1, 0.1)
  m <- runif(n.spp, 1, 3) / 1000
  w <- runif(n.spp, 1, 3) # weighting of resources 
  
  # niche-based consumption to ensure coexistence
  pref <- sample(1:n.res, n.spp, replace = FALSE)
  #pref <- seq(1, n.res, length.out = n.spp + 1)[1:n.spp]
  
  u <- sapply(1:n.spp, function(i) {
    u.sd <- runif(1, 0.7, 1.1) * n.res / n.spp
    u.pref <- pref[i]
    dist <- dnorm(1:n.res, u.pref, sd = u.sd)
    dist.2 <- dnorm(1:n.res, u.pref + n.res, sd = u.sd)
    dist.3 <- dnorm(1:n.res, u.pref - n.res, sd = u.sd)
    sapply(1:n.res, function(j) {max(dist[j],dist.2[j],dist.3[j])})
  }, simplify = 'array')
  
  rescale <- array(rep(0.01 / colSums(u), each = n.res), dim = dim(u))
  u <- u * rescale
  
  # resource par
  r <- K <- rep(0.25, n.res) #rowSums(u) * 100 #rnorm(n.res, 100, 30)
  
  w <- runif(n.spp, 10, 20) #* 10 * c(m / r %*% u)
  
  vars <- list('w' = w, 'b' = b, 'u' = u, 'm' = m, 'r' = r, 'K' = K)
  
  equi.isol <- sapply(1:n.spp, function(i) {
    vars.isol <- list('b' = b[i],
                      'u' = as.matrix(u[,i], ncol = 1),
                      'w' = w[i],
                      'm' = m[i],
                      'r' = r,
                      'K' = K)
    tryCatch(uniroot(macarthur, interval = c(10, 1e3), vars = vars.isol)$root, 
             error = function(i) {NA})
  })
  
  if (!anyNA(equi.isol)) {
    equi <- nleqslv(x = equi.isol / 10, fn = macarthur, vars = vars,
                    control = list(maxit = 1e5, allowSingular = TRUE))$x
    
    if (all(equi >= 0.01 * equi.isol & equi <= equi.isol)) {break}
  }
}
