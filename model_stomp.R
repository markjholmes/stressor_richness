#==============================================================================
# set up stomp sim 
#==============================================================================

# import required functions
library('nleqslv')
library('fda.usc')

# growth function
stomp <- function(X, vars) {
  with(vars, {
    # version that models depths explicitly
    #z <- seq(0, zm, by = 0.1) # vector of depth
    #I.z <- I.in * exp(-abs.spec %*% (X %o% z) - I.abs %o% z)
    #gamma <- t(I.z) %*% abs.spec # should be simpsons
    #dXdt <- phi / zm * apply(gamma, 2, int.simpson2, x = z) - mort
    
    # version using jurg's reformulation of the model
    abs <- c(abs.spec %*% X + I.abs) * zm
    gamma <- int.simpson(t(abs.spec / abs * I.in * (1 - exp(-abs))))
    dXdt <- phi * gamma - mort # growth rate
    return(dXdt)
  })
}

zm <- 50#100 

I.in <- read.csv('Data/sunlight_from_graphreader.csv')[,2]
I.in <- I.in / sum(I.in) * 40 #/ 3 #* 10000 # convert to units per m^2 
I.abs.water <- read.csv('Data/water_abs_graphreader.csv')[,2] # water absorption
I.abs.gt <- read.csv('Data/gil_tri_abs_graphreader.csv')[,2] # other absorption
I.abs <- (I.abs.water + I.abs.gt) / 100 # total absorption per cm

# photosynthetic info
pigs.base <- read.csv('Data/pigments.csv')[,-1] # pigment abs from 
pigs <- matrix(nrow = 301, ncol = ncol(pigs.base)) # make length same
for (i in 1:ncol(pigs.base)) {pigs[,i] <- approx(pigs.base[,i], n = 301)$y}
colnames(pigs) <- colnames(pigs.base)
pig.spp <- read.csv('Data/pigment_algae_table.csv')[,-c(1,2)] # phyto taxa pigs

repeat {
  # growth parameters
  spp.id <- sample(1:ncol(pig.spp), n.spp, replace = FALSE) # choose phyto taxa
  
  phi <- runif(n.spp, 1, 2) * 1e6 # similar to stomp 2004
  mort <- 0.015 / 5#00 # jurg version divided by 5000
  abs.spec <- pigs %*% as.matrix(pig.spp[,spp.id] * runif(9 * n.spp, 1, 2))
  for (i in 1:n.spp) {abs.spec[,i] <- abs.spec[,i] / sum(abs.spec[,i])}
  abs.spec <- abs.spec * rep(rnorm(n.spp, 2e-7, 2e-8), each = 301)
  
  vars <- list('phi' = phi, 'mort' = mort, 'abs.spec' = abs.spec, 
               'I.abs' = I.abs, 'I.in' = I.in, 'zm' = zm)

  equi.isol <- sapply(1:n.spp, function(i) {
    vars.isol <- list('phi' = phi[i],
                      'mort' = mort,
                      'abs.spec' = array(abs.spec[,i], 
                                         dim = c(nrow(abs.spec), 1)),
                      'I.abs' = I.abs,
                      'I.in' = I.in,
                      'zm' = zm)
    tryCatch(uniroot(stomp, interval = c(1e6, 1e14), vars = vars.isol)$root, 
             error = function(i) {NA})
  }) 
  
  if (!anyNA(equi.isol)) {
    equi <- nleqslv(x = equi.isol, fn = stomp, vars = vars, 
                    control = list(maxit = 1e5))$x 
    
    if (all(equi >= 0.01 * equi.isol & equi <= equi.isol)) {break}
  }
}
