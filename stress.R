###########################################
### FUNCTION USED TO GENERATE STRESSORS ###
###########################################

if (!require('nleqslv')) install.packages('nleqslv'); library('nleqslv')


stress.func <- function(x, 
                        n.stress = 20, 
                        interactions = 0, 
                        control = TRUE, 
                        d = 0.1) {
  
  n.spp <- length(x) # get n.spp 
  
  # check inputs 
  if (n.spp <= 1 | n.spp %% 1 != 0) { # species richness
    stop('Number of species must be an integer greater than 1')
  }
  if (n.stress <= 1 | !is.numeric(n.stress) | # stressor richness
      is.nan(n.stress) | n.stress %% 1 != 0) {
    stop('Number of stressors must be an integer greater than 1')
  }
  if (!is.numeric(interactions) | interactions < 0) { # stressor interactions
    stop('Give a positive numeric value for the interactions')
  }
  if (!is.logical(control)) {
    stop('State whether stressor action is controlled')
  }
  
  # generate stressors
  stress.eff <- lapply(1:n.stress, function(i) {
    stress <- # generate stressor values
      array(exp(-rgamma(n.spp * i, shape = 1, scale = 0.2)), dim = c(n.spp, i))
    if (control == TRUE) { # controlled stressor action
      d <- array(d)
      init <- prod(stress)
      trans <- array(rep(log(d) / log(init), each = n.spp * i), 
                     dim = dim(stress))
      stress.c <- stress ^ trans
    } else {
      stress.c <- stress # uncontrolled stressor action
    }
    return(stress.c)
  })
  
  # generate stressor interactions
  int.eff <- lapply(1:n.stress, function(i) {
    eta <- array(rnorm(n.spp * i * i, mean = 0, sd = interactions),
                 dim = c(n.spp, i, i)) # generate interaction values
    for (j in 1:i) {eta[,j,j] <- 0} # remove self-interactions
    eta <- sapply(1:i, function(l) { # account for  stressor intensity
      sapply(1:i, function(k) {
        sapply(1:n.spp, function(j) {
          f1 <- log(stress.eff[[i]][j,k])
          f2 <- log(stress.eff[[i]][j,l])
          return(exp(f1 * f2 * eta[j,k,l]))
        }, simplify = 'array')
      }, simplify = 'array')
    }, simplify = 'array')
    return(apply(eta, 1, prod)) # return species-specific values
  })
  
  # get total effect
  t.eff <- lapply(1:n.stress, function(i) {
    return(apply(stress.eff[[i]], 1, prod) * int.eff[[i]])
  })
  
  # apply stress to traits
  x.out <- sapply(t.eff, function(i) {x * i}, simplify = 'array')
  
  # calculate stressor diversity
  div <- sapply(1:n.stress, function(i) {
    stressor.diversity(recalc.tot(stress.eff[[i]], int.eff[[i]]))
  })
  
  return(list('x.out' = x.out, 'stress' = 1:n.stress, 'div' = div,
              'interactions' = interactions, 'control' = control))
}