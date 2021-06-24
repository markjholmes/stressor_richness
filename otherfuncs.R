### OTHER FUNCTIONS USED TO CALCULATE INFORMATION ABOUT COMMUNITIES ###

if (!require('vegan')) install.packages('vegan'); library('vegan')

# coefficient of variation function
CV <- function(x) {
  (sd(x) / mean(x)) * 100
}

# covariance function
E <- function(i, j) {
  mean(i * j) - (mean(i) * mean(j))
}

# function for recalculating total stress using interactive and non-interactive
recalc.tot <- function(s.eff, i.eff) {
  if(is.vector(s.eff)) {
    s.eff <- matrix(s.eff, nrow = length(a), ncol = 1)
  }
  i.eff <- array(i.eff, dim = dim(s.eff))
  return(s.eff * i.eff)
}

# function to compute bray-curtis similarity index 
bray.curtis.sim <- function(equis, init.equis) {
  equi.df <- rbind(init.equis, equis)
  return(1 - unname(vegdist(equi.df, na.rm = TRUE)))
}

# selection and complementarity =====

# function to calculate complementarity and selection effects
comp.sel.func <- function(Mi, Yoi) {
  # skip bad cases 
  if(any(Mi == 0)) {
    return(list('deltaY' = NA, 
                'complementarity' = NA,
                'selection' = NA,
                'n.spp' = NA))} 
  
  N <- length(Yoi) # number of species
  Yo <- sum(Yoi) # observed total yield in mix
  RYei <- rep(1 / N, N) # expected relative yield in mix
  RYoi <- Yoi / Mi # observed relative yield in mix
  Yei <- RYei * Mi # expected yield in mix
  Ye <- sum(Yei) # expected total yield in mix
  deltaY <- Yo - Ye # deviation from total expected yield in mix
  deltaRYi <- RYoi - RYei # deviation from expected relative yield in mix
  complementarity <- N * mean(deltaRYi) * mean(Mi) 
  selection <- N * (E(deltaRYi, Mi))  
  complementarity <- complementarity / deltaY # divide by delta Y
  selection <- selection / deltaY
  return(list('complementarity' = complementarity,
              'selection' = selection
  ))
}