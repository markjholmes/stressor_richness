#######################################################################
### OTHER FUNCTIONS USED TO CALCULATE INFORMATION ABOUT COMMUNITIES ###
#######################################################################

if (!require('vegan')) install.packages('vegan'); library('vegan')

# =============================================================================
# stressor diversity functions

stressor.diversity <- function(a) {
  if(is.vector(a)) {
    a <- matrix(a, nrow = length(a), ncol = 1)
  }
  rescale <- 1 / prod(dim(a)[2], 2)
  d <- (abs(det(t(a) %*% a)) ^ rescale) / norm(a, type = '2')
  return(d) 
}

# recalculating total stress
recalc.tot <- function(s.eff, i.eff) {
  if(is.vector(s.eff)) {
    s.eff <- matrix(s.eff, nrow = length(a), ncol = 1)
  }
  i.eff <- array(i.eff, dim = dim(s.eff))
  return(s.eff * i.eff)
}

# =============================================================================
# b-c similarity index

bray.curtis.sim <- function(equis, init.equis) {
  equi.df <- rbind(init.equis, equis)
  return(1 - unname(vegdist(equi.df, na.rm = TRUE)))
}

# =============================================================================
# selection and complementarity

# function to calculate complementarity and selection effects
comp.sel.func <- function(Mi, Yoi) {
  N <- length(Yoi) # number of species
  Yo <- sum(Yoi) # observed total yield in mix
  RYei <- Mi / sum(Mi) #rep(1 / N, N) # expected relative yield in mix
  RYoi <- Yoi / Mi # observed relative yield in mix
  Yei <- RYei * Mi # expected yield in mix
  Ye <- sum(Yei) # expected total yield in mix
  deltaY <- Yo - Ye # deviation from total expected yield in mix
  deltaRYi <- RYoi - RYei # deviation from expected relative yield in mix
  complementarity <- N * mean(deltaRYi) * mean(Mi) / deltaY
  selection <- N * cov(deltaRYi, Mi) / deltaY
  return(list('deltaY' = deltaY, 
              'complementarity' = complementarity,
              'selection' = selection, 
              'n.spp' = N))
}