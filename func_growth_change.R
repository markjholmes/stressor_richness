#==============================================================================
# function for producing growth rate changes 
#==============================================================================

# coefficient of variation function
CV <- function(x) {
  (sd(x) / mean(x)) * 100
}

growth.change <- function(x, n.stress = 20, interactions = 0, control = NA) {
  
  # INPUTS:
  # x: a vector containing the growth parameters to be stressed 
  
  # n.stress: an integer for the maximum number of stressors to be simulated,
  #           all stressor richnesses up to this number will be simulated
  
  # interactions: positive number representing the standard deviation of the
  #               stressor interactions that will be generated
  
  # control: indicating whether stressor action is controlled, if not NA it is 
  #          community-level remaining functioning of the growth parameter, as
  #          a number between 0 and 1 
  
  # OUTPUT LIST:
  # x.out: matrix with rows for number of species and columns for number of
  #        stressors, showing the stressed values for the growth parameters
  #        at all stressor richnesses
  
  # stress.eff: list with length n.stress, with each element containing the
  #             species-stressor matrix at that list levels stressor richness
  
  # OUTPUT LIST:
  # x.out: matrix with rows for number of species and columns for number of
  #        stressors, showing the stressed values for the growth parameters
  #        at all stressor richnesses
  
  # stress.eff: list with length n.stress, with each element containing the
  #             species-stressor matrix at that list levels stressor richness
  
  # int.eff: 
  # t.eff:
  # stress: 
  # control:
  # div:
  
  # determine species richness
  n.spp <- length(x)
  
  stress.eff <- lapply(1:n.stress, function(i) {
    # generate species-specific stressor action using gamma distribution, 
    # producing a matrix with rows representing species and columns representing
    # stressors
    # using gamma
    #stress <- exp(-rgamma(n.spp * i, shape = 1, scale = 0.1)) # change the scale parameter to alter mean stress intensity
    # using beta
    stress <- rbeta(n.spp * i, shape1 = 6.5, shape2 = 0.25)
    stress <- array(stress, dim = c(n.spp, i))
    # control stressor action
    if (!is.na(control)) {
      control <- array(control)
      init <- prod(stress)
      trans <- array(rep(log(control) / log(init), each = n.spp * i), 
                     dim = dim(stress))
      stress.c <- stress ^ trans
    } else {
      stress.c <- stress
    }
    return(stress.c)
  })
  
  # generate stressor interactions using normal distribution, producting an 
  # array of three dimensions: stressor 1, stressor 2, and species. Stressor 
  # interaction strength depends on the initial intensity of the stressors 
  int.eff <- lapply(1:n.stress, function(i) {
    # generate initial interactions
    eta <- array(rnorm(n.spp * i * i, mean = 0, sd = interactions),
                 dim = c(n.spp, i, i))
    # remove self-interactions
    for (j in 1:i) {eta[,j,j] <- 0}
    # factor in stressor intensity
    eta <- sapply(1:i, function(l) {
      sapply(1:i, function(k) {
        sapply(1:n.spp, function(j) {
          f1 <- log(stress.eff[[i]][j,k]) # stressor 1 effect
          f2 <- log(stress.eff[[i]][j,l]) # stressor 2 effect
          return(exp(f1 * f2 * eta[j,k,l])) # combine stress effects and eta
        }, simplify = 'array')
      }, simplify = 'array')
    }, simplify = 'array') 
    eta <- array(eta, dim = c(n.spp, i, i)) # ensure correct shape
    return(apply(eta, 1, prod)) # return species-specific values
  })
  
  # get total stressor effect on growth parameter, including interactions
  t.eff <- lapply(1:n.stress, function(i) {
      return(apply(stress.eff[[i]], 1, prod) * int.eff[[i]])
    })
  
  # apply stress to growth parameter x
  x.out <- sapply(t.eff, function(i) {x * i}, simplify = 'array')

  # calculate coefficient of variation
  s.cv <- sapply(t.eff, CV)
  
  return(list('x.out' = x.out, 'stress.eff' = stress.eff, 'int.eff' = int.eff, 
              't.eff' = t.eff, 'stress' = as.list(1:n.stress), 
              'control' = control, 's.cv' = s.cv))
}
