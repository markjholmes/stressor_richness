# =============================================================================
# unified function for subjecting communities to stress
# =============================================================================

# import required functions 
library('nleqslv')
library('vegan')

# b-c similarity index
bray.curtis.sim <- function(equis, init.equis) {
  equi.df <- rbind(init.equis, equis)
  return(c(1 - vegdist(equi.df, na.rm = TRUE)))
}

# covariance function
E <- function(i, j) {
  mean(i * j) - (mean(i) * mean(j))
}

# function to calculate complementarity and selection effects
comp.sel.func <- function(Mi, Yoi) {
  if(any(Mi == 0)) {
    return(list('deltaY' = NA, 
                'complementarity' = NA,
                'selection' = NA,
                'n.spp' = NA))} # TRY THIS?
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
  return(list(#'deltaY' = deltaY, 
              'complementarity' = complementarity,
              'selection' = selection
              ))
}

stress.func <- function(n.spp = 4, n.stress = 20, interactions = 0, n.comm = 1, 
                        model = c('stomp', 'macarthur', 'lotka-volterra'),
                        control = NA, parallel.id = NULL) {
  
  # ensure method is correct
  match.arg(model, c('stomp', 'macarthur', 'lotka-volterra'))
  
  # ensure args are okay
  if (n.spp %% 1 != 0 | n.stress %% 1 != 0| !is.numeric(interactions) | 
      n.comm %% 1 != 0 | interactions < 0 | interactions > 1) {
    stop('Input incorrect')
  }
  
  # ensure output folder exists and create new subdir
  model.names <- list('stomp' = 'Stomp', 
                      'macarthur' = 'MacArthur', 
                      'lotka-volterra' = 'Lotka-Volterra')
  modelname <- model.names[model]
  
  sim.pars <- list(spp = n.spp, stress = n.stress)
  newdir <- paste0('Data/', modelname, '/int_', interactions, '/',
                   paste0(names(sim.pars), sim.pars, collapse = '_'),
                   '/control_', control)
  
  dir.create(file.path(newdir), recursive = TRUE)
  
  # get growth change function
  source('func_growth_change.R')
  
  # Stomp version
  if (model == 'stomp') {
    if (n.spp > 4) {
      stop('Stomp model can have at maximum 4 species')
    }
    
    # generate communities 
    comms <- lapply(1:n.comm, function(id) {
      source('model_stomp.R', local = TRUE)
      
      # get pars
      comm <- list('phi' = phi,
                   'mort' = mort,
                   'abs.spec' = abs.spec,
                   'I.abs' = I.abs,
                   'I.in' = I.in,
                   'zm' = zm,
                   'equi' = equi,
                   'equi.isol' = equi.isol)
      
      rm(phi, mort, abs.spec, equi, equi.isol)
      
      # generate stress
      stress.info <- growth.change(x = comm$phi, 
                                   n.stress = n.stress,
                                   interactions = interactions, 
                                   control = control)
      
      # get new equi
      stressed.communities <- lapply(1:n.stress, function(i) {
        # base pars
        s.phi <- c(stress.info$x.out[,i])
        s.vars <- list('phi' = s.phi,
                       'mort' = comm$mort,
                       'abs.spec' = comm$abs.spec,
                       'I.abs' = comm$I.abs,
                       'I.in' = comm$I.in,
                       'zm' = comm$zm)
        
        # compute monoculture equi 
        equi.isol <- sapply(1:n.spp, function(j) {
          vars.isol <- list('phi' = s.phi[j],
                            'mort' = comm$mort,
                            'abs.spec' = array(comm$abs.spec[,j], 
                                               dim = c(nrow(comm$abs.spec),1)),
                            'I.abs' = comm$I.abs,
                            'I.in' = comm$I.in,
                            'zm' = comm$zm)
          tryCatch(uniroot(stomp, 
                           interval = c(0, max(comm$equi.isol)),
                           vars = vars.isol)$root, 
                   error = function(i) {0})
        })
        
        # compute equi 
        x0 <- comm$equi * s.vars$phi / comm$phi
        equi <- nleqslv(x = x0, fn = stomp, vars = s.vars, 
                        control = list(maxit = 1e5, allowSingular = TRUE))$x
        
        # try combinations with extinctions
        while (any(equi < 0)) {
          extant <- (1:n.spp)[equi > 0]
          vars.ext <- list('phi' = s.phi[extant],
                           'mort' = comm$mort,
                           'abs.spec' = array(comm$abs.spec[,extant], 
                                              dim = c(nrow(comm$abs.spec),
                                                      length(extant))),
                           'I.abs' = comm$I.abs,
                           'I.in' = comm$I.in,
                           'zm' = comm$zm)
          equi <- rep(0, n.spp)
          
          if (length(extant) == 0) {
          } else if (length(extant) == 1) {
            equi[extant] <- tryCatch(uniroot(
              f = stomp, 
              interval = c(0, max(comm$equi.isol)),
              vars = vars.ext)$root,
              error = function(i) {0})
          } else {
            equi[extant] <- nleqslv(x = x0[extant], fn = stomp, 
                                    vars = vars.ext, 
                                    control = list(maxit = 1e5, 
                                                   allowSingular = TRUE))$x
          }
        }
        return(list('model' = model, 'interactions' = interactions, 
                    'control' = control, 'n.spp' = n.spp, 
                    'id' = id, 'n.stress' = i, 
                    'equi' = equi, 'equi.isol' = equi.isol))
      })
      
      print(paste(model, 'model,', n.spp, 'species,', 
                  interactions, 'interactions,', 
                  control, 'control,', id, 'of', n.comm, 
                  'complete at', Sys.time()))
      
      return(list('id' = id, 
                  'comm' = comm, 
                  'stress.info' = stress.info,
                  'stressed.communities' = stressed.communities))
    })
  } else if (model == 'macarthur') { # MacArthur version
    
    # generate communities 
    n.res <- 16
    comms <- lapply(1:n.comm, function(id) { 
      source('model_macarthur.R', local = TRUE)
      
      # list pars
      comm <- list('b' = b,
                   'u' = u,
                   'w' = w,
                   'm' = m,
                   'r' = r,
                   'K' = K,
                   'equi' = equi,
                   'equi.isol' = equi.isol)
      
      rm(b, u, w, m, r, K, equi, equi.isol)
      
      # generate stress
      stress.info <- growth.change(x = comm$w,
                                   n.stress = n.stress,
                                   interactions = interactions, 
                                   control = control)
      
      # get new equi
      stressed.communities <- lapply(1:n.stress, function(i) {
        # base pars
        s.w <- stress.info$x.out[,i]
        s.vars <- list('b' = comm$b,
                       'u' = comm$u,
                       'w' = s.w,
                       'm' = comm$m,
                       'r' = comm$r,
                       'K' = comm$K)
        
        # compute monoculture equi 
        equi.isol <- sapply(1:n.spp, function(j) {
          vars.isol <- list('b' = comm$b[j],
                            'u' = as.matrix(comm$u[,j], ncol = 1),
                            'w' = s.w[j],
                            'm' = comm$m[j],
                            'r' = comm$r,
                            'K' = comm$K)
          tryCatch(uniroot(macarthur, 
                           interval = c(0, max(comm$equi.isol)),
                           vars = vars.isol)$root, 
                   error = function(i) {0})
        })
        
        # compute equi
        x0 <- comm$equi * s.vars$w / comm$w
        equi <- nleqslv(x = x0, fn = macarthur, vars = s.vars, 
                        control = list(maxit = 1e5))$x
        
        # check for extinctions and calculate equi without them
        while (any(equi < 0)) {
          extant <- (1:n.spp)[equi > 0]
          vars.ext <- list('b' = comm$b[extant],
                           'u' = as.matrix(comm$u[,extant], 
                                           ncol = length(extant)),
                           'w' = s.w[extant],
                           'm' = comm$m[extant],
                           'r' = comm$r,
                           'K' = comm$K)
          equi <- rep(0, n.spp)
          
          if (length(extant) == 0) {
          } else if (length(extant) == 1) {
            equi[extant] <- tryCatch(
              uniroot(macarthur, 
                      interval = c(0, max(comm$equi.isol)),
                      vars = vars.ext)$root, 
              error = function(i) {0})
          } else {
            equi[extant] <- nleqslv(x = x0[extant], 
                                    fn = macarthur,
                                    vars = vars.ext, 
                                    control = list(maxit = 1e5, 
                                                   allowSingular = TRUE))$x
          }
        }
        return(list('model' = model, 'interactions' = interactions, 
                    'control' = control, 'n.spp' = n.spp, 
                    'id' = id, 'n.stress' = i, 
                    'equi' = equi, 'equi.isol' = equi.isol))
      })
      
      print(paste(model, 'model,', n.spp, 'species,', 
                  interactions, 'interactions,', 
                  control, 'control,', id, 'of', n.comm, 
                  'complete at', Sys.time()))

      return(list('id' = id, 
                  'comm' = comm, 
                  'stress.info' = stress.info,
                  'stressed.communities' = stressed.communities))
    })
  } else if (model == 'lotka-volterra') { # Lotka-Volterra version
    
    # generate communities     
    comms <- lapply(1:n.comm, function(id) {
      source('model_lotka-volterra.R', local = TRUE)
      
      # get pars
      comm <- list('id' = id, 'n.spp' = n.spp, 
                   'mu' = mu,
                   'alphas' = alphas,
                   'equi' = equi,
                   'equi.isol' = equi.isol)
      
      rm(mu, alphas, equi, equi.isol)
      
      # generate stress
      stress.info <- growth.change(x = comm$mu, 
                                   n.stress = n.stress,
                                   interactions = interactions, 
                                   control = control)
      
      # get new equi
      stressed.communities <- lapply(1:n.stress, function(i) {
        # base pars
        s.mu <- c(stress.info$x.out[,i])
        s.vars <- list('mu' = s.mu,
                       'alphas' = comm$alphas)
        
        # compute monoculture equi        
        equi.isol <- sapply(1:n.spp, function(j) {
          vars.isol <- list('mu' = s.mu[j],
                            'alphas' = comm$alphas[j,j])
          tryCatch(uniroot(lotkavolterra, 
                           interval = c(0, max(comm$equi.isol)), 
                           vars = vars.isol)$root, 
                   error = function(i) {0})
        })
        
        # compute equi 
        x0 <- comm$equi * s.vars$mu / comm$mu
        equi <- nleqslv(x = x0, fn = lotkavolterra, vars = s.vars, 
                        control = list(maxit = 1e5))$x
        
        # try combinations with extinctions
        while (any(equi < 0)) {
          extant <- (1:n.spp)[equi > 0]
          vars.ext <- list('mu' = s.mu[extant],
                           'alphas' = comm$alphas[extant,extant])
          equi <- rep(0, n.spp)
          
          if (length(extant) == 0) {
          } else if (length(extant) == 1) {
            equi[extant] <- tryCatch(uniroot(
              lotkavolterra, 
              interval = c(0, max(comm$equi.isol)), 
              vars = vars.ext)$root, 
              error = function(i) {0})
            
          } else {
            equi[extant] <- nleqslv(x = x0[extant], fn = lotkavolterra,
                                    vars = vars.ext, 
                                    control = list(maxit = 1e5, 
                                                   allowSingular = TRUE))$x
          }
        }

        return(list('model' = model, 'interactions' = interactions, 
                    'control' = control, 'n.spp' = n.spp, 
                    'id' = id, 'n.stress' = i, 
                    'equi' = equi, 'equi.isol' = equi.isol))
      })
      
      print(paste(model, 'model,', n.spp, 'species,', 
                  interactions, 'interactions,', 
                  control, 'control,', id, 'of', n.comm, 
            'complete at', Sys.time()))
      
      return(list('id' = id, 
                  'comm' = comm, 
                  'stress.info' = stress.info,
                  'stressed.communities' = stressed.communities))
    })
  }
  
  # summary info
  # simulation settings
  model <- sapply(1:n.comm, function(i) {
    sapply(1:n.stress, function(j) {
      comms[[i]]$stressed.communities[[j]]$model
    })
  })
  interactions <- sapply(1:n.comm, function(i) {
    sapply(1:n.stress, function(j) {
      comms[[i]]$stressed.communities[[j]]$interactions
    })
  })
  control <- sapply(1:n.comm, function(i) {
    sapply(1:n.stress, function(j) {
      comms[[i]]$stressed.communities[[j]]$control
    })
  })
  
  # initial community parameters
  rich <- sapply(1:n.comm, function(i) {
    sapply(1:n.stress, function(j) {
      sum(comms[[i]]$comm$equi > 0)
    })
  })
  pop <- sapply(1:n.comm, function(i) {
    sapply(1:n.stress, function(j) {
      sum(comms[[i]]$comm$equi)
    })
  })
  comp <- sapply(1:n.comm, function(i) {
    sapply(1:n.stress, function(j) {
      comp.sel.func(Mi = comms[[i]]$comm$equi.isol,
                    Yoi = comms[[i]]$comm$equi)$comp
    })
  })
  sel <- sapply(1:n.comm, function(i) {
    sapply(1:n.stress, function(j) {
      comp.sel.func(Mi = comms[[i]]$comm$equi.isol,
                    Yoi = comms[[i]]$comm$equi)$sel
    })
  })
  #deltaY <- sapply(1:n.comm, function(i) {
  #  sapply(1:n.stress, function(j) {
  #    comp.sel.func(Mi = comms[[i]]$comm$equi.isol,
  #                  Yoi = comms[[i]]$comm$equi)$deltaY
  #  })
  #})
  # stressed parameters
  stress <- sapply(1:n.comm, function(i) {
    as.numeric(comms[[i]]$stress.info$stress)
  })
  t.eff <- sapply(1:n.comm, function(i) {
    sapply(1:n.stress, function(j) {
      prod(comms[[i]]$stress.info$t.eff[[j]])
    })
  })
  s.cv <- sapply(1:n.comm, function(i) {
    comms[[i]]$stress.info$s.cv
  })
  s.rich <- sapply(1:n.comm, function(i) {
    sapply(1:n.stress, function(j) {
      sum(comms[[i]]$stressed.communities[[j]]$equi > 0)
    })
  })
  s.pop <- sapply(1:n.comm, function(i) {
    sapply(1:n.stress, function(j) {
      sum(comms[[i]]$stressed.communities[[j]]$equi)
    })
  })
  s.comp <- sapply(1:n.comm, function(i) {
    sapply(1:n.stress, function(j) {
      comp.sel.func(Mi = comms[[i]]$stressed.communities[[j]]$equi.isol,
                    Yoi = comms[[i]]$stressed.communities[[j]]$equi)$comp
    })
  })
  s.sel <- sapply(1:n.comm, function(i) {
    sapply(1:n.stress, function(j) {
      comp.sel.func(Mi = comms[[i]]$stressed.communities[[j]]$equi.isol,
                    Yoi = comms[[i]]$stressed.communities[[j]]$equi)$sel
    })
  })
  #s.deltaY <- sapply(1:n.comm, function(i) {
  #  sapply(1:n.stress, function(j) {
  #    comp.sel.func(Mi = comms[[i]]$stressed.communities[[j]]$equi.isol,
  #                  Yoi = comms[[i]]$stressed.communities[[j]]$equi)$deltaY
  #  })
  #})
  bc <- sapply(1:n.comm, function(i) {
    sapply(1:n.stress, function(j) {
      bray.curtis.sim(comms[[i]]$stressed.communities[[j]]$equi,
                      comms[[i]]$comm$equi)
    })
  })

  cov.st.sp <- sapply(1:n.comm, function(i) {
    sapply(1:n.stress, function(j) {
      cov(comms[[i]]$comm$equi,
          comms[[i]]$stress.info$t.eff[[j]])
    })
  })
  df.all <- data.frame('model' = c(model),
                       'interactions' = c(interactions),
                       'control' = c(control),
                       'n.stress' = c(stress),
                       't.eff' = c(t.eff),
                       's.cv' = c(s.cv),
                       'i.rich' = c(rich),
                       'i.pop' = c(pop),
                       'i.sel' = c(sel),
                       'i.comp' = c(comp),
                       's.rich' = c(s.rich),
                       's.pop' = c(s.pop),
                       's.sel' = c(s.sel),
                       's.comp' = c(s.comp),
                       'bc.sim' = c(bc),
                       'cov.st.sp' = c(cov.st.sp))
  
  # save communities, if being performed in parallel, save parallel id
  if (!is.null(parallel.id)) {
    saveRDS(comms, paste0(newdir, '/all_data_', parallel.id, '.RData'))
    saveRDS(df.all, paste0(newdir, '/summary_', parallel.id, '.RData'))
  } else {
    saveRDS(comms, paste0(newdir, '/all_data.RData'))
    saveRDS(df.all, paste0(newdir, '/summary.RData'))
  }
}
