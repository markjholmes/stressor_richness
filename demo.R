### DEMONSTRATION  ###

# source functions needed
if (!require('tidyverse')) install.packages('tidyverse'); library('tidyverse')
source('stress.R')
source('otherfuncs.R')

# factorial design settings to generate communities and stressors =====

reps <- 5 # how many repeats at each factor combination, low is fast, high gives
          # better results. remove geom_point() if you make this high
models <- c('stomp', 'macarthur', 'lotkavolterra') # community models
stress.control <- c(FALSE, TRUE) # stressor control
interactions <- c(0, 1) # interaction presence
n.stress <- 20 # max stressor richness
init.spp <- c(4, 8, 16) # initial species richness 

# combine factors
combinations <- expand.grid('rep' = 1:reps,
                            'model' = models, 
                            'control' = stress.control, 
                            'interactions' = interactions, 
                            'n.spp' = init.spp,
                            'n.stress' = n.stress)

# remove impossible simulations due to stomp limitations on sp richness =====

combinations[combinations$model == 'stomp' & combinations$n.spp > 4,] <- NA
combinations <- na.omit(combinations)

# perform factorial design (this may take several minutes per rep) =====
factorial.out <- lapply(1:nrow(combinations), function(comb) {
  
  # get factor details 
  rep <- combinations$rep[comb]
  model <- combinations$model[comb]
  control <- combinations$control[comb]
  interactions <- combinations$interactions[comb]
  n.stress <- combinations$n.stress[comb]
  n.spp <- combinations$n.spp[comb]
  n.stress <- n.stress
  
  # generate community
  source(paste0(model, '.R'))
  community <- generate(n.spp) 
  
  # stress first parameter in community 
  stressors <- stress.func(
    community[[1]], 
    n.stress = n.stress, 
    interactions = interactions, 
    control = control) 
  
  # calculate new stressed equilibrium
  stressed.comm <- lapply(1:n.stress, function(s) {
    s.vars <- community
    s.vars[[1]] <- stressors$x.out[,s]
    equi.isol <- monocalc(n.spp, s.vars)
    equi <- polycalc(n.spp, equi.isol, s.vars, TRUE)
    return(list('equi.isol' = equi.isol, 'equi' = equi))
  }) 
  
  # summarise community changes
  df.all <- cbind('rep' = rep,
                  'control' = control,
                  'interactions' = interactions,
                  'model' = model,
                  'spp.rich' = n.spp,
                  'n.stress' = rep(1:n.stress, each = n.spp),
                  'scv' = rep(stressors$scv, each = n.spp),
                  do.call(rbind.data.frame, stressed.comm),
                  'i.equi.isol' = community$equi.isol,
                  'i.equi' = community$equi)
  print(paste(comb, 'of', nrow(combinations), 'complete'))
  return(df.all)
})

# cleanup output
factorial.out <- do.call(rbind, factorial.out)

# summarising data =====

factorial.summary <- by(factorial.out, factorial.out[,1:6], function(x) {
  
  rep <- x$rep[1]
  control <- x$control[1]
  interactions <- x$interactions[1]
  model <- as.character(x$model[1])
  spp.rich <- x$spp.rich[1]
  n.stress <- x$n.stress[1]
  scv <- x$scv[1]
  
  d.pop <- sum(x$equi) / sum(x$i.equi)
  d.rich <- sum(x$equi > 0) / sum(x$i.equi > 0)
  bc.sim <- bray.curtis.sim(x$equi, x$i.equi) # may produce warnings
  sel <- comp.sel.func(x$equi.isol, x$equi)$selection # may produce warnings
  comp <-  comp.sel.func(x$equi.isol, x$equi)$complementarity
  
  return(list('rep' = rep, 'control' = control, 'interactions' = interactions, 
              'model' = model, 'spp.rich' = as.factor(spp.rich), 
              'n.stress' = n.stress, 'scv' = scv,
              'd.pop' = d.pop, 'd.rich' = d.rich, 'bc.sim' = bc.sim,
              'sel' = sel, 'comp' = comp))
})

factorial.summary <- do.call(rbind.data.frame, factorial.summary)

# plotting results =====

# stressor richness against stressor diversity
ggplot(factorial.summary, aes(x = n.stress, y = scv, 
                              lty = spp.rich, pch = spp.rich)) +
  facet_grid(cols = vars(control), rows = vars(interactions)) +
  stat_summary(fun = mean, geom = 'point') +
  stat_summary(fun = mean, geom = 'line') +
  scale_y_log10() 

# effect of stressor richness
ggplot(factorial.summary, aes(x = n.stress, y = d.pop, col = model,
                              lty = spp.rich, pch = spp.rich)) +
  facet_grid(cols = vars(control), rows = vars(interactions)) +
  stat_summary(fun = mean, geom = 'point', alpha = 1/reps) +
  stat_summary(fun = mean, geom = 'line') 

ggplot(factorial.summary, aes(x = n.stress, y = d.rich, col = model,
                              lty = spp.rich, pch = spp.rich)) +
  facet_grid(cols = vars(control), rows = vars(interactions)) +
  stat_summary(fun = mean, geom = 'point', alpha = 1/reps) +
  stat_summary(fun = mean, geom = 'line')

ggplot(factorial.summary, aes(x = n.stress, y = bc.sim, col = model,
                              lty = spp.rich, pch = spp.rich)) +
  facet_grid(cols = vars(control), rows = vars(interactions)) +
  stat_summary(fun = mean, geom = 'point', alpha = 1/reps) +
  stat_summary(fun = mean, geom = 'line') 

ggplot(factorial.summary, aes(x = n.stress, y = sel, col = model,
                              lty = spp.rich, pch = spp.rich)) +
  facet_grid(cols = vars(control), rows = vars(interactions)) +
  stat_summary(fun = mean, geom = 'point', alpha = 1/reps) +
  stat_summary(fun = mean, geom = 'line') 

ggplot(factorial.summary, aes(x = n.stress, y = comp, col = model,
                              lty = spp.rich, pch = spp.rich)) +
  facet_grid(cols = vars(control), rows = vars(interactions)) +
  stat_summary(fun = mean, geom = 'point', alpha = 1/reps) +
  stat_summary(fun = mean, geom = 'line') 

# effect of stressor diversity
# here geom_smooth is used for simplicity, rather than dividing stressor 
# diversity into bins as shown in paper, which can be achieved using cut()
ggplot(factorial.summary, aes(x = scv, y = d.pop, col = model,
                              lty = spp.rich, pch = spp.rich)) +
  facet_grid(cols = vars(control), rows = vars(interactions)) +
  geom_point(alpha = 1/reps) + 
  scale_x_log10() +
  geom_smooth(se = F)

ggplot(factorial.summary, aes(x = scv, y = d.rich, col = model,
                              lty = spp.rich, pch = spp.rich)) +
  facet_grid(cols = vars(control), rows = vars(interactions)) +
  geom_point(alpha = 1/reps) +
  scale_x_log10() +
  geom_smooth(se = F)

ggplot(factorial.summary, aes(x = scv, y = bc.sim, col = model,
                              lty = spp.rich, pch = spp.rich)) +
  facet_grid(cols = vars(control), rows = vars(interactions)) +
  geom_point(alpha = 1/reps) + 
  scale_x_log10() +
  geom_smooth(se = F)

ggplot(factorial.summary, aes(x = scv, y = sel, col = model,
                              lty = spp.rich, pch = spp.rich)) +
  facet_grid(cols = vars(control), rows = vars(interactions)) +
  geom_point(alpha = 1/reps) + 
  scale_x_log10() +
  geom_smooth(se = F)

ggplot(factorial.summary, aes(x = scv, y = comp, col = model,
                              lty = spp.rich, pch = spp.rich)) +
  facet_grid(cols = vars(control), rows = vars(interactions)) +
  geom_point(alpha = 1/reps) + 
  scale_x_log10() +
  geom_smooth(se = F)
