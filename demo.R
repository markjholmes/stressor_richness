### DEMONSTRATION  ###

# source functions needed
if (!require('tidyverse')) install.packages('tidyverse'); library('tidyverse')
source('stress.R')
source('otherfuncs.R')

# factorial design settings to generate communities and stressors =====

reps <- 1 # how many repeats at each factor combination
models <- c('stomp', 'macarthur', 'lotkavolterra') # community models
stress.control <- c(FALSE, TRUE) # stressor control
interactions <- c(0, 1) # interaction presence
n.stress <- 20 # max stressor richness
init.spp <- c(4, 8) # initial species richnesses 
d <- 0.1 # 1 - total stressor intensity. ignored for uncontrolled stress action

# combine factors
combinations <- expand.grid('rep' = 1:reps,
                            'model' = models, 
                            'control' = stress.control, 
                            'interactions' = interactions, 
                            'n.spp' = init.spp,
                            'n.stress' = n.stress,
                            'd' = d)

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
  d <- combinations$d[comb]
  
  # generate community
  source(paste0(model, '.R'))
  community <- generate(n.spp) 
  
  # stress first parameter in community 
  stressors <- stress.func(
    community[[1]], 
    n.stress = n.stress, 
    interactions = interactions, 
    control = control, 
    d = d) 
  
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
                  'd' = d,
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

factorial.summary <- factorial.out %>%
  group_by(rep, control, d, interactions, model, spp.rich, n.stress) %>%
  summarise(d.pop = sum(equi) / sum(i.equi),
            d.rich = sum(equi > 0) / sum(i.equi > 0),
            sel = comp.sel.func(equi.isol, equi)$selection,
            comp = comp.sel.func(equi.isol, equi)$complementarity,
            scv = scv[1]) %>%
  mutate(spp.rich = as.factor(spp.rich))

# plotting results =====
# all results are facetted by interaction and stressor control

# stressor richness against stressor diversity
ggplot(factorial.summary, aes(x = n.stress, y = scv, 
                              lty = spp.rich, pch = spp.rich)) +
  facet_grid(interactions ~ control, labeller = 'label_both') +
  stat_summary(fun = mean, geom = 'point') +
  stat_summary(fun = mean, geom = 'line')

# effects of stressor richness
ggplot(factorial.summary, aes(x = n.stress, y = d.pop, col = model,
                              lty = spp.rich, pch = spp.rich)) +
  facet_grid(interactions ~ control, labeller = 'label_both') +
  stat_summary(fun = mean, geom = 'point', alpha = 1/reps) +
  stat_summary(fun = mean, geom = 'line') 

ggplot(factorial.summary, aes(x = n.stress, y = d.rich, col = model,
                              lty = spp.rich, pch = spp.rich)) +
  facet_grid(interactions ~ control, labeller = 'label_both') +
  stat_summary(fun = mean, geom = 'point', alpha = 1/reps) +
  stat_summary(fun = mean, geom = 'line')

ggplot(factorial.summary, aes(x = n.stress, y = bc.sim, col = model,
                              lty = spp.rich, pch = spp.rich)) +
  facet_grid(interactions ~ control, labeller = 'label_both') +
  stat_summary(fun = mean, geom = 'point', alpha = 1/reps) +
  stat_summary(fun = mean, geom = 'line') 

ggplot(factorial.summary, aes(x = n.stress, y = sel, col = model,
                              lty = spp.rich, pch = spp.rich)) +
  facet_grid(interactions ~ control, labeller = 'label_both') +
  stat_summary(fun = mean, geom = 'point', alpha = 1/reps) +
  stat_summary(fun = mean, geom = 'line') 

ggplot(factorial.summary, aes(x = n.stress, y = comp, col = model,
                              lty = spp.rich, pch = spp.rich)) +
  facet_grid(interactions ~ control, labeller = 'label_both') +
  stat_summary(fun = mean, geom = 'point', alpha = 1/reps) +
  stat_summary(fun = mean, geom = 'line') 

# effects of stressor coefficient of variation
# here geom_smooth is used for simplicity, rather than dividing stressor 
# diversity into bins as shown in paper, using cut or stat_summary_bin
ggplot(factorial.summary, aes(x = scv, y = d.pop, col = model,
                              lty = spp.rich, pch = spp.rich)) +
  facet_grid(interactions ~ control, labeller = 'label_both') +
  geom_point(alpha = 1/reps) + 
  scale_x_log10() +
  geom_smooth(se = F)

ggplot(factorial.summary, aes(x = scv, y = d.rich, col = model,
                              lty = spp.rich, pch = spp.rich)) +
  facet_grid(interactions ~ control, labeller = 'label_both') +
  geom_point(alpha = 1/reps) +
  scale_x_log10() +
  geom_smooth(se = F)

ggplot(factorial.summary, aes(x = scv, y = bc.sim, col = model,
                              lty = spp.rich, pch = spp.rich)) +
  facet_grid(interactions ~ control, labeller = 'label_both') +
  geom_point(alpha = 1/reps) + 
  scale_x_log10() +
  geom_smooth(se = F)

ggplot(factorial.summary, aes(x = scv, y = sel, col = model,
                              lty = spp.rich, pch = spp.rich)) +
  facet_grid(interactions ~ control, labeller = 'label_both') +
  geom_point(alpha = 1/reps) + 
  scale_x_log10() +
  geom_smooth(se = F)

ggplot(factorial.summary, aes(x = scv, y = comp, col = model,
                              lty = spp.rich, pch = spp.rich)) +
  facet_grid(interactions ~ control, labeller = 'label_both') +
  geom_point(alpha = 1/reps) + 
  scale_x_log10() +
  geom_smooth(se = F)
