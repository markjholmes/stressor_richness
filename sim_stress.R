# =============================================================================
# script used to generate and stress communities using clusters
# =============================================================================

source('func_stress.R')

# get array id
parallel.id <- Sys.getenv('SLURM_ARRAY_TASK_ID')
if (parallel.id == '') {parallel.id <- NULL}

# parameters
n.stress <- 20
n.comm <- 10

# generate stressed.communities
for (interactions in c(0, 1)) {
  for (control in c(NA, 0.1, 0.5, 0.9)) {
    for (model in c('lotka-volterra', 'macarthur')) {
      for (n.spp in c(4, 8, 16)) {
        stress.func(n.spp = n.spp, 
                    n.stress = n.stress, 
                    interactions = interactions, 
                    n.comm = n.comm, 
                    model = model, 
                    parallel.id = parallel.id,
                    control = control)
      }
    }
    model <- 'stomp'
    n.spp <- 4
    stress.func(n.spp = n.spp, 
                n.stress = n.stress, 
                interactions = interactions, 
                n.comm = n.comm, 
                model = model, 
                parallel.id = parallel.id,
                control = control)
  }      
}

