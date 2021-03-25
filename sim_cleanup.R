# =============================================================================
# script used to clean up stressed community files using clusters
# =============================================================================

file.names <- list.files('Data', pattern = '*summary*', 
                         full.names = TRUE, recursive = TRUE)
all.obj <- lapply(file.names, readRDS)
all.data <- do.call(rbind, all.obj)

saveRDS(all.data, 'Data/all_summary.RData')
