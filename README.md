# stressor_richness

`abs_bg.csv`, `pigment_algae_table.csv`, `pigments.csv`, and `sunlight.csv` all contain information needed to generate Stomp phytoplankton communities

`lotkavolterra.R`, `macarthur.R`, and `stomp.R` contain functions used to generate communities of the three community models

`stress.R` contains functions used to generate stressors

`otherfuncs.R` contains functions used in interpretation of results

`demo.R` is the only file that needs user interactions. Set up to perform a factorial design as demonstrated in the paper and plot the results in a similar fashion. Factorial design takes up to a few minutes per repetition. Editing the parameters in the first few lines (reps, n.stress, init.spp, etc.) will change how long the simulations take and how clear the results are. In order to replicate the publication results, use the following parameterisations:

| Parameter | Meaning | Values |
| --- | --- | --- |
| reps | Numer of replicates per parameterisation | 1000 |
| init.spp | Initial species richness | c(4, 8, 16) |
| n.stress | Maximum number of stressors | 20 |
| models | Community model used | c('stomp', 'macarthur', 'lotkavolterra') |
| stress.control | Whether or not stressor action is fixed | c(FALSE, TRUE) |
| d | 1 - total stressor intensity, if stressor action is fixed | c(0.1, 0.5, 0.9) |
| interactions | SD of stressor interactions, generated with normal distribution (mean 0) | c(0, 1) |

Written with R 4.0.3 and requires `nleqslv`, `vegan`, `fda.usc`, and `tidyverse` packages.

If you have any questions, please contact me at mark.holmes@unamur.be
