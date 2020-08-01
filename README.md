# stressor_richness

`abs_bg.csv`, `pigment_algae_table.csv`, `pigments.csv`, and `sunlight.csv` all contain information needed to generate Stomp phytoplankton communities

`lotkavolterra.R`, `macarthur.R`, and `stomp.R` contain functions used to generate communities of the three community models

`stress.R` contains functions used to generate stressors

`otherfuncs.R` contains functions used in interpretation of results

`demo.R` is the only file that needs user-end interactions. Set up to perform a factorial design as demonstrated in the paper and plot the results in a similar fashion. Factorial design takes up to a few minutes per repetition.

Written with R 3.6.0 and requires three packages: "nleqslv", "vegan", and "ggplot2".

If you have any questions, please contact me at mark.holmes@unamur.be
