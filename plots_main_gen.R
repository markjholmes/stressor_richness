# =============================================================================
# initial setup

# load packages and clear workspace
library('tidyverse')
library('lemon')
library('ggpubr')
library('scales')

# set figure scale
size.scale <- 1

# import data and fix
all.data <- readRDS('Data/all_summary.RData') %>% # import
  mutate(control = replace_na(control, 'Variable intensity')) %>% # rename NAs
  filter(interactions == 0, i.rich == 4, control != '0.5', control != '0.9') %>%
  mutate(control = recode(control, `0.1` = 'Fixed intensity'),
         model = recode(model, 'lotka-volterra' = 'Lotka-Volterra',
                        'macarthur' = 'MacArthur',
                        'stomp' = 'Stomp')) %>% 
  mutate(control = fct_relevel(control, 'Variable intensity')) %>%
  mutate(s.intensity = 1 - t.eff, 
         d.rich = s.rich / i.rich, 
         d.pop = s.pop / i.pop) %>% 
  dplyr::select(-t.eff, -i.rich, -i.pop, -s.pop) %>%
  pivot_longer(cols = i.sel:d.pop) %>%
  mutate(name = recode(name, 
                       'bc.sim' = 'Compositional resistance',
                       's.intensity' = 'Stressor intensity',
                       'd.rich' = 'Species persistence',
                       'd.pop' = 'Ecosystem functioning',
                       's.sel' = 'Selection effect',
                       's.comp' = 'Complementarity effect'))

# =============================================================================
# Summarise by stressor richness

# plot 1 frames
data.1 <- all.data %>%
  dplyr::filter(name == 'Stressor intensity')

means.1 <- data.1 %>% 
  group_by(n.stress, control, name) %>% 
  summarise(value = mean(value),
            s.cv = mean(s.cv))

labs.1 <- data.frame(
  label = letters[1:2],
  model = factor('Stomp'),
  control = factor(unique(means.1$control)),
  name = factor(rep(unique(means.1$name), each = 2)),
  s.cv = 0, value = 0,
  x = -Inf, y = Inf)

# plot 2 frames
data.2 <- all.data %>%
  filter(name == 'Ecosystem functioning' | 
           name == 'Species persistence' | 
           name == 'Compositional resistance' | 
           name == 'Selection effect' | 
           name == 'Complementarity effect' |
           name == 'i.sel' | name == 'i.comp') %>%
  pivot_wider(names_from = name, values_from = value) %>%
  filter(is.finite(`Selection effect`), 
         is.finite(`Complementarity effect`)) %>%
  filter(#`Complementarity effect` < 
         #  quantile(`Complementarity effect`, 0.9, na.rm = T),
         #`Complementarity effect` > 
         #  quantile(`Complementarity effect`, 0.1, na.rm = T),
         #`Selection effect` <
         #  quantile(`Selection effect`, 0.9, na.rm = T),
         #`Selection effect` > 
         #  quantile(`Selection effect`, 0.1, na.rm = T),
         `Complementarity effect` > 0, `Selection effect` > 0,
         i.comp > 0, i.sel > 0) %>%
  group_by(n.stress, control, model) %>%
  pivot_longer(cols = i.sel:`Ecosystem functioning`) %>%
  filter(!is.na(value)) 

# initial values

i.sel <- data.2[data.2$name == 'i.sel',] %>% 
  group_by(model) %>% 
  summarise(value = mean(value)) %>% 
  add_column(s.cv = 0, name = 'Selection effect') %>% 
  uncount(20, .id = 'n.stress') %>%
  mutate(name = factor(name, levels = c('Ecosystem functioning', 
                                        'Species persistence', 
                                        'Compositional resistance',
                                        'Selection effect',
                                        'Complementarity effect')))

i.comp <- data.2[data.2$name == 'i.comp',] %>% 
  group_by(model) %>% 
  summarise(value = mean(value)) %>% 
  add_column(s.cv = 0, name = 'Complementarity effect') %>% 
  uncount(20, .id = 'n.stress') %>%
  mutate(name = factor(name, levels = c('Ecosystem functioning', 
                                        'Species persistence', 
                                        'Compositional resistance',
                                        'Selection effect',
                                        'Complementarity effect')))

data.2 <- data.2 %>%
  filter(name != 'i.sel', name != 'i.comp')

i.other <- i.comp %>%
  slice(rep(1:n(), each = 3)) %>%
  mutate(value = 1, name = NA) %>%
  mutate(name = rep(c('Ecosystem functioning', 
                      'Species persistence', 
                      'Compositional resistance'), 60)) %>%
  mutate(name = factor(name, levels = c('Ecosystem functioning', 
                                        'Species persistence', 
                                        'Compositional resistance',
                                        'Selection effect',
                                        'Complementarity effect')))

means.2 <- data.2 %>%
  filter(name != 'i.sel', name != 'i.comp') %>% 
  group_by(n.stress, control, model, name) %>% 
  summarise(value = mean(value),
            s.cv = mean(s.cv)) %>%
  mutate(name = factor(name, levels = c('Ecosystem functioning', 
                                        'Species persistence', 
                                        'Compositional resistance',
                                        'Selection effect',
                                        'Complementarity effect')))

data.2 <- data.2 %>% 
  mutate(name = factor(name, levels = c('Ecosystem functioning', 
                                        'Species persistence', 
                                        'Compositional resistance',
                                        'Selection effect',
                                        'Complementarity effect')))

labs.2 <- data.frame(
  label = letters[1:10],
  model = factor('Stomp'),
  control = factor(unique(means.2$control)),
  name = factor(rep(levels(means.2$name), each = 2)),
  s.cv = 0, value = 0,
  x = -Inf, y = Inf)

# =============================================================================
# plotting theme
theme.main <- theme_classic() + 
  theme(legend.position = 'bottom', 
        legend.box = 'horizontal',
        legend.title = element_text(size = 10, colour = 'black'), 
        text = element_text(size = 10, colour = 'black'),
        strip.background = element_blank(), 
        strip.placement = 'outside',
        strip.text = element_text(size = 10, colour = 'black'), 
        axis.title = element_text(size = 10, colour = 'black'),
        axis.ticks = element_line(colour = 'black'), 
        axis.line = element_line(),
        plot.title = element_text(size = 10, hjust = 0.5),
        plot.margin = unit(c(1, 1, 1, 1), 'pt'),
        panel.background = element_blank())

# =============================================================================
# plots

# 1
fig.1 <- ggplot(means.1, 
                aes(x = n.stress, y = s.cv, fill = 100*value)) +
  facet_grid(. ~ control, scales = 'free_y', switch = 'y') + 
  geom_point(data = data.1,
             alpha = 0.01, size = 0.5, stroke = 0, pch = 21) +
  geom_point(size = 1.5, pch = 21) + 
  geom_text(data = labs.1, 
            mapping = aes(x = x, y = y, label = label),
            hjust = -1.25, vjust = 1.25, col = 'black') +
  theme.main +
  scale_x_continuous(breaks = c(1, 5, 10, 15, 20)) +
  scale_fill_viridis_c(
    option = 'D',
    guide = guide_colorbar(title.position = 'left', direction = 'horizontal'),
    limits = range(means.1$value * 100), 
    oob = squish) +
  labs(x = 'Stressor richness',
       fill = 'Stressor intensity (%)',
       y = 'SCV') +
  annotate("segment", size = 1.1, x = -Inf, xend = Inf, y = -Inf, yend = -Inf) +
  annotate("segment", size = 1.1, x = -Inf, xend = -Inf, y = -Inf, yend = Inf)

ggsave('Figures/fig_Srich_SCV.jpeg', fig.1, device = 'jpeg', 
       width = 80, height = 70, scale = size.scale, units = 'mm')

# 2
fig.2 <- ggplot(means.2, 
                aes(x = n.stress, lty = model, pch = model, 
                    y = value, fill = s.cv)) +
  facet_grid(name ~ control, switch = 'y') +
  geom_line(data = i.other) +
  geom_line(data = i.sel) +
  geom_line(data = i.comp) +
  geom_point(data = data.2,
             alpha = 0.01, size = 0.5, stroke = 0, 
             position = position_dodge(1)) +
  geom_point(size = 1.5, position = position_dodge(1)) + 
  geom_text(data = labs.2, 
            mapping = aes(x = x, y = y, label = label),
            hjust = -1.25, vjust = 1, col = 'black') +
  theme.main +
  coord_cartesian(xlim = c(1, 20), ylim = c(0, 1.1)) +
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  scale_x_continuous(breaks = c(1, 5, 10, 15, 20)) +
  scale_fill_viridis_c(
    option = 'C',
    guide = guide_colorbar(title.position = 'top', direction = 'horizontal'),
    limits = range(means.2$s.cv), 
    oob = squish) +
  scale_linetype_manual(
    name = 'Model',
    labels = unique(means.2$model),
    values = c(2, 3, 4),
    guide = guide_legend(direction = 'vertical', reverse = TRUE)) +
  scale_shape_manual(
    name = 'Model',
    labels = unique(means.2$model),
    values = c(21, 22, 23),
    guide = guide_legend(direction = 'vertical', reverse = TRUE)) +
  labs(x = 'Stressor richness',
       y = NULL,
       fill = 'SCV') +
  theme(strip.text = element_text(size = 8.5, colour = 'black')) + 
  annotate("segment", size = 1.1, x = -Inf, xend = Inf, y = -Inf, yend = -Inf) +
  annotate("segment", size = 1.1, x = -Inf, xend = -Inf, y = -Inf, yend = Inf)

ggsave('Figures/fig_factorial.jpeg', fig.2, device = 'jpeg', 
       width = 80, height = 220, scale = size.scale, units = 'mm')
