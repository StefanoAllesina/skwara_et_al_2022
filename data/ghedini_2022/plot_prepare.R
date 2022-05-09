library(tidyverse)
dt <- read_csv("Biovolume_by species.csv")
# plotting density cells per microliter
dt <- dt %>% 
  mutate(days_old = ifelse(days_old == "target", -1, days_old)) %>% 
  mutate(days_old = as.numeric(days_old)) %>% rename(experiment = sample_code) %>% 
  select(days_old, experiment, species_counted, cells_ul) %>% 
  rename(time = days_old, species = species_counted, density = cells_ul) %>% 
    filter(time >0) %>% 
  # rename species
  mutate(species = ifelse(species == "Amphidinium", "A", species)) %>% 
  mutate(species = ifelse(species == "Dunaliella", "D", species)) %>% 
  mutate(species = ifelse(species == "Synechococcus", "S", species)) %>% 
  mutate(species = ifelse(species == "Tetraselmis", "T", species)) %>% 
  mutate(species = ifelse(species == "Tisochrysis", "Ti", species))
  
  
pl <- ggplot(dt, aes(x = time, y = density, colour = species, group = interaction(species, experiment))) + 
  geom_rect(xmin = 8 - 0.5, xmax = 8 + 0.5, ymin = -Inf, ymax = Inf, colour = NA,fill = 'lightgrey', alpha = 0.1) +
  geom_point() + geom_line() + 
  scale_y_log10("cells/\u00B5l") + facet_wrap(~str_match(experiment, "[a-zA-Z]+"), scales = "free") + 
  theme_bw() + xlab("days")  + theme(legend.position = "bottom") + scale_color_brewer(palette = "Set1")
ggsave(pl, filename = "ghedini_2022_time_series.pdf", width = 8, height = 6)
# select time point day 8
dt2 <- dt %>% filter(time == 8) %>% pivot_wider(names_from = species, values_from = density, values_fill = 0) %>% select(-time, -experiment)
write_csv(dt2, file = "../ghedini_2022.csv")
