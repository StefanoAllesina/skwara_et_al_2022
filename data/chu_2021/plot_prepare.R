library(tidyverse)
dt <- read_tsv("raw_data.csv")
dt <- dt %>% select(Transfer, ID, no.P, no.S, no.O, no.A, no.V)
dt <- dt %>% rename(time = Transfer, experiment = ID, 
                    P = no.P, S = no.S, O = no.O, A = no.A, V = no.V)
dt <- dt %>% pivot_longer(names_to = "species", values_to = "density", cols = -c("time", "experiment")) %>% 
  filter(!is.na(density))
# plotting colony count
pl <- ggplot(dt, aes(x = time, y = density, colour = species, group = interaction(species, experiment))) + 
  geom_rect(
    data = dt %>% filter(str_match(experiment, "[a-zA-Z]+") %in% c("P", "PA", "PO", "POA", "PS", "PSA")),
    xmin = 3 - 0.25, xmax = 3 + 0.25, ymin = -Inf, ymax = Inf, colour = NA,fill = 'lightgrey', alpha = 0.1) +
  geom_point() + geom_line() + geom_hline(yintercept = 0, linetype = 2, colour = "black")+ 
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 50), "colonies",
                     breaks = c(0,10,100,250,1000,2000)) + facet_wrap(~str_match(experiment, "[a-zA-Z]+"), scales = "free") + 
  theme_bw() + xlab("transfer")  + theme(legend.position = "bottom") + scale_color_brewer(palette = "Set1")
ggsave(pl, filename = "chu_2021_time_series.pdf", width = 8, height = 6)
# select experiments excluding V and those yielding extinctions
dt2 <- dt %>% filter(time == 3) %>% filter(str_match(experiment, "[a-zA-Z]+") %in% c("P", "PA", "PO", "POA", "PS", "PSA")) %>% 
  pivot_wider(names_from = species, values_from = density, values_fill = 0) %>% select(-time, -experiment)
write_csv(dt2, file = "../chu_2021.csv")
