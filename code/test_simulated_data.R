# Generate data
n <- 4
rep <- 25
seed <- 2
#source("../data/simulated_data/generate_gamma_distributed_data.R")
source("../data/simulated_data/generate_inversegauss_distributed_data.R")
#realpars <- generate_gamma_simulated(seed, n, rep)
#filename <- paste0("../data/simulated_data/gamma_n_", n, "_rep_", rep, ".csv")
realpars <- generate_ig_simulated(seed, n, rep)
filename <- paste0("../data/simulated_data/ig_n_", n, "_rep_", rep, ".csv")
source("general.R")
# fit simulations using SSQ
resSSQ <- run_model(datafile = filename, model = "full", goalf = "SSQ", 
                    skipEM = FALSE)
resWLS <- run_model(datafile = filename, model = "full", goalf = "WLS", 
                    skipEM = FALSE)
resGamma <- run_model(datafile = filename, model = "full", goalf = "IG", 
                      skipEM = TRUE)
# store observed and predicted
# compute communities
communities <- community <- apply((resSSQ$observed > 0) * 1, 1, paste, collapse = "")
observed <- resSSQ$observed %>% as_tibble() %>% 
  mutate(experiment = row_number()) %>% add_column(community = community) %>% 
  pivot_longer(names_to = "species", values_to = "abundance", cols = -c(experiment, community))
# extract predictions SSQ
predictedSSQ <- resSSQ$predicted %>% as_tibble() %>% 
  mutate(experiment = row_number()) %>% add_column(community = community) %>% 
  pivot_longer(names_to = "species", values_to = "SSQ", cols = -c(experiment, community))
# extract predictions WLS
predictedWLS <- resWLS$predicted %>% as_tibble() %>% 
  mutate(experiment = row_number()) %>% add_column(community = community) %>% 
  pivot_longer(names_to = "species", values_to = "WLS", cols = -c(experiment, community))
# extract predictions Gamma
predictedGamma <- resGamma$predicted %>% as_tibble() %>% 
  mutate(experiment = row_number()) %>% add_column(community = community) %>% 
  pivot_longer(names_to = "species", values_to = "Gamma", cols = -c(experiment, community))

toplot <- observed %>% inner_join(predictedSSQ) %>% inner_join(predictedWLS)%>% inner_join(predictedGamma)
toplot <- toplot %>% pivot_longer(names_to = "goal function", values_to = "prediction", cols = c(SSQ, WLS, Gamma))
toplotavg <- toplot %>% group_by(community, species, `goal function`) %>% 
  summarise(prediction = mean(prediction), abundance = mean(abundance), .groups = "drop") %>% 
  filter(abundance > 0)

pl <- ggplot(toplotavg) + 
  geom_abline(slope = 1, intercept = 0, linetype = 2) + 
  aes(x = abundance, y= prediction, colour = species) + 
  geom_point() + scale_x_log10() + scale_y_log10() + 
  theme_bw() + theme(legend.position = "bottom")+
  facet_grid(~`goal function`) 
show(pl)
