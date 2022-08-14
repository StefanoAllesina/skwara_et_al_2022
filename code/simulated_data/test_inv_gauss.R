# Goal:
# simulate data in which the "observed" abundances are noisy estimates of a true mean
# x_i^(k) following InvGauss(x_i^(k), 1)
# - recover the approximate parameters by maximizing the likelihood using the Gamma density function
# - contrast with results obtained using SSQ (Normal) or WLS (also Normal)

library(statmod) # for inverse Gaussian
set.seed(7)
# get prediction given matrix B, number of species, number of replicates for community
getpred <- function(B, n, nrep){
  results <- matrix(0, 0, n)
  for (i in 1:(2^n - 1)){
    presence <- as.numeric(intToBits(i)[1:n])
    Bk <- B[presence > 0, presence > 0, drop = FALSE]
    x <- rep(0, n)
    x[presence > 0] <- rowSums(solve(Bk))
    results <- rbind(results, x)
  }
  results <- results[sort(rep(1:nrow(results), nrep)),]
  colnames(results) <- letters[1:n]
  return(results)
}

# create labels for community
getcomm <- function(X){
  comm <- apply((X > 0) * 1, 1, paste, collapse = "")
  return(comm)
}

# 1) Choose matrix B yielding coexistence for all combinations
#    and with a large spread between abundances (to test SSQ)      
n <- 3 # number of species
nrep <- 25 # number of replicates for each community (for better results choose many reps)
B <- matrix(c(10,0.5,-1,
              -4,1,-0.7,
              0.5,0.5,0.01), 3, 3, byrow = TRUE) / 10
truemeans <- getpred(B, n, nrep)

# save the "observed" data
samples <- apply(truemeans, c(1,2), function(x) rinvgauss(1, x, 1))
samples[is.na(samples)] <- 0 
write.csv(samples, file = "invgauss_samples.csv", row.names = FALSE)

# 2) Now fit using SSQ and WLS
setwd("../")
source("general.R")
tmp_ssq <- run_model("simulated_data/invgauss_samples.csv", model = "full", skipEM = FALSE)
tmp_wls <- run_model("simulated_data/invgauss_samples.csv", model = "full", goalf = "WLS", skipEM = FALSE)
setwd("simulated_data/")

# 3) Recover parameters using Gamma distribution
# compute likelihood given parameters
getlik <- function(pars, n, nrep){
  B1 <- matrix(pars, n, n)
  newmeans <- getpred(B1, n, nrep)
  penalization <- -sum(newmeans < 0) * 1000
  newmeans <- abs(newmeans)
  ll <- dinvgauss(samples, newmeans, 1, log = TRUE)
  ll[is.infinite(ll)] <- 0
  return(-sum(ll) - penalization)
}

# likelihood with true parameters
print(getlik(as.vector(B), n, nrep))
tmp <- optim(par = as.vector(B), n = n, nrep = nrep, fn = getlik, method = "Nelder-Mead", control = list(maxit = 5000, trace = TRUE))
tmp <- optim(par = tmp$par, n = n, nrep = nrep, fn = getlik, method = "BFGS", control = list(maxit = 5000, trace = TRUE))
tmp <- optim(par = tmp$par, n = n, nrep = nrep, fn = getlik, method = "Nelder-Mead", control = list(maxit = 5000, trace = TRUE))
tmp <- optim(par = tmp$par, n = n, nrep = nrep, fn = getlik, method = "BFGS", control = list(maxit = 5000, trace = TRUE))
tmp <- optim(par = tmp$par, n = n, nrep = nrep, fn = getlik, method = "Nelder-Mead", control = list(maxit = 5000, trace = TRUE))

Binvgauss <- matrix(tmp$par, n, n)
pred_invgauss <- getpred(Binvgauss, n, nrep)
plot(as.vector(B), as.vector(Binvgauss)); abline(c(0,1))

# now compute means for observed and predictions
library(tidyverse)
getmeans <- function(X){
  comms <- getcomm(X)
  X <- as_tibble(X) %>% mutate(experiment = row_number()) %>% add_column(comm = comms)
  means <- X %>% pivot_longer(names_to = "species", values_to = "abundance", 
                              cols = -c(experiment, comm))
  means <- means %>% group_by(species, comm) %>% 
    summarize(`mean abundance` = mean(abundance), .groups = "drop") %>% filter(`mean abundance` > 0)
  return(means)
}

# observed means
toplot <- getmeans(samples)
toplot <- toplot %>% rename(`observed mean abundance` = `mean abundance`)
# means according to SSQ
ssqm <- getmeans(tmp_ssq$predicted)
ssqm <- ssqm %>% rename(`SSQ` = `mean abundance`)
# means according to WLS
wlsm <- getmeans(tmp_wls$predicted)
wlsm <- wlsm %>% rename(`WLS` = `mean abundance`)
# means according to Inverse Gaussian
invgaussm <- getmeans(pred_invgauss)
invgaussm <- invgaussm %>% rename(`Inverse Gaussian` = `mean abundance`)
# combine data
toplot <- toplot %>% inner_join(ssqm, by = c("species", "comm")) %>% 
  inner_join(wlsm,  by = c("species", "comm")) %>% 
  inner_join(invgaussm, by = c("species", "comm"))
# make tidy
toplot <- toplot %>% pivot_longer(names_to = "method", 
                                  values_to = " predicted mean abundance", 
                                  cols = c(SSQ, WLS, `Inverse Gaussian`))

pl <- ggplot(toplot, aes(x = `observed mean abundance`, y = ` predicted mean abundance`, colour = method)) + 
  geom_abline(slope = 1, intercept = 0) + 
  geom_point() +  facet_wrap(~method, scales = "free") + 
  theme_bw() + scale_x_log10() + scale_y_log10() + 
  scale_color_brewer(palette = "Set1") + 
  theme(legend.position = "none")
show(pl)

ggsave(pl, filename = "InvGauss_predabundance_v_observed.pdf", width = 7, height = 3)

# FIGURE RECOVER PARAMETERS
Bssq <- solve(matrix(tmp_ssq$pars, n, n))
Bwls <- solve(matrix(tmp_wls$pars, n, n))

toplot2 <- tibble(`true parameters` = as.vector(B),
                  `inferred parameters` = as.vector(Binvgauss), 
                  `goalf` = "InvGauss")
toplot2 <- toplot2 %>% bind_rows(
  tibble(`true parameters` = as.vector(B),
         `inferred parameters` = as.vector(Bssq), 
         `goalf` = "SSQ"))
toplot2 <- toplot2 %>% bind_rows(
  tibble(`true parameters` = as.vector(B),
         `inferred parameters` = as.vector(Bwls), 
         `goalf` = "WLS"))


pl2 <- ggplot(toplot2) + 
  aes(x = `true parameters`, y = `inferred parameters`, colour = goalf) + 
  geom_abline(slope = 1, intercept = 0) + 
  geom_point() +  facet_wrap(~goalf, scales = "free") + 
  theme_bw() +
  scale_color_brewer(palette = "Set1") + 
  theme(legend.position = "none")
show(pl2)

ggsave(pl2, filename = "InvGauss_pars.pdf", width = 7, height = 3)
