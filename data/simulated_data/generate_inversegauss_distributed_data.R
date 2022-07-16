# Goal:
# generate synthetic data in which 
# \tilde{x}_i^(k, r) ~ inversegauss(x_i^(k), 1)

# n number of species in the pool
# rep number of replicates per community
generate_ig_simulated <- function(seed, n, rep){
# choose a random seed for reproducibility
set.seed(seed)
# now find random matrix giving rise to all positive predictions
# matrix of means
success <- FALSE
while(!success){
  B <- matrix(rnorm(n * n), n, n)
  # make the monoculture differ of orders of magnitude
  diag(B) <-1/(2^((-1):(n-2)))
  X <- matrix(0,0,n)
  for (i in 1:(2^n - 1)){
    presence <- as.integer(intToBits(i)[1:n]) > 0
    Bk <- B[presence, presence, drop = FALSE]
    xk <- rep(-1, nrow(Bk))
    try({xk <- rowSums(solve(Bk))}, silent = TRUE)
    x <- rep(0, n)
    x[presence] <- xk
    X <- rbind(X, x)
  }
  if (all(X[X!=0] >= 0.1 )) success <- TRUE
}
community <- apply((X > 0) * 1, 1, paste, collapse = "")
X <- X[order(community),]


# Now sample from IG distributions with mean x^(k)_i and dispersion 1
Xs <- matrix(0, 0, n)
Xg <- matrix(0, 0, n)
for (i in 1:nrow(X)){
  x <- X[i, ]
  xs <- x[x > 0]
  for (j in 1:rep){
    xsr <- statmod::rinvgauss(n = length(xs), mean = xs, dispersion = 1)
    tmp <- x
    tmp[x > 0] <- xsr
    Xs <- rbind(Xs, tmp)
    Xg <- rbind(Xg, x)
  }
}

write.csv(Xs, file = paste0("../data/simulated_data/ig_n_", n, "_rep_", rep, ".csv"), row.names = FALSE)
return(list(Xtrue = X, Xg = Xg, Xs = Xs, Btrue = B))
}