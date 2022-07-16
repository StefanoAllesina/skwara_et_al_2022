goal_type <- "IG"

goal_function <- function(X, E, Evar, return_pred = FALSE, use_penalization = TRUE){
  if (return_pred) return(X)
  penalization <- 0
  if (use_penalization){
    penalization <- sum(X < 0) * max(E) * 1000 - sum(X[X < 0])
    X <- abs(X)
  }
  X <- abs(X)
  loglikelihood <- 0
  for (i in 1:nrow(E)){
    pres <- E[i,] > 0
    loglikelihood <- loglikelihood - sum(statmod::dinvgauss(x = E[i, pres], mean = X[i, pres], log = TRUE))
  }
  return(loglikelihood + penalization)
}
