
# negative log-liklihood of Ricker model with Poisson error
ricker_Pois_neg_ll <- function(theta, y, X){
  
  n <- length(y)
  # compute means
  eta <- vector(mode = "double", length = n)
  eta[1] <- log(y[1])
  for(t in 2:n){
    eta[t] <- log(y[t - 1]) + X[t - 1, ] %*% theta
  }
  # return the negative log-likelihood
  return(-sum(
    dpois(x = y[2:n], lambda = exp(eta[2:n]), log = T)
  ))
}

# check the function by simulating a Ricker model
r <- log(1.5)
alpha <- 0.01

steps <- 100

y <- vector(mode = "double", length = steps)

y[1] <- rpois(1, 50)

for(t in 2:steps){
  
  mu_t <- y[t - 1] * exp(r - alpha * y[t - 1])
  
  y[t] <- rpois(1, lambda = mu_t)
  
}

# provide X matrix
X <- cbind(
  rep(1, steps),
  y
)


# Find MLEs
fit <- optim(par = c(1, 0), fn = neg_ll, y = y, X = X, method = "BFGS", hessian = T)

# standard errors of the two parameters can be found using the square root of the diagonal
# of the inverse of the Hessian matrix
ses <- sqrt(diag(solve(fit$hessian)))

# This is the same as fitting a glm with an offset, but more flexible for when we start deleting
# observations in the time series (Though the function still needs to be adapted to handle that)
summary(
  glm(
    y[2:steps] ~ y[1:(steps - 1)], 
    family = "poisson", 
    offset = log(y[1:(steps-1)])
  )
)








