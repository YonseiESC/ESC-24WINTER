### Setup
rm(list=ls())
set.seed(2021122084)
library(coda); library(mvtnorm); library(invgamma)
library(adaptMCMC); library(nimble) ; library(LaplacesDemon)

## link function
sigmoid <- function(x){
  return(exp(x)/(1+exp(x)))
}

# Data Generating
n=100
X = matrix(NA, nrow=n, ncol=4)
X[,1] = 1
X[,2] = rnorm(n, mean = 0, sd = 1)
X[,3] = rpois(n, lambda = 5)
X[,4] = rbern(n, 0.7)
true_beta = c(-1, 2, 0.3, 4)
Y = rbern(n, sigmoid(X%*%true_beta))
X
Y


# Declare
n.iter = 20000
beta = matrix(NA, nrow = 4, ncol = n.iter)
mu = rep(NA, n.iter)

# Initialize
beta[,1] = rep(0, 4)

tuningS=c(0.2) #beta


## link function
sigmoid <- function(x){
  return(exp(x)/(1+exp(x)))
}

## Kernel
# Prior Layer ; log-prior
# 생략

# Likelihood ; Y | mu(beta, X)
logL <- function(Y, X, beta){
  mu = sigmoid( X %*% beta )
  return( sum( Y*log(mu) + (1-Y)*log(1-mu) ) )
}

# Posterior ; beta | Y, X
## Use prior beta ~ N(0, 100)
logPost <- function(Y, X, beta){
  return( logL(Y, X, beta) +
            dnorm(beta[1], mean = 0, sd = 10, log=T) + dnorm(beta[2], mean = 0, sd = 10, log=T) +
            dnorm(beta[3], mean = 0, sd = 10, log=T) + dnorm(beta[4], mean = 0, sd = 10, log=T) )
}


## RUN MCMC
for(i in 2:n.iter){
  if(i%%100==0) print(i)
  
  # beta
  betaPrps = rnorm(4, mean = beta[, i-1], sd = tuningS[1])
  accept_rate = exp( logPost(Y, X, betaPrps) - logPost(Y, X, beta[,i-1]) )
  
  if( !is.na(accept_rate) & runif(1) < accept_rate ){
    beta[,i] = betaPrps
  } else {
    beta[,i] = beta[,i-1]
  }
  
  mu[i] = sigmoid(X*beta[i])
  
}

#
par(mfrow=c(1,1))
ts.plot(beta[1,])
ts.plot(beta[2,])
ts.plot(beta[3,])
ts.plot(beta[4,])

AcceptanceRate(t(beta))
rowMeans(beta)

# burn-in
nburnin = 5000
beta = beta[, (nburnin+1):n.iter]

ts.plot(beta[1,])
ts.plot(beta[2,])
ts.plot(beta[3,])
ts.plot(beta[4,])
rowMeans(beta)


plot(density(beta[1,]))
abline(v=true_beta[1])

plot(density(beta[2,]))
abline(v=true_beta[2])

plot(density(beta[3,]))
abline(v=true_beta[3])

plot(density(beta[4,]))
abline(v=true_beta[4])
