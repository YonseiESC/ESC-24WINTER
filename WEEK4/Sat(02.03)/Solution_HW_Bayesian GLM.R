### Setup
rm(list=ls())
set.seed(100)
##### Before run this library code, install these libraries. Use "install.packages("name")"
library(coda); library(mvtnorm); library(invgamma)
library(adaptMCMC); library(nimble) ; library(LaplacesDemon)

## Use the canonical link function for Poisson r.v.
inv_link <- function(x){
  ##### TODO
  ##### Complete this inv_link function.
  return( exp(x) )
}

# Data Generating
n=100
X = matrix(NA, nrow=n, ncol=4)
X[,1] = 1
X[,2] = rnorm(n, mean = 0, sd = 1)
X[,3] = rpois(n, lambda = 2)
X[,4] = rbern(n, 0.7)

## Note that value of true_beta is different from Bayesian Logistic Regression example.
true_beta = c(-1, 2, 0.3, -0.4)
true_lambda = exp(X%*%true_beta)
Y = rpois(n, true_lambda)
X
Y


# Declare
n.iter = 20000
beta = matrix(NA, nrow = 4, ncol = n.iter)

# Initialize
beta[,1] = rep(0, 4)

##### TODO (Optional)
##### Research the "Metropolis-Hastings Algorithm" and consider using alternative tuning parameters.
##### Each time you adjust the tuning parameter value, examine the acceptance rate and the trace plot.
##### These are the key points.
##### If you do this, also attempt to set "n.burnin" independently based on your observations.
tuningS=c(0.04) #for-beta



## Kernel
# Prior Layer ; log-prior
# 생략. Use the function logPost directly without specifying it as a log-prior function.

# Likelihood ; Y | lambda(beta, X)
logL <- function(Y, X, beta){
  ##### TODO
  ##### Complete this log-likelihood function.
  lamb = inv_link( X %*% beta )
  return( sum( -lamb + Y*log(lamb) - log(factorial(Y)) ) )
}

# Posterior ; beta | Y, X
##### TODO
##### Complete this log-posterior function.
##### Use prior beta_i ~ N(0, 100) i=1,2,3,4
logPost <- function(Y, X, beta){
  return( logL(Y, X, beta) +
            dnorm(beta[1], mean = 0, sd = 10, log=T) + dnorm(beta[2], mean = 0, sd = 10, log=T) +
            dnorm(beta[3], mean = 0, sd = 10, log=T) + dnorm(beta[4], mean = 0, sd = 10, log=T) )
}


## RUN MCMC
for(i in 2:n.iter){
  if(i%%100==0) print(i)
  
  # beta
  betaPrps = rnorm(4, mean = beta[, i-1], sd = tuningS[1]) #Prps means Proposal.
  ##### TODO
  ##### Fill in the blanks (...), and try to understand it with reference to the code of Bayesian Logistic Regression example.
  ##### I STRONGLY recommend researching the "Metropolis-Hastings Algorithm." It is a fundamentally important algorithm in statistics.
  accept_rate = exp( logPost(Y, X, betaPrps) - logPost(Y, X, beta[,i-1]) )
  
  if( !is.na(accept_rate) & runif(1) < accept_rate ){
    beta[,i] = betaPrps
  } else {
    beta[,i] = beta[,i-1]
  }
  
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
##### TODO (Optional)
##### Modify the tuning parameter value and decide on "n.burnin" yourself, guided by the trace plot. There is no definitive answer.
n.burnin = 5000

beta = beta[,(n.burnin+1):n.iter]

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
