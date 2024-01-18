# Data Sampling
n_sample <- 100
X_par <- c(0, 2) # mean, sd

set.seed(2024)
X <- rnorm(n_sample, X_par[1], X_par[2])
MLE <- c(mean(X), sqrt(sum((X-mean(X))^2)/length(X)))
plot(MLE[1], 0.3, xlim = c(X_par[1]-X_par[2]-0.1, X_par[1]+X_par[2]+0.1), ylim = c(-.1, .5), col='red', pch=19, cex=1,
    main = 'Robbins-Monro Algorithm', xlab = 'MLE for Mean', ylab = 'Estimates')

# Robbins-Monro
mean_update <- function(x, mean_old, sd_old, coef){
    var_old <- sd_old^2
    result <- mean_old + coef * (1/var_old) * (x - mean_old)
    return(result)
}

sd_update <- function(x, mean, sd_old, coef){
    z_fct <- function(x, mean, sd_old){
        temp <- ((x-mean)^2-sd_old^2)/(sd_old^3)
        return(temp)
    }
    result <- sd_old + coef * z_fct(x, mean, sd_old)
    return(result)
}

mean_history <- c()
sd_history <- c()

for(i in 1:n_sample){
    mean_old <- ifelse(length(mean_history) == 0, 0, mean_history[length(mean_history)])
    sd_old <- ifelse(length(sd_history) == 0, 1, sd_history[length(sd_history)])

    mean_new <- mean_update(X[i], mean_old = mean_old, sd_old = sd_old, coef = 1/i)
    sd_new <- sd_update(X[i], mean = mean_new, sd_old = sd_old, coef = 1/i)
    
    mean_history <- c(mean_history, mean_new)
    update_sd <- c(sd_history, sd_new)

    col_mean <- hsv((i/n_sample)^0.7, 0.8, 0.8)
    col_sd <- hsv(0.8, (i/n_sample)^0.7, 0.8)
    if(i != n_sample){
        points(mean_new, 0, col = col_mean, pch = 19, cex = .5)
    } else{
        abline(v = mean_new, col = 'red')
    }
}

# TS plot
title <- expression(paste(mu[ML], ' Update History'))
plot(1:length(mean_history), mean_history, type='n', xlim=c(1, length(mean_history)), ylim=c(min(mean_history), max(mean_history)), xlab='Iteration', ylab='MLE for Mean', main=title)

lines(1:length(mean_history), mean_history, col='black')
abline(h = MLE[1], col = 'red', lwd=1.5)
