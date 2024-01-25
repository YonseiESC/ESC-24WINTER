# setup
library(MASS)
library(tidyverse)
library(gridExtra)
library(emdbook)
set.seed(2024)

a0 <- -.3; a1 <- .5
target <- function(x, coef=a0, slope=a1){
  noise <- rnorm(1, 0, 0.2)
  target <- a0 + a1 * x + noise
  return(target)
}
x <- seq(-1, 1, length.out = 1e3)
df <- data.frame(x = x, y = (a0 + x * a1))

alpha <- 2; beta <- 25

# 1st row
theme_update(plot.title = element_text(hjust = 0.5))
mean_vector <- c(0, 0)
covariance_matrix <- diag(2) / alpha

df.1_3 <- data.frame(idx = character(), x = numeric(), y = numeric())

for(i in 1:6){
  w <- mvrnorm(n = 1, mu = mean_vector, Sigma = covariance_matrix)
  df.1_3 <- rbind(df.1_3, data.frame(idx = rep(paste("y", i, sep = "_"), 100), x = x, y = w[1] + x * w[2]))
}
df.1_3$idx <- as.factor(df.1_3$idx)

# (1,2)
grid.1_2 <- expand.grid(w0 = seq(-1, 1, length.out = 200), w1 = seq(-1, 1, length.out = 200))
df.1_2  <- cbind(grid.1_2, prob = dmvnorm(as.matrix(grid.1_2), mean_vector, covariance_matrix))
breaks.1_2 <- seq(min(df.1_2$prob), max(df.1_2$prob), by = .01)
p.1_2 <- ggplot(df.1_2, aes(x = w0, y = w1, z = prob)) +
  geom_contour_filled(breaks = breaks.1_2, show.legend = F) +
  coord_fixed(xlim = c(-1, 1), ylim = c(-1, 1), ratio = 1) +
  theme(legend.position = "none") +
  labs(x = "w0", y = "w1", title = "prior/posterior")

# (1,3)
p.1_3 <- ggplot(df.1_3, aes(x = x, y = y, color = idx))+
  geom_line() +
  labs(x = "x", y = "y", title="data space") +
  scale_color_manual(values = rep("red", 6)) +
  theme(legend.position = "none") +
  xlim(-1,1) + ylim(-1, 1) +
  geom_line(data = df, aes(x = x, y = y), color = "black", size = 1.5)

# lay <- matrix(c(NA, 1, 2), nrow = 1, byrow = T)
# grid.arrange(p.1_2, p.1_3, layout_matrix = lay,
             # widths = c(1,1,1), heights = c(1,1,1))

x.hist <- c(0)
t.hist <- c(0)

posterior.precision <- list(alpha*diag(2))
posterior.cov <- list(solve(alpha*diag(2)))
posterior.mean <- list(c(0,0))

likelihood <- list(0)
posterior <- list()
data.space <- list(df.1_3)


update <- function(posterior.precision, posterior.cov, posterior.mean, x){
  # sampling
  x.new <- runif(1, -1, 1)
  t.new <- target(x.new)
  x.hist <<- c(x.hist, x.new)
  t.hist <<- c(t.hist, t.new)
  
  # likelihood
  update_likelihood <- function(t, x, w0, w1, beta){
    mu <- w0 + x*w1
    likelihood <- dnorm(t, mean = mu, sd = 1/sqrt(beta))
    return(likelihood)
  }
  
  grid <- expand.grid(w0 = seq(-1, 1, length.out = 200), w1 = seq(-1, 1, length.out = 200))
  grid <- cbind(t = t.new, x = x.new, grid, beta = beta)
  likelihood <- update_likelihood(grid$t, grid$x, grid$w0, grid$w1, grid$beta)
  likelihood <- cbind(grid, likelihood)
  breaks.likelihood <- seq(min(likelihood$likelihood), max(likelihood$likelihood), by = .01)
  
  # posterior
  phi <- matrix(c(1, x.new), nrow = 1)
  precision <- posterior.precision[[length(posterior.precision)]] + (beta * t(phi) %*% phi)
  posterior.precision[[length(posterior.precision)+1]] <<- precision
  posterior.cov[[length(posterior.cov)+1]] <<- solve(precision)
  mean <- as.vector(solve(precision) %*% (beta * t(phi) %*% target(x.new)))
  posterior.mean[[length(posterior.mean)+1]] <<- mean
  
  grid <- expand.grid(w0 = seq(-1, 1, length.out = 200), w1 = seq(-1, 1, length.out = 200))
  posterior <- cbind(grid, prob = dmvnorm(as.matrix(grid), mean, solve(precision)))
  breaks.posterior <- seq(min(posterior$prob), max(posterior$prob), by = .01)
  
  # data space
  data.space <- data.frame(idx = character(), x = numeric(), y = numeric())
  for(i in 1:6){
    w <- mvrnorm(n = 1, mean, solve(precision))
    data.space <- rbind(data.space, data.frame(idx = rep(paste("y", i, sep = "_"), 100), x = x, y = w[1] + x * w[2]))
  }
  data.space$idx <- as.factor(data.space$idx)
  
  result <- list("likelihood" = likelihood,
                 "breaks.likelihood" = breaks.likelihood,
                 "posterior" = posterior,
                 "breaks.posterior" = breaks.posterior,
                 "data.space" = data.space)
  return(result)
}

result.2 <- update(posterior.precision, posterior.cov, posterior.mean, x)
points <- data.frame(x = x.hist, t = t.hist)

p.2_1 <- ggplot(result.2$likelihood, aes(x = w0, y = w1, z = likelihood)) +
  geom_contour_filled(breaks = result.2$breaks.likelihood, show.legend = F) +
  geom_point(x = x.hist[2], y = t.hist[2], shape = 18, size = 3, col = "white") +
  coord_fixed(xlim = c(-1, 1), ylim = c(-1, 1), ratio = 1) +
  theme(legend.position = "none") +
  labs(x = "w0", y = "w1")
p.2_2 <- ggplot(result.2$posterior, aes(x = w0, y = w1, z = prob)) +
  geom_contour_filled(breaks = result.2$breaks.posterior, show.legend = F) +
  geom_point(x = x.hist[2], y = t.hist[2], shape = 18, size = 3, col = "white") +
  coord_fixed(xlim = c(-1, 1), ylim = c(-1, 1), ratio = 1) +
  theme(legend.position = "none") +
  labs(x = "w0", y = "w1")
p.2_3 <- ggplot() +
  geom_line(data = result.2$data.space, aes(x = x, y = y, color = idx)) +
  geom_point(data = points, aes(x = x, y = t), color = "blue", shape = 1, size = 3, ) +
  labs(x = "x", y = "y") +
  scale_color_manual(values = rep("red", 6)) +
  theme(legend.position = "none") +
  xlim(-1, 1) + ylim(-1, 1) +
  geom_line(data = df, aes(x = x, y = y), color = "black", size = 1.5)

# result.6
result.3 <- update(posterior.precision, posterior.cov, posterior.mean, x)
result.4 <- update(posterior.precision, posterior.cov, posterior.mean, x)
result.5 <- update(posterior.precision, posterior.cov, posterior.mean, x)
result.6 <- update(posterior.precision, posterior.cov, posterior.mean, x)
points <- data.frame(x = x.hist, t = t.hist)

p.3_1 <- ggplot(result.6$likelihood, aes(x = w0, y = w1, z = likelihood)) +
  geom_contour_filled(breaks = result.6$breaks.likelihood, show.legend = F) +
  geom_point(x = x.hist[length(x.hist)], y = t.hist[length(t.hist)], shape = 18, size = 3, col = "white") +
  coord_fixed(xlim = c(-1, 1), ylim = c(-1, 1), ratio = 1) +
  theme(legend.position = "none") +
  labs(x = "w0", y = "w1")
p.3_2 <- ggplot(result.6$posterior, aes(x = w0, y = w1, z = prob)) +
  geom_contour_filled(breaks = result.6$breaks.posterior, show.legend = F) +
  geom_point(x = x.hist[length(x.hist)], y = t.hist[length(t.hist)], shape = 18, size = 3, col = "white") +
  coord_fixed(xlim = c(-1, 1), ylim = c(-1, 1), ratio = 1) +
  theme(legend.position = "none") +
  labs(x = "w0", y = "w1")
p.3_3 <- ggplot() +
  geom_line(data = result.6$data.space, aes(x = x, y = y, color = idx)) +
  geom_point(data = points, aes(x = x, y = t), color = "blue", shape = 1, size = 3, ) +
  labs(x = "x", y = "y") +
  scale_color_manual(values = rep("red", 6)) +
  theme(legend.position = "none") +
  xlim(-1, 1) + ylim(-1, 1) +
  geom_line(data = df, aes(x = x, y = y), color = "black", size = 1.5)

# result.12
result.7 <- update(posterior.precision, posterior.cov, posterior.mean, x)
result.8 <- update(posterior.precision, posterior.cov, posterior.mean, x)
result.9 <- update(posterior.precision, posterior.cov, posterior.mean, x)
result.10 <- update(posterior.precision, posterior.cov, posterior.mean, x)
result.11 <- update(posterior.precision, posterior.cov, posterior.mean, x)
result.12 <- update(posterior.precision, posterior.cov, posterior.mean, x)
points <- data.frame(x = x.hist, t = t.hist)


p.4_1 <- ggplot(result.12$likelihood, aes(x = w0, y = w1, z = likelihood)) +
  geom_contour_filled(breaks = result.12$breaks.likelihood, show.legend = F) +
  geom_point(x = x.hist[length(x.hist)], y = t.hist[length(t.hist)], shape = 18, size = 3, col = "white") +
  coord_fixed(xlim = c(-1, 1), ylim = c(-1, 1), ratio = 1) +
  theme(legend.position = "none") +
  labs(x = "w0", y = "w1")
p.4_2 <- ggplot(result.12$posterior, aes(x = w0, y = w1, z = prob)) +
  geom_contour_filled(breaks = result.12$breaks.posterior, show.legend = F) +
  geom_point(x = x.hist[length(x.hist)], y = t.hist[length(t.hist)], shape = 18, size = 3, col = "white") +
  coord_fixed(xlim = c(-1, 1), ylim = c(-1, 1), ratio = 1) +
  theme(legend.position = "none") +
  labs(x = "w0", y = "w1")
p.4_3 <- ggplot() +
  geom_line(data = result.12$data.space, aes(x = x, y = y, color = idx)) +
  geom_point(data = points, aes(x = x, y = t), color = "blue", shape = 1, size = 3, ) +
  labs(x = "x", y = "y") +
  scale_color_manual(values = rep("red", 6)) +
  theme(legend.position = "none") +
  xlim(-1, 1) + ylim(-1, 1) +
  geom_line(data = df, aes(x = x, y = y), color = "black", size = 1.5)

lay <- matrix(c(NA, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11), nrow = 4, byrow = T)
grid.arrange(p.1_2, p.1_3,
             p.2_1, p.2_2, p.2_3,
             p.3_1, p.3_2, p.3_3,
             p.4_1, p.4_2, p.4_3,
             layout_matrix = lay, widths = c(1,1,1), heights = c(1,1,1,1))
