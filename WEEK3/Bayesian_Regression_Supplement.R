# setup
library(MASS)
library(tidyverse)
library(gridExtra)
library(emdbook)
set.seed(2024)

a0 <- -.3; a1 <- .5
target <- function(x, a0, a1){
    noise <- rnorm(1, 0, 0.2)
    target <- a0+a1*x+noise
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

lay <- matrix(c(NA, 1, 2), nrow = 1, byrow = T)
grid.arrange(p.1_2, p.1_3, layout_matrix = lay,
             widths = c(1,1,1), heights = c(1,1,1))

# 2nd row
x.2 <- runif(1, -1, 1)
t.2 <- target(x.2, a0, a1)

# (2,1)
likelihood <- function(t, x, w0, w1, beta){
  mu <- w0 + x*w1
  likelihood <- dnorm(t, mean = mu, sd = 1/sqrt(beta))
  return(likelihood)
}
grid.2_1 <- expand.grid(w0 = seq(-1, 1, length.out = 200), w1 = seq(-1, 1, length.out = 200))
grid.2_1 <- cbind(t = t.2, x = x.2, grid.2_1, beta = beta)
likelihood <- likelihood(grid.2_1$t, grid.2_1$x, grid.2_1$w0, grid.2_1$w1, grid.2_1$beta)
df.2_1  <- cbind(grid.2_1, likelihood)

breaks.2_1 <- seq(min(df.2_1$likelihood), max(df.2_1$likelihood), by = .01)
p.2_1 <- ggplot(df.2_1, aes(x = w0, y = w1, z = likelihood)) +
  geom_contour_filled(breaks = breaks.2_1, show.legend = F) +
  geom_point(x = x.2, y = t.2, shape = 18, size = 3) +
  coord_fixed(xlim = c(-1, 1), ylim = c(-1, 1), ratio = 1) +
  theme(legend.position = "none") +
  labs(x = "w0", y = "w1")

# (2,2)
phi <- matrix(c(1, x.2), nrow = 1)
posterior_precision <- (alpha * diag(2)) + (beta * t(phi) %*% phi)
posterior_cov <- solve(posterior_precision)
posterior_mean <- as.vector(beta * (posterior_cov %*% t(phi)) * t.2)

grid.2_2 <- expand.grid(w0 = seq(-1, 1, length.out = 200), w1 = seq(-1, 1, length.out = 200))
df.2_2 <- cbind(grid.2_2, prob = dmvnorm(as.matrix(grid.2_2), posterior_mean, posterior_cov))
breaks.2_2 <- seq(min(df.2_2$prob), max(df.2_2$prob), by = .01)
p.2_2 <- ggplot(df.2_2, aes(x = w0, y = w1, z = prob)) +
  geom_contour_filled(breaks = breaks.2_2, show.legend = F) +
  geom_point(x = x.2, y = t.2, shape = 18, size = 3) +
  coord_fixed(xlim = c(-1, 1), ylim = c(-1, 1), ratio = 1) +
  theme(legend.position = "none") +
  labs(x = "w0", y = "w1")

# (2,3)
df.2_3 <- data.frame(idx = character(), x = numeric(), y = numeric())

for(i in 1:6){
  w <- mvrnorm(n = 1, posterior_mean, posterior_cov)
  df.2_3 <- rbind(df.2_3, data.frame(idx = rep(paste("y", i, sep = "_"), 100), x = x, y = w[1] + x * w[2]))
}
df.2_3$idx <- as.factor(df.2_3$idx)

p.2_3 <- ggplot(df.2_3, aes(x = x, y = y, color = idx)) +
  geom_line() +
  geom_point(x = x.2, y = t.2, color = "blue", shape = 1, size = 3, ) +
  labs(x = "x", y = "y", title = "data space") +
  scale_color_manual(values = rep("red", 6)) +
  theme(legend.position = "none") +
  xlim(-1, 1) + ylim(-1, 1) +
  geom_line(data = df, aes(x = x, y = y), color = "black", size = 1.5)

lay <- matrix(c(NA, 1, 2, 3, 4, 5), nrow = 2, byrow = T)
grid.arrange(p.1_2, p.1_3, p.2_1, p.2_2, p.2_3, layout_matrix = lay,
             widths = c(1,1,1), heights = c(1,1,1))

# 3rd row
x.3 <- runif(1, -1, 1)
t.3 <- target(x.3, a0, a1)

# (3,1)
grid.3_1 <- expand.grid(w0 = seq(-1, 1, length.out = 200), w1 = seq(-1, 1, length.out = 200))
grid.3_1 <- cbind(t = t.3, x = x.3, grid.3_1, beta = beta)
likelihood <- likelihood(grid.3_1$t, grid.3_1$x, grid.3_1$w0, grid.3_1$w1, grid.3_1$beta)
df.3_1  <- cbind(grid.3_1, likelihood)

breaks.3_1 <- seq(min(df.3_1$likelihood), max(df.3_1$likelihood), by = .01)
p.3_1 <- ggplot(df.3_1, aes(x = w0, y = w1, z = likelihood)) +
  geom_contour_filled(breaks = breaks.3_1, show.legend = F) +
  geom_point(x = x.3, y = t.3, shape = 18, size = 3) +
  coord_fixed(xlim = c(-1, 1), ylim = c(-1, 1), ratio = 1) +
  theme(legend.position = "none") +
  labs(x = "w0", y = "w1")
  
# (3,2)
phi.3 <- matrix(c(1, x.3), nrow = 1)
posterior_precision.3 <- posterior_precision + (beta * t(phi.3) %*% phi.3)
posterior_cov.3 <- solve(posterior_precision.3)
posterior_mean.3 <- as.vector(posterior_cov.3 %*% (beta * t(phi.3) %*% t.3))

grid.3_2 <- expand.grid(w0 = seq(-1, 1, length.out = 200), w1 = seq(-1, 1, length.out = 200))
df.3_2 <- cbind(grid.3_2, prob = dmvnorm(as.matrix(grid.3_2), posterior_mean.3, posterior_cov.3))
breaks.3_2 <- seq(min(df.3_2$prob), max(df.3_2$prob), by = .01)
p.3_2 <- ggplot(df.3_2, aes(x = w0, y = w1, z = prob)) +
  geom_contour_filled(breaks = breaks.3_2, show.legend = F) +
  geom_point(x = x.3, y = t.3, shape = 18, size = 3) +
  coord_fixed(xlim = c(-1, 1), ylim = c(-1, 1), ratio = 1) +
  theme(legend.position = "none") +
  labs(x = "w0", y = "w1")

# (3,3)
df.3_3 <- data.frame(idx = character(), x = numeric(), y = numeric())

for(i in 1:6){
  w <- mvrnorm(n = 1, posterior_mean.3, posterior_cov.3)
  df.3_3 <- rbind(df.3_3, data.frame(idx = rep(paste("y", i, sep = "_"), 100), x = x, y = w[1] + x * w[2]))
}
df.3_3$idx <- as.factor(df.3_3$idx)

p.3_3 <- ggplot(df.3_3, aes(x = x, y = y, color = idx)) +
  geom_line() +
  geom_point(x = x.2, y = t.2, color = "blue", shape = 1, size = 3, ) +
  geom_point(x = x.3, y = t.3, color = "blue", shape = 1, size = 3, ) +
  labs(x = "x", y = "y", title = "data space") +
  scale_color_manual(values = rep("red", 6)) +
  theme(legend.position = "none") +
  xlim(-1, 1) + ylim(-1, 1) +
  geom_line(data = df, aes(x = x, y = y), color = "black", size = 1.5)

lay <- matrix(c(NA, 1, 2, 3, 4, 5, 6, 7, 8), nrow = 3, byrow = T)
grid.arrange(p.1_2, p.1_3,
             p.2_1, p.2_2, p.2_3,
             p.3_1, p.3_2, p.3_3,
             layout_matrix = lay, widths = c(1,1,1), heights = c(1,1,1))
