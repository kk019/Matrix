library(mvtnorm)
n <- 5000
x <- rmvnorm(n, rep(0, 10), cov)
x1 <- x[,1:3]
x2 <- x[,-(1:3)]
x3 <- x[,1]
x4 <- x[,-1]
# projection matrix (I-p2)*x1 = x1 - x2 * Gamma
# Gamma is block of inv(cov)
gamma <- c(solve(t(x4) %*% x4) %*% t(x4) %*% x3)  # x1 - x2 * Gamma
I <- diag(n)
p4 <- x4 %*% solve(t(x4) %*% x4) %*% t(x4)
tau <- solve(t(x3) %*% (I - p4) %*% x3 / n)
tau
Theta[1,1]                  # Theta[i,i] = inv(x'(I-P)x)
Theta <- solve(cov)   
Theta[1,]
tau * gamma * (-1)          # Theta[i,-i] = -Theta[i,i] * gamma
