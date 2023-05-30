rm(list = ls())

invlink <- function(X, y, beta){
  X %*% beta
}

ftrl <- function(X, y, beta.ini, lambda1, lambda2, alpha, gamma, z, nvec){
  # storage
  maxiter <- 10
  tol <- 0.01
  for (r in 1:maxiter){
    mu <- invlink(X, y, beta.ini)
    g <- t(X) %*% (mu - y)
    step.sigma <- (sqrt(nvec + g**2) - sqrt(nvec)) / alpha
    z.tmp <- z + g - step.sigma * beta.ini
    nvec.tmp <- nvec + g**2
    
    # update
    beta <- -1 * (abs(z.tmp) > lambda1) * (z.tmp - sign(z.tmp) * lambda1) / ((gamma + sqrt(nvec.tmp)) / alpha + lambda2)
    
    if (sum((beta - beta.ini)**2) < tol){
      break
    }
    beta.ini <- beta
  }
  
  z <- z.tmp
  nvec <- nvec.tmp
  # beta[abs(beta) < 0.2] <- 0
  return(list(z = z, nvec = nvec, beta = beta))
}

ftrl.hard <- function(X, y, beta.ini, lambda1, lambda2, alpha, gamma, z, nvec){
  # storage
  maxiter <- 10
  tol <- 0.01
  for (r in 1:maxiter){
    mu <- invlink(X, y, beta.ini)
    g <- t(X) %*% (mu - y)
    step.sigma <- (sqrt(nvec + g**2) - sqrt(nvec)) / alpha
    z.tmp <- z + g - step.sigma * beta.ini
    nvec.tmp <- nvec + g**2
    
    # update
    beta <- -1 * (abs(z.tmp) > lambda1) * (z.tmp - sign(z.tmp) * lambda1) / ((gamma + sqrt(nvec.tmp)) / alpha + lambda2)
    
    if (sum((beta - beta.ini)**2) < tol){
      break
    }
    beta.ini <- beta
  }
  
  z <- z.tmp
  nvec <- nvec.tmp
  beta[abs(beta) < 0.2] <- 0
  return(list(z = z, nvec = nvec, beta = beta))
}



library(mvtnorm)
library(glmnet)

p <- 200
n <- 200

# for ftrl
lambda1 = 1.2
lambda2 = 1
alpha = 0.1
gamma = 1

sigma.x <- matrix(0, p, p)
for (i in 1:p) {
  for (j in 1:p) {
    sigma.x[i, j] <- 0.5 ^ abs(i - j)
  }
}
xB <- rmvnorm(n, mean=rep(0,p), sigma=sigma.x)
gamma.pro <- matrix(0, p-1, p)
z.pro <- matrix(0, p-1, p)
nvec.pro <- matrix(0, p-1, p)

XtX <- t(xB) %*% xB

b <- 1
var.j <- rep(0, p)
for (r in 1:p){
  fit <- ftrl(xB[, -r], xB[, r], gamma.pro[,r], lambda1=lambda1* sqrt(n*b), lambda2=lambda2, alpha=alpha, gamma=gamma, z.pro[,r], nvec.pro[,r])
  z.pro[, r] <- fit$z
  nvec.pro[, r] <- fit$nvec
  gamma.pro[,r] <- fit$beta
  var.j[r] <- 1 / (XtX[r, r] - t(gamma.pro[,r]) %*% XtX[-r, r]) 
}
tmp <- -1 * matrix(c(0, c(rbind(matrix(c(gamma.pro %*% diag(var.j)) , p, p-1), 0))), p, p)                                    # precision matrix
tmp <- (abs(tmp) <= abs(t(tmp))) * tmp + (abs(tmp) > abs(t(tmp))) * t(tmp)                        # symmetrilization
Theta <- tmp + diag(var.j)


a <- xB %*% Theta %*% t(xB)


library("lattice")
# Dummy data
x.axis <- seq(1,n, length.out=n)
y.axis <- seq(1,n, length.out=n)
data <- expand.grid(X=x.axis, Y=y.axis)
data$Z <- c(xB %*% Theta %*% t(xB))

## Try it out
levelplot(Z ~ X*Y, data=data  ,xlab="X",
          main="")
