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
n <- 100
d <- 6


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

gamma.pro <- matrix(0, p-1, p)
z.pro <- matrix(0, p-1, p)
nvec.pro <- matrix(0, p-1, p)
XtX <- matrix(0, p, p)

knots.L <- seq(0, 1, length.out=round(n ^ 0.2 + 2))


c1 <- c2 <- c3 <- 0
for (b in 1:10){
  xB <- rmvnorm(n, mean=rep(0,p), sigma=sigma.x)
  tB <- matrix(runif(3*n, 0, 1), n, 3)
  XtX <- XtX + t(xB) %*% xB
  piB <- matrix(0, n, 3*d)
  piB[,1:d] <- bs(tB[,1], df=NULL, knots=knots.L[-c(1, length(knots.L))],
                  degree=3, intercept=F, Boundary.knots=range(tB[,1]))
  c1 <- (c1 * (b-1) + apply(piB[,1:d], 2, mean)) / b
  piB[,1:d] <- scale(piB[,1:d], center=T, scale=F)
  piB[,(d+1):(2*d)] <- bs(tB[,2], df=NULL, knots=knots.L[-c(1, length(knots.L))],
                          degree=3, intercept=F, Boundary.knots=range(tB[,2]))
  c2 <- (c2 * (b-1) + apply(piB[,(d+1):(2*d)], 2, mean)) / b
  piB[,(d+1):(2*d)] <- scale(piB[,(d+1):(2*d)], center=T, scale=F)
  piB[,(2*d+1):(3*d)] <- bs(tB[,3], df=NULL, knots=knots.L[-c(1, length(knots.L))],
                            degree=3, intercept=F, Boundary.knots=range(tB[,3]))
  c3 <- (c3 * (b-1) + apply(piB[,(2*d+1):(3*d)], 2, mean))/ b
  piB[,(2*d+1):(3*d)] <- scale(piB[,(2*d+1):(3*d)], center=T, scale=F)
  
  
  W <- diag(n) - piB %*% solve(t(piB) %*% piB) %*% t(piB)
  
  
  var.j <- rep(0, p)
  for (r in 1:p){
    fit <- ftrl(W %*% xB[, -r], W %*% xB[, r], gamma.pro[,r], lambda1=lambda1* sqrt(n*b), lambda2=lambda2, alpha=alpha, gamma=gamma, z.pro[,r], nvec.pro[,r])
    z.pro[, r] <- fit$z
    nvec.pro[, r] <- fit$nvec
    gamma.pro[,r] <- fit$beta
    var.j[r] <- 1 / (XtX[r, r] - t(gamma.pro[,r]) %*% XtX[-r, r]) 
  }
  tmp <- -1 * matrix(c(0, c(rbind(matrix(c(gamma.pro %*% diag(var.j)) , p, p-1), 0))), p, p)                                    # precision matrix
  tmp <- (abs(tmp) <= abs(t(tmp))) * tmp + (abs(tmp) > abs(t(tmp))) * t(tmp)                        # symmetrilization
  Theta <- tmp + diag(var.j)
}



a <- W %*% xB %*% Theta %*% t(W %*% xB)


library("lattice")
# Dummy data
x.axis <- seq(1,n, length.out=n)
y.axis <- seq(1,n, length.out=n)
data <- expand.grid(X=x.axis, Y=y.axis)
data$Z <- c(a)

## Try it out
levelplot(Z ~ X*Y, data=data  ,xlab="X",
          main="")
