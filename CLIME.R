# node-wise lasso for high-dimensional sparse corrrelation matrix estimation
library(mvtnorm)
library(glmnet)

p <- 200
n <- 200
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

var.j <- rep(0, p)
for (r in 1:p){
  tmp <- cv.glmnet(xB[, -r], xB[, r], family = 'gaussian',intercept = FALSE)
  gamma.pro[,r] <- as.matrix(tmp$glmnet.fit$beta[,which(tmp$cvm == min(tmp$cvm))])[,1]
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
