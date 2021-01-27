#test the posterior distribution of two methods (three columns)
rm(list=ls())
library(mice, warn.conflicts = FALSE)
library(miceadds)
library(MASS)
library(magrittr)
library(purrr)
library(MCMCpack)
library(norm)


#the correlation is 1/3
#mu.0 = c(0, 0, 0)
#tau = 1
#m = 3
#lambda = matrix(c(400, 0, 0, 0, 400, 0, 0, 0, 400), 3, 3)


set.seed(123)
rngseed(1234567)
sample.size <- 200
#niter <- 20
nsim <- 1000
priorr <- list(tau = 1, m = 3, mu0 = c(0, 0, 0),lambdainv = matrix(c(1/60, 0, 0, 0, 1/60, 0, 0, 0, 1/60), nrow = 3, ncol = 3))
set.seed(123)
rngseed(1234567)
impute.single <- function(incomplete.data, nit) {
  incomplete.data <- as.matrix(incomplete.data)
  s <- prelim.norm(incomplete.data)
  theta.ini <- em.norm(s, showits = FALSE)
  theta <- da.norm(s, theta.ini, priorr, steps = nit)
  imp.data <- imp.norm(s, theta, incomplete.data)
  return(as.data.frame(imp.data))
  
}



rmvnorm <- function(n, mu, Sigma){
  E <- matrix(rnorm(n * length(mu)), n, length(mu))
  t(t(E %*% chol(Sigma)) + c(mu))
}

bs.draw <- function(y, x, nit, beta.0, Sigma.0, nu.0, sigma2.0, ...){
  invSigma.0 <- solve(Sigma.0)
  #starting value
  sigma2 <- var(residuals(lm(y ~ 0 + x)))
  for (i in 1 : nit) {
    #updat beta
    beta.var <- solve(invSigma.0 + t(x) %*% x / sigma2)
    beta.mean <- beta.var %*% (invSigma.0 %*% beta.0 + t(x) %*% y / sigma2)
    beta <- t(rmvnorm(1, beta.mean, beta.var))
    #update sigma2
    nu.n <- nu.0 + length(y)
    ss.n <- nu.0 * sigma2.0 + sum((y - x %*% beta)^2)
    sigma2 <- 1 / rgamma(1, nu.n / 2, ss.n / 2)
  }
  return(list(beta = beta, sigma = sqrt(sigma2)))
}

mice.impute.norm.bs <- function(y, ry, x, nit, beta.0, Sigma.0, nu.0, sigma2.0, ...){
  wy <- !ry
  x <- cbind(1, as.matrix(x))
  parm <- bs.draw(y[ry], x[ry, ], nit, beta.0, Sigma.0, nu.0, sigma2.0)
  y[wy] <- x[wy, ] %*% parm$beta + rnorm(sum(wy)) * parm$sigma
  return(y[wy]) 
}


OUT.1 <- rep(0, nsim)
OUT.2 <- rep(0, nsim)
complete.data <- data.frame(mvrnorm(n = 200, mu = c(1, 4, 9), Sigma = matrix(c(4, 2, 2, 2, 4, 2, 2, 2, 9), nrow = 3)))
colnames(complete.data) <- c("x", "y", "z")
mypattern <- rbind(c(1, 0, 0), c(0, 1, 0), c(0, 0, 1))
incomplete.data <- ampute(complete.data, prop = 0.5, patterns = mypattern, mech = "MCAR")$amp
wx <- is.na(incomplete.data$x)
wy <- is.na(incomplete.data$y)
wz <- is.na(incomplete.data$z)
imp.data <- incomplete.data
imp.data$x[wx] <- mice.impute.sample(incomplete.data$x, !wx)
imp.data$y[wy] <- mice.impute.sample(incomplete.data$y, !wy)
pb <- txtProgressBar(min = 0, max = nsim, style = 3)
for (i in 1 : nsim) {
  OUT.1[i] <- mean(impute.single(incomplete.data, i)$y)
    imp.data$z[wz] <- mice.impute.norm.bs(incomplete.data$z, !is.na(incomplete.data$z), cbind(imp.data$x, imp.data$y),
                                          50, c(0, 0, 0), matrix(c(60, 0, 0, 0, 3600, 0, 0, 0, 3600), 3, 3),
                                          1.5, 30)
    imp.data$x[wx] <- mice.impute.norm.bs(incomplete.data$x, !is.na(incomplete.data$x), cbind(imp.data$y, imp.data$z),
                                          50, c(0, 0, 0), matrix(c(60, 0, 0, 0, 3600, 0, 0, 0, 3600), 3, 3),
                                          1.5, 30)
    imp.data$y[wy] <- mice.impute.norm.bs(incomplete.data$y, !is.na(incomplete.data$y), cbind(imp.data$x, imp.data$z), 
                                          50, c(0, 0, 0), matrix(c(60, 0, 0, 0, 3600, 0, 0, 0, 3600), 3, 3),
                                          1.5, 30)
  OUT.2[i] <- mean(imp.data$y)
  setTxtProgressBar(pb, i)
  }
close(pb)  

qqplot(OUT.1, OUT.2, xlab = "JM", ylab = "FCS", main = expression("posterior distribution of coefficient B"[1]))
abline(coef = c(0,1), col = "red")
