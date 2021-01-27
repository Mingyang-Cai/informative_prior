#test the order effect (three columns)
rm(list=ls())
library(mice, warn.conflicts = FALSE)
library(miceadds)
library(MASS)
library(magrittr)
library(purrr)
library(rstanarm)



#the correlation is 1/3
#mu.0 = c(0, 0, 0)
#tau = 1
#m = 3
#lambda = matrix(c(60, 0, 0, 0, 60, 0, 0, 0, 60), 3, 3)


#prior:
#nu.0 <- 1
#sigma2.0 <- 60 / 2
#beta.0 <- c(0, 0, 0)
#Sigma.0 <- matrix(c(60, 0, 0, 0, 3600, 0, 0, 0, 3600), 3, 3)



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
set.seed(123)
count <- rep(0, 500)
beta.tar.mmean <- rep(0, 500)
pb <- txtProgressBar(min = 0, max = 500, style = 3)
for (j in 1 : 500) {
  complete.data <- data.frame(mvrnorm(n = 200, mu = c(1, 4, 9), Sigma = matrix(c(4, 2, 2, 2, 4, 2, 2, 2, 9), nrow = 3)))
  colnames(complete.data) <- c("x", "y", "z")
  mypattern <- rbind(c(1, 0, 0), c(0, 1, 0), c(0, 0, 1))
  incomplete.data <- ampute(complete.data, prop = 0.5, patterns = mypattern, mech = "MCAR")$amp
  beta.tar <- rep(0, 1010)
  beta.1 <- rep(0, 1010)
  beta.2 <- rep(0, 1010)
  batch.mean <- rep(0, 50)
  wx <- is.na(incomplete.data$x)
  wy <- is.na(incomplete.data$y)
  wz <- is.na(incomplete.data$z)
  imp.data <- incomplete.data
  imp.data$x[wx] <- mice.impute.sample(incomplete.data$x, !wx)
  imp.data$y[wy] <- mice.impute.sample(incomplete.data$y, !wy)
  for (i in 1 : 1010) {
    imp.data$z[wz] <- mice.impute.norm.bs(incomplete.data$z, !is.na(incomplete.data$z), cbind(imp.data$x, imp.data$y),
                                          50, c(0, 0, 0), matrix(c(60, 0, 0, 0, 3600, 0, 0, 0, 3600), 3, 3),
                                          1, 30)
    beta.1[i] <- lm(imp.data$y~imp.data$x + imp.data$z)$coef[2]
    imp.data$x[wx] <- mice.impute.norm.bs(incomplete.data$x, !is.na(incomplete.data$x), cbind(imp.data$y, imp.data$z),
                                          50, c(0, 0, 0), matrix(c(60, 0, 0, 0, 3600, 0, 0, 0, 3600), 3, 3),
                                          1, 30)
    beta.2[i] <- lm(imp.data$y~imp.data$x + imp.data$z)$coef[2]
    imp.data$y[wy] <- mice.impute.norm.bs(incomplete.data$y, !is.na(incomplete.data$y), cbind(imp.data$x, imp.data$z), 
                                          50, c(0, 0, 0), matrix(c(60, 0, 0, 0, 3600, 0, 0, 0, 3600), 3, 3),
                                          1, 30)
    
  }
  beta.tar <- beta.1 - beta.2
  for (i in 1 : 50) {
    batch.mean[i] <- mean(beta.tar[seq(10+i, 960+i, 50)])
  }
  beta.tar.se <- sd(batch.mean) / sqrt(50)
  beta.tar.mmean[j] <- beta.tar.mean <- mean(beta.tar)
  count[j] <- (beta.tar.mean - 1.96 * beta.tar.se) < 0 & 0 < (beta.tar.mean + 1.96 * beta.tar.se)
  setTxtProgressBar(pb, j)
}
close(pb)
sum(count)/500
