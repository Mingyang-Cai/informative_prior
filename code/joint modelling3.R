#test the joint modelling (three columns)
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
nsim <- 500
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


impute.multiple <- function(incomplete.data, nit){
  multiple.imputed.dataset <- list()
  for (i in 1 : 5) {
    multiple.imputed.dataset[[i]] <- impute.single(incomplete.data, nit) 
  }
  return(multiple.imputed.dataset)
}


OUT <- list()
pb <- txtProgressBar(min = 0, max = nsim, style = 3)
for (i in 1 : nsim) {
  complete.data <- data.frame(mvrnorm(n = 200, mu = c(1, 4, 9), Sigma = matrix(c(4, 2, 2, 2, 4, 2, 2, 2, 9), nrow = 3)))
  colnames(complete.data) <- c("x", "y", "z")
  mypattern <- rbind(c(1, 0, 0), c(0, 1, 0), c(0, 0, 1))
  incomplete.data <- ampute(complete.data, prop = 0.5, patterns = mypattern, mech = "MCAR")$amp
  OUT[[i]] <- impute.multiple(incomplete.data, 100)
  
  setTxtProgressBar(pb, i)
}
close(pb)  

# Evaluate mean y
evaluate.sims <- function(sims, truth = 4){
  POOL <- list()
  pb <- txtProgressBar(min = 0, max = nsim, style = 3)
  for (i in 1:length(sims)){
    #Extract means and variances
    Q            <- c(mean(sims[[i]][[1]]$y), mean(sims[[i]][[2]]$y),
                      mean(sims[[i]][[3]]$y), mean(sims[[i]][[4]]$y),
                      mean(sims[[i]][[5]]$y))
    
    U            <- c(var(sims[[i]][[1]]$y), var(sims[[i]][[2]]$y),
                      var(sims[[i]][[3]]$y), var(sims[[i]][[4]]$y), 
                      var(sims[[i]][[5]]$y)) / 200
    #Pool the regular way
    pool         <- mice::pool.scalar(Q, U, n = 1000) # A really large number
    pool$lower   <- pool$qbar - qt(0.975, pool$df) * sqrt(pool$t)
    pool$upper   <- pool$qbar + qt(0.975, pool$df) * sqrt(pool$t)
    pool$coverage <- pool$lower < mean(truth) & mean(truth) < pool$upper
    POOL[[i]]     <- unlist(pool)
    setTxtProgressBar(pb, i)
  }
  return(POOL)
  close(pb)
}
EVAL <- evaluate.sims(OUT)
# Summarize
AVG.EVAL <- Reduce("+", EVAL) / length(EVAL)
round(AVG.EVAL[c(12, 13, 19, 20, 21)], 2)
  